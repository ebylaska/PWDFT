#!/usr/bin/env python3
"""
pwdft-python.py

Outer workflow driver for PWDFT (.nwx) with per-task lowering into atomic units.

Core semantics (as we agreed):
- The .nwx file is a declarative workflow.
- Geometry + initial nwpw are latent until the first native PWDFT task.
- The first task MUST be engine=pwdft.
- Subsequent nwpw blocks are overlays consumed by the *next* native PWDFT task only.
- Python-engine tasks operate on (restart + global JSON), and may call PWDFT only via
  atomic single-task invocations (never by re-running the whole workflow).
- Only PWDFT writes the GLOBAL truth JSON: <calc>.json
- Python writes SIDE truths: aux/<op>-python[-N].json

Usage:
  python3 pwdft-python.py benzene.nwx
  python3 pwdft-python.py benzene.nwx --dry-run --keep-units

Environment:
  PWDFT_EXE: optional explicit path/name to pwdft binary.
            If unset, driver searches: pwdft, pwdft.x
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# -----------------------------------------------------------------------------
# Logging / small utils
# -----------------------------------------------------------------------------

def info(msg: str) -> None:
    print(f"[pwdft-python] {msg}", file=sys.stdout)

def warn(msg: str) -> None:
    print(f"[pwdft-python] WARNING: {msg}", file=sys.stderr)

def die(msg: str, code: int = 2) -> None:
    print(f"[pwdft-python] ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)

def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

def atomic_write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write("\n")
    tmp.replace(path)

def read_json(path: Path) -> Dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        die(f"Missing JSON file: {path}")
    except json.JSONDecodeError as e:
        die(f"Invalid JSON in {path}: {e}")

def find_pwdft_exe() -> str:
    """
    Resolve PWDFT executable.

    Search order:
      1) $PWDFT_EXE (explicit override)
      2) pwdft
      3) pwdft.x
    """
    env = os.environ.get("PWDFT_EXE", "").strip()
    if env:
        p = Path(env)
        if p.exists():
            return str(p.resolve())
        found = shutil.which(env)
        if found:
            return found
        die(f"PWDFT_EXE is set but not found/executable: {env}")

    for name in ("pwdft", "pwdft.x"):
        found = shutil.which(name)
        if found:
            return found

    die("Cannot find PWDFT executable (tried: $PWDFT_EXE, pwdft, pwdft.x)")
    return ""  # unreachable


# -----------------------------------------------------------------------------
# IR (intermediate representation)
# -----------------------------------------------------------------------------

@dataclass
class Block:
    kind: str            # "geometry" or "nwpw"
    lines: List[str]     # includes opening keyword and trailing 'end'

@dataclass
class TaskIR:
    module: str          # e.g., "pspw"
    operation: str       # e.g., "energy"
    engine: str          # "pwdft" or "python"
    options: Dict[str, str] = field(default_factory=dict)
    raw_line: str = ""

@dataclass
class ProgramIR:
    calc_name: str
    start_or_restart: str            # "start" or "restart" (as declared in the .nwx)
    header_lines: List[str]          # Title/echo/memory/dirs/comments/blank lines, BUT NO geometry/nwpw/task lines
    geometry: Optional[Block]
    initial_nwpw: Optional[Block]
    tasks: List[TaskIR]
    overlays_before: List[List[Block]]  # overlays_before[i] are nwpw blocks that appear between task i-1 and task i


@dataclass
class ExecUnit:
    kind: str               # "pwdft" or "python"
    calc_name: str
    module: str
    operation: str
    engine: str             # "pwdft" or "python"
    unit_file: Optional[Path] = None   # for pwdft units
    task_json: Optional[Dict[str, Any]] = None  # for python units
    note: str = ""


# -----------------------------------------------------------------------------
# Parser (line-based, preserves block text exactly)
# -----------------------------------------------------------------------------

_START_RE = re.compile(r"^\s*(start|restart)\s+(\S+)\s*$", re.IGNORECASE)
_TASK_RE = re.compile(r"^\s*task\s+(\S+)\s+(\S+)(.*)$", re.IGNORECASE)
_ENGINE_KV_RE = re.compile(r"(?:^|\s)engine\s*=\s*([A-Za-z0-9_\-\.]+)", re.IGNORECASE)
_BLOCK_START_RE = re.compile(r"^\s*(geometry|nwpw)\s*$", re.IGNORECASE)
_BLOCK_END_RE = re.compile(r"^\s*end\s*$", re.IGNORECASE)


def _indent(s: str) -> int:
    return len(s) - len(s.lstrip(" "))

def _is_blank_or_comment(s: str) -> bool:
    t = s.strip()
    return (t == "") or t.startswith("#")

_SINGLE_TOKEN_RE = re.compile(r"^\s*([A-Za-z0-9_\-]+)\s*$")

def read_block(lines: List[str], i: int) -> Tuple[Block, int]:
    """
    Read a top-level geometry or nwpw block starting at lines[i].

    geometry:
      - ends at first 'end'

    nwpw:
      - may contain nested sub-blocks with their own 'end'
      - closes only when the *outer* nwpw 'end' is reached
    """
    m = _BLOCK_START_RE.match(lines[i])
    if not m:
        die("Internal error: read_block called on non-block start")

    kind = m.group(1).lower()
    blk_lines = [lines[i]]
    i += 1

    # geometry has no nested blocks in your usage
    if kind == "geometry":
        while i < len(lines):
            blk_lines.append(lines[i])
            if _BLOCK_END_RE.match(lines[i]):
                return Block(kind=kind, lines=blk_lines), i + 1
            i += 1
        die("Unterminated geometry block")

    # nwpw block: track nested sub-block depth
    depth = 0
    while i < len(lines):
        line = lines[i]
        blk_lines.append(line)

        if _BLOCK_END_RE.match(line):
            if depth > 0:
                depth -= 1
            else:
                # closes the nwpw block
                return Block(kind=kind, lines=blk_lines), i + 1
            i += 1
            continue

        # detect nested sub-blocks by indentation
        m1 = _SINGLE_TOKEN_RE.match(line)
        if m1:
            tok = m1.group(1).lower()
            if tok not in ("nwpw", "end", "task", "start", "restart", "geometry"):
                curr_indent = _indent(line)
                j = i + 1
                while j < len(lines) and _is_blank_or_comment(lines[j]):
                    j += 1
                if j < len(lines) and _indent(lines[j]) > curr_indent:
                    depth += 1

        i += 1

    die("Unterminated nwpw block")



def parse_keyvals(tail: str) -> Dict[str, str]:
    """
    Parse simple whitespace-separated key=value tokens from the tail of a task line.
    (We keep it intentionally simple.)
    """
    out: Dict[str, str] = {}
    for tok in tail.strip().split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k.strip()] = v.strip()
    return out

def parse_nwx(path: Path) -> ProgramIR:
    text_lines = path.read_text(encoding="utf-8").splitlines(True)  # keep newlines

    calc_name: Optional[str] = None
    start_or_restart: Optional[str] = None

    header_lines: List[str] = []       # will exclude geometry/nwpw/task lines
    geometry: Optional[Block] = None
    initial_nwpw: Optional[Block] = None

    tasks: List[TaskIR] = []
    overlays_before: List[List[Block]] = []
    pending_overlays: List[Block] = []

    seen_first_task = False

    i = 0
    while i < len(text_lines):
        line = text_lines[i]

        # start/restart
        m_sr = _START_RE.match(line)
        if m_sr:
            if start_or_restart is not None:
                die("Multiple start/restart lines are not allowed.")
            start_or_restart = m_sr.group(1).lower()
            calc_name = m_sr.group(2)
            header_lines.append(line)
            i += 1
            continue

        # geometry / nwpw block
        m_bs = _BLOCK_START_RE.match(line)
        if m_bs:
            blk, i = read_block(text_lines, i)

            if not seen_first_task:
                if blk.kind == "geometry":
                    if geometry is not None:
                        die("Multiple geometry blocks before first task are not supported.")
                    geometry = blk
                elif blk.kind == "nwpw":
                    if initial_nwpw is None:
                        initial_nwpw = blk
                    else:
                        pending_overlays.append(blk)
            else:
                if blk.kind == "nwpw":
                    pending_overlays.append(blk)
                else:
                    die("geometry block after tasks is not supported (restart carries geometry).")
            continue


        # task line
        m_task = _TASK_RE.match(line)
        if m_task:
            module = m_task.group(1)
            operation = m_task.group(2)
            tail = m_task.group(3) or ""
            opts = parse_keyvals(tail)

            eng = "pwdft"
            m_eng = _ENGINE_KV_RE.search(tail)
            if m_eng:
                eng = m_eng.group(1).lower()
            if eng not in ("pwdft", "python"):
                die(f"Unknown engine '{eng}' in task line: {line.strip()}")

            tasks.append(TaskIR(
                module=module,
                operation=operation,
                engine=eng,
                options=opts,
                raw_line=line.rstrip("\n"),
            ))

            overlays_before.append(pending_overlays)
            pending_overlays = []
            seen_first_task = True
            i += 1
            continue

        # otherwise: header/comment/blank
        header_lines.append(line)
        i += 1

    if start_or_restart is None or calc_name is None:
        die(f"{path}: missing required 'start <name>' or 'restart <name>'.")

    if not tasks:
        die(f"{path}: no tasks found.")

    if pending_overlays:
        die("Dangling nwpw block at end of file: overlay must precede a PWDFT task.")

    # Add the pre-task extra-nwpw overlays to the first task's overlay list (already stored in overlays_before[0])
    # because they were pending at the moment we hit first task, they are in overlays_before[0].

    return ProgramIR(
        calc_name=calc_name,
        start_or_restart=start_or_restart,
        header_lines=header_lines,
        geometry=geometry,
        initial_nwpw=initial_nwpw,
        tasks=tasks,
        overlays_before=overlays_before,
    )


# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

def validate_ir(ir: ProgramIR) -> None:
    if len(ir.tasks) != len(ir.overlays_before):
        die("Internal error: overlays_before length mismatch.")

    if ir.tasks[0].engine != "pwdft":
        die("First task must be engine=pwdft (native PWDFT).")

    if ir.start_or_restart == "start" and ir.geometry is None:
        die("For 'start' workflows, geometry must be provided before the first PWDFT task.")

    # Minimal: only pspw supported for now
    for t in ir.tasks:
        if t.module.lower() != "pspw":
            die(f"Only module 'pspw' is supported by this driver (got '{t.module}').")


# -----------------------------------------------------------------------------
# Lowering
# -----------------------------------------------------------------------------

def unit_dir(workdir: Path) -> Path:
    d = workdir / ".pwdft_units"
    d.mkdir(parents=True, exist_ok=True)
    return d

def emit_pwdft_unit_file(
    *,
    ir: ProgramIR,
    workdir: Path,
    unit_index: int,
    sr: str,                           # "start" or "restart"
    task: TaskIR,
    include_geometry: bool,
    include_initial_nwpw: bool,
    overlay_nwpw: List[Block],
) -> Path:
    """
    Emit a minimal, syntactically clean .nw file for a single PWDFT task.

    Canonical ordering:
      header (with sr line rewritten)
      geometry (optional)
      initial nwpw (optional)
      overlay nwpw blocks
      single task line
    """
    udir = unit_dir(workdir)
    f = udir / f"{ir.calc_name}.unit{unit_index:03d}.{task.operation}.nw"

    out: List[str] = []

    # Header, but rewrite start/restart line deterministically
    for ln in ir.header_lines:
        if _START_RE.match(ln):
            out.append(f"{sr} {ir.calc_name}\n")
        else:
            out.append(ln)

    # Canonical blocks
    if include_geometry and ir.geometry is not None:
        out.append("\n")
        out.extend(ir.geometry.lines)

    if include_initial_nwpw and ir.initial_nwpw is not None:
        out.append("\n")
        out.extend(ir.initial_nwpw.lines)

    for blk in overlay_nwpw:
        if blk.kind != "nwpw":
            die("Internal error: non-nwpw overlay encountered.")
        out.append("\n")
        out.extend(blk.lines)

    # Single task line (force no engine=)
    out.append("\n")
    out.append(f"task {task.module} {task.operation}\n")

    f.write_text("".join(out), encoding="utf-8")
    return f

def lower_to_units(ir: ProgramIR, workdir: Path) -> List[ExecUnit]:
    units: List[ExecUnit] = []
    instantiated = False
    pwdft_unit_index = 0
    python_seq = 0

    for i, t in enumerate(ir.tasks):
        overlays = ir.overlays_before[i] or []

        if t.engine == "pwdft":
            # Determine sr and inclusion of geometry/nwpw:
            if not instantiated:
                sr = ir.start_or_restart  # may be "start" or "restart"
                include_geom = True
                include_nwpw = True
            else:
                sr = "restart"
                include_geom = False
                include_nwpw = False

            unit_file = emit_pwdft_unit_file(
                ir=ir,
                workdir=workdir,
                unit_index=pwdft_unit_index,
                sr=sr,
                task=t,
                include_geometry=include_geom,
                include_initial_nwpw=include_nwpw,
                overlay_nwpw=overlays,
            )

            units.append(ExecUnit(
                kind="pwdft",
                calc_name=ir.calc_name,
                module=t.module,
                operation=t.operation,
                engine="pwdft",
                unit_file=unit_file,
                note=("instantiation" if not instantiated else "continuation"),
            ))
            instantiated = True
            pwdft_unit_index += 1

        else:
            if not instantiated:
                die(f"Python task '{t.operation}' appears before instantiation.")
            python_seq += 1
            units.append(ExecUnit(
                kind="python",
                calc_name=ir.calc_name,
                module=t.module,
                operation=t.operation,
                engine="python",
                task_json={
                    "restart": ir.calc_name,
                    "module": t.module,
                    "operation": t.operation,
                    "engine": "python",
                    "global_json": f"{ir.calc_name}.json",
                    "options": t.options,
                    "workdir": str(workdir),
                    "invocation": {"mode": "driver"},
                    "sequence": python_seq,
                },
                note="python-engine",
            ))

    return units


# -----------------------------------------------------------------------------
# PWDFT execution
# -----------------------------------------------------------------------------

def run_pwdft(unit_file: Path, workdir: Path) -> None:
    exe = find_pwdft_exe()
    cmd = [exe, str(unit_file)]
    info(f"PWDFT: {' '.join(cmd)}")

    proc = subprocess.Popen(
        cmd,
        cwd=str(workdir),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,          # line-buffered
        universal_newlines=True,
    )

    assert proc.stdout is not None
    assert proc.stderr is not None

    # Stream stdout live
    for line in proc.stdout:
        print(line, end="")

    # Drain stderr at the end (or also stream if you prefer)
    err = proc.stderr.read()
    if err.strip():
        warn("PWDFT stderr:\n" + err)

    rc = proc.wait()
    if rc != 0:
        die(f"PWDFT unit failed (rc={rc}) for {unit_file.name}")



# -----------------------------------------------------------------------------
# Python engine implementations
# -----------------------------------------------------------------------------

def aux_side_truth_path(workdir: Path, op: str, seq: int) -> Path:
    d = workdir / "aux"
    d.mkdir(parents=True, exist_ok=True)
    suffix = f"-{seq}" if seq > 1 else ""
    return d / f"{op}-python{suffix}.json"



def python_call_pwdft_atomic(
    workdir: Path,
    calc: str,
    module: str,
    operation: str,
    nwpw_block: List[str],
):
    udir = unit_dir(workdir)
    f = udir / f"{calc}.pycall.{module}.{operation}.nw"

    lines = [
        f"restart {calc}\n\n",
        "permanent_dir ./perm\n",
        "scratch_dir   ./perm\n\n",
    ]

    lines.extend(nwpw_block)
    lines.append("\n")
    lines.append(f"task {module} {operation}\n")

    f.write_text("".join(lines), encoding="utf-8")
    run_pwdft(f, workdir)




def python_engine_optimize(task_json: Dict[str, Any]) -> Dict[str, Any]:
    """
    Minimal optimize engine:
    - Reads global truth (if exists)
    - Delegates to native PWDFT optimize via atomic call (default)
    - Reads global truth after
    - Writes side truth
    """
    workdir = Path(task_json["workdir"]).resolve()
    calc = task_json["restart"]
    module = task_json.get("module", "pspw")
    global_json = workdir / task_json.get("global_json", f"{calc}.json")
    options = task_json.get("options", {}) or {}
    seq = int(task_json.get("sequence", 1))

    before = read_json(global_json) if global_json.exists() else None

    delegate = str(options.get("delegate", "pwdft")).lower()
    did_delegate = False
    delegate_error = None

    info(f"PYTHON optimize: delegating to PWDFT optimize")

    if delegate in ("pwdft", "1", "true", "yes", "on"):
        try:
            python_call_pwdft_atomic(workdir, calc, module, "optimize")
            did_delegate = True
        except Exception as e:
            delegate_error = str(e)

    after = read_json(global_json) if global_json.exists() else None

    return {
        "operation": "optimize",
        "engine": "python",
        "restart": calc,
        "module": module,
        "timestamp": now_iso(),
        "sequence": seq,
        "delegate": {"mode": delegate, "did_delegate": did_delegate, "error": delegate_error},
        "before": {
            "exists": before is not None,
            "status": before.get("status") if before else None,
            "energy": before.get("energy") if before else None,
            "geometry_hash": (before.get("geometry", {}).get("hash") if before else None),
        },
        "after": {
            "exists": after is not None,
            "status": after.get("status") if after else None,
            "energy": after.get("energy") if after else None,
            "geometry_hash": (after.get("geometry", {}).get("hash") if after else None),
        },
        "notes": [
            "Minimal python optimize engine.",
            "Currently delegates to PWDFT optimize using an atomic restart+task call.",
        ],
    }

def python_engine_freq(task_json: Dict[str, Any]) -> Dict[str, Any]:
    """
    Placeholder freq engine (side truth only).
    Later: finite differences calling PWDFT energy/forces atomically.
    """
    workdir = Path(task_json["workdir"]).resolve()
    calc = task_json["restart"]
    module = task_json.get("module", "pspw")
    global_json = workdir / task_json.get("global_json", f"{calc}.json")
    options = task_json.get("options", {}) or {}
    seq = int(task_json.get("sequence", 1))

    g = read_json(global_json) if global_json.exists() else None

    return {
        "operation": "freq",
        "engine": "python",
        "restart": calc,
        "module": module,
        "timestamp": now_iso(),
        "sequence": seq,
        "options": options,
        "global": {
            "exists": g is not None,
            "status": g.get("status") if g else None,
            "geometry_hash": (g.get("geometry", {}).get("hash") if g else None),
        },
        "notes": [
            "freq engine not implemented yet.",
            "Plan: finite-difference displacements and PWDFT force evaluations.",
        ],
    }

def run_python_engine(task_json: Dict[str, Any]) -> Path:

    op = str(task_json["operation"]).lower()
    workdir = Path(task_json["workdir"]).resolve()
    seq = int(task_json.get("sequence", 1))

    info(f"PYTHON: begin {op}")

    if op == "optimize":
        side = python_engine_optimize(task_json)
    elif op == "freq":
        side = python_engine_freq(task_json)
    else:
        die(f"Python engine does not implement operation '{op}'")

    out = aux_side_truth_path(workdir, op, seq)
    atomic_write_json(out, side)
    info(f"Python side truth written: {out}")

    info(f"PYTHON: end   {op}")

    return out


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

def cleanup_units(workdir: Path, calc_name: str) -> None:
    d = unit_dir(workdir)
    for p in d.glob(f"{calc_name}.unit*.nw"):
        try:
            p.unlink()
        except Exception:
            pass
    for p in d.glob(f"{calc_name}.pycall.*.nw"):
        try:
            p.unlink()
        except Exception:
            pass

def run_workflow(nwx_path: Path, dry_run: bool, keep_units_flag: bool) -> int:
    workdir = nwx_path.parent.resolve()

    ir = parse_nwx(nwx_path)
    validate_ir(ir)
    units = lower_to_units(ir, workdir)

    info(f"Parsed calc='{ir.calc_name}' mode='{ir.start_or_restart}' tasks={len(ir.tasks)} units={len(units)}")
    for k, u in enumerate(units):
        if u.kind == "pwdft":
            info(f"unit[{k}] PWDFT  op={u.operation:10s} file={u.unit_file.name if u.unit_file else None} note={u.note}")
        else:
            info(f"unit[{k}] PYTHON op={u.operation:10s} note={u.note}")

    if dry_run:
        info("Dry run: not executing.")
        return 0

    for u in units:
        if u.kind == "pwdft":
            assert u.unit_file is not None
            run_pwdft(u.unit_file, workdir)
        else:
            assert u.task_json is not None
            run_python_engine(u.task_json)

    if not keep_units_flag:
        cleanup_units(workdir, ir.calc_name)

    return 0

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("nwx", help="Input .nwx workflow file")
    ap.add_argument("--dry-run", action="store_true", help="Parse/lower only; do not execute")
    ap.add_argument("--keep-units", action="store_true", help="Keep generated unit files under .pwdft_units/")
    args = ap.parse_args()

    nwx_path = Path(args.nwx).resolve()
    if not nwx_path.exists():
        die(f"Missing input file: {nwx_path}")

    return run_workflow(nwx_path, dry_run=args.dry_run, keep_units_flag=args.keep_units)

if __name__ == "__main__":
    raise SystemExit(main())

