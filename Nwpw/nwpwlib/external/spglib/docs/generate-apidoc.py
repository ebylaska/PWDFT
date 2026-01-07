# /// script
# dependencies = [
#   "sphinx",
# ]
# ///
"""Helper for running sphinx-apidoc."""

from __future__ import annotations

from pathlib import Path

from sphinx.ext import apidoc

ROOT_DIR = Path(__file__).parent.parent
PYTHON_API_DIR = ROOT_DIR / "python/spglib"
AOUTODOC_DIR = ROOT_DIR / "docs/api/autodoc"

apidoc.main(
    [
        "--force",
        "--remove-old",
        "--no-toc",
        "--module-first",
        "--separate",
        "--private",  # So that we document spglib._internal also
        # Output
        f"--output-dir={AOUTODOC_DIR}",
        # Module path
        str(PYTHON_API_DIR),
        # Exclude patterns
        str(PYTHON_API_DIR / "_compat"),
        str(PYTHON_API_DIR / "_version*"),
    ]
)
