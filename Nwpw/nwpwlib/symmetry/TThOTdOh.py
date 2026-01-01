import itertools
import numpy as np

# ---------- helpers ----------

def det_int(M):
    return int(round(np.linalg.det(M)))

def flattenF(M):
    return M.flatten(order="F")

def print_group(name, mats):
    print(name)
    print(len(mats))
    for M in mats:
        flat = flattenF(M)
        print(" ".join(f"{int(x):2d}" for x in flat))
    print()

# ---------- build O (proper octahedral group) ----------

O = []
for perm in itertools.permutations([0,1,2]):
    for signs in itertools.product([1,-1], repeat=3):
        M = np.zeros((3,3), dtype=int)
        for i,j in enumerate(perm):
            M[i,j] = signs[i]
        if det_int(M) == 1:
            O.append(M)

assert len(O) == 24

# ---------- tetrahedron inside cube ----------

tet = [
    np.array([ 1, 1, 1]),
    np.array([ 1,-1,-1]),
    np.array([-1, 1,-1]),
    np.array([-1,-1, 1]),
]
tetset = {tuple(v) for v in tet}

def preserves_tet(M):
    return {tuple(M @ v) for v in tet} == tetset

# ---------- T and Th ----------

T = [M for M in O if preserves_tet(M)]
assert len(T) == 12

Th = T + [-M for M in T]
assert len(Th) == 24

# ---------- inversion ----------

I = -np.eye(3, dtype=int)

# ---------- Oh ----------

Oh = O + [-M for M in O]
assert len(Oh) == 48

# ---------- Td ----------
# Td = all operations mapping tetrahedron to itself,
# including det = -1 (but excluding inversion-only elements)

Td = []
for M in itertools.chain(O, [-M for M in O]):
    if preserves_tet(M) and not np.array_equal(M, I):
        Td.append(M)

# remove duplicates
Td_unique = []
for M in Td:
    if not any(np.array_equal(M,N) for N in Td_unique):
        Td_unique.append(M)

Td = Td_unique
assert len(Td) == 24

# ---------- print everything ----------

print_group("T",  T)
print_group("Th", Th)
print_group("O",  O)
print_group("Td", Td)
print_group("Oh", Oh)

