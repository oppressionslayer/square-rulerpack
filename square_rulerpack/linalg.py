"""GF(2) linear algebra for the relation matrix: singleton pruning (as the Rust
``prune_singletons``) and dense bit-packed Gaussian elimination in numpy.

Each relation is a row; column c is set when that relation has an odd exponent at c. A
dependency is a set of rows whose XOR is zero: multiplying those relations gives a square on
the right-hand side.
"""
from __future__ import annotations

from collections import deque

import numpy as np

__all__ = ["prune_singletons", "dependencies", "dense_matrix"]


def prune_singletons(odd: list[list[int]], cols: int) -> list[int]:
    """Indices of the rows that survive removing, repeatedly, every row that is the only one
    with a set bit in some column (such a row can be in no dependency)."""
    count = [0] * cols
    owner = [0] * cols
    for i, r in enumerate(odd):
        for c in r:
            count[c] += 1
            owner[c] ^= i
    alive = [True] * len(odd)
    queue = deque(c for c in range(cols) if count[c] == 1)
    while queue:
        c = queue.popleft()
        if count[c] != 1:
            continue
        i = owner[c]
        if not alive[i]:
            continue
        alive[i] = False
        for d in odd[i]:
            count[d] -= 1
            owner[d] ^= i
            if count[d] == 1:
                queue.append(d)
    return [i for i in range(len(odd)) if alive[i]]


def dense_matrix(rows: list[list[int]], cols: int) -> np.ndarray:
    """The 0/1 matrix (rows × cols) as uint8, for looking at small cases."""
    m = np.zeros((len(rows), cols), dtype=np.uint8)
    for i, r in enumerate(rows):
        m[i, r] = 1
    return m


def dependencies(rows: list[list[int]], cols: int, max_deps: int = 64) -> list[list[int]]:
    """Up to ``max_deps`` dependencies (lists of row indices) of the sparse 0/1 matrix given
    as per-row column lists. Dense elimination over 64-bit words: the matrix is packed as
    [columns | identity]; eliminating the column part column by column leaves, in every row
    that never became a pivot, the identity part naming the rows that XOR to zero."""
    n = len(rows)
    if n == 0:
        return []
    W = -(-cols // 64)
    Wi = -(-n // 64)
    M = np.zeros((n, W + Wi), dtype=np.uint64)
    one = np.uint64(1)
    for i, r in enumerate(rows):
        for c in r:
            M[i, c >> 6] |= one << np.uint64(c & 63)
        M[i, W + (i >> 6)] |= one << np.uint64(i & 63)
    used = np.zeros(n, dtype=bool)
    # columns that appear at all, in order (columns nobody has need no pass)
    present = np.zeros(cols, dtype=bool)
    for r in rows:
        present[r] = True
    for c in np.flatnonzero(present):
        w, bit = int(c) >> 6, np.uint64(int(c) & 63)
        has = ((M[:, w] >> bit) & one).astype(bool)
        with_bit = np.flatnonzero(has)
        if len(with_bit) == 0:
            continue
        cand = with_bit[~used[with_bit]]
        if len(cand) == 0:
            continue
        piv = cand[0]
        used[piv] = True
        others = with_bit[with_bit != piv]
        if len(others):
            M[others] ^= M[piv]
    deps = []
    for i in np.flatnonzero(~used):
        assert not M[i, :W].any(), "a non-pivot row kept a column bit"
        words = M[i, W:]
        members = []
        for wi in np.flatnonzero(words):
            word = int(words[wi])
            while word:
                low = word & -word
                members.append(int(wi) * 64 + low.bit_length() - 1)
                word ^= low
        deps.append(members)
        if len(deps) >= max_deps:
            break
    return deps
