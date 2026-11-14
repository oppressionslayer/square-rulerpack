# square-rulerpack in Python

The self-initialising quadratic sieve with double large primes of
[square-rulerpack 0.2.0](https://crates.io/crates/square-rulerpack) (Rust, by Lars Rocha),
ported stage by stage to Python so that every step can be run by hand in IPython and looked
at as a DataFrame or a numpy array. Same parameter tables, same Knuth–Schroeppel multiplier,
same family selection (family *i* picks the same primes for *a* as the Rust engine's family *i*),
same Gray-code polynomial switching, same Q1/Q2 forms, same exact and plain thresholds, same
large-prime cycle finder and square root. The sieve kernel is numba when it is installed and
numpy slices otherwise; relation collection runs across processes.

```sh
pip install square-rulerpack            # from PyPI; add '[fast]' for the numba kernel

or

pip install square-rulerpack[native]

or

pip install -e .                        # or from this checkout
python -m square_rulerpack 1009732533765211 --workers 8
```

```
1009732533765211 = 11344301 × 89007911
2 factors, product checks: true, all probable prime: true, 0.34 s
```

## Factor a number

```python
import square_rulerpack as sq
sq.factor_completely(1009732533765211)                 # [11344301, 89007911]
sq.factor_completely(n, workers=8, verbose=True)       # progress on stderr, like the Rust CLI
sq.factor(n, workers=4)                                # one nontrivial factor or None
```

`Settings` carries what the Rust engine reads from the environment:

| Settings field | Rust variable | default |
|---|---|---|
| `exact`   | `SIQS_EXACT` | exact threshold from 215 bits (about 65 digits), plain below |
| `small`   | `SIQS_SMALL` | 512 unsieved primes under the exact threshold, 30 under the plain |
| `q2`      | `SIQS_Q2`    | Q2 polynomials when k·n ≡ 1 (mod 8) |
| `dlp`     | `SIQS_DLP`   | two large primes |
| `close_a` | `SIQS_CLOSE_A` | repair *a* toward √(2kn)/M |
| `ecap`, `margin` | `SIQS_ECAP`, `SIQS_MARGIN` | 32, 1 |
| `multiplier` | `SIQS_MULT` | choose by score |

```python
sq.factor_completely(n, settings=sq.Settings(exact=False, q2=False))
```

## Walk through a factorisation by hand

`Walkthrough` runs the stages one at a time. Each method prints a one-line summary on stderr
and returns the stage's data; the object keeps everything (`w.k`, `w.kn`, `w.base`, `w.fam`,
`w.fs`, `w.full`, `w.cycles`, `w.deps` …) for further digging. `WALKTHROUGH.md` is a full
session on a 40-digit number; `examples_walkthrough.py` reproduces it.

```python
import square_rulerpack as sq
w = sq.Walkthrough(5857505067604849028619727330686829482307, workers=4)  # python only implementation not the fatest, but it will factor
#  sq.factor_completely(a, workers=8, verbose=True) is way faster, it uses the rust engine and will factor a 100 digit number in 12 minutes dependent on your machine
w.summary()          # digits, bits, trial division, square and primality checks, parameters
w.multiplier()       # Knuth–Schroeppel scores for every candidate k; the chosen one
fb = w.factor_base() # DataFrame: p, sqrt(kn) mod p, both roots, log2 p, byte log, sieved?
w.family()           # a = q1 q2 … q_s, the B_l, how many polynomials (2^(s-1))
w.polynomials()      # every b of the family in Gray-code order, lin and c of g(x)
w.roots(0)           # the two roots of g mod each prime, their sieve starts, g(root) mod p == 0
s = w.sieve(0)       # s.acc: the int16 sieve array over x in [-M, M]; s.hits: the survivors
s.table(20)          # positions by sieve sum: x, g(x), log2|g|, sum, threshold, survives
w.trial_divide(0)    # each survivor: factor-base part, cofactor, kind (full / one-large /
                     #   two-large / drop: why), y, and the check y² ≡ product (mod kn)
w.collect()          # the rest with worker processes until enough relations
w.relations(); w.partials(); w.cycle_summary(); w.relation(0)
w.matrix()           # 0/1 DataFrame (odd exponents) for small cases, sparse rows otherwise
w.dependencies()     # row sets that XOR to zero
w.square_root(0)     # X, Y, X² ≡ Y² (mod kn), gcd(X − Y, n), gcd(X + Y, n), the factor
w.factors()          # try the dependencies in order
```

Or `w.run()` for everything at once.

What the stages are, in the engine's own words (module docstrings say more):

1. **multiplier** — pick k so that k·n is a square modulo many small primes.
2. **factor base** — the first `base` primes with (k·n / p) = 1, each with √(k·n) mod p.
3. **family** — a = q₁⋯q_s from the primes just above 1000 with a ≈ √(2kn)/M; B_l with
   B_l ≡ √kn (mod q_l) and ≡ 0 modulo the other q's; every b = Σ ±B_l gives b² ≡ kn (mod a).
4. **polynomial** — g(x) = a·x² + 2b·x + (b² − kn)/a (Q1) or a·x² + b·x + (b² − kn)/(4a)
   (Q2, for kn ≡ 1 mod 8); (a·x + b)² ≡ a·g(x) (mod kn), so a smooth g(x) is a relation.
   Consecutive b's differ by one B_l (Gray code), so every root moves by a precomputed delta.
5. **sieve** — each prime adds its byte log2 p at x ≡ root (mod p) along the interval; a
   position whose sum reaches the threshold is a survivor. The exact threshold compares each
   survivor with its own log2|g(x)| and counts the unsieved small primes exactly.
6. **trial division** — p | g(x) ⇔ x ≡ root (mod p), so the dividing primes are read off the
   roots; the cofactor is 1 (full), one prime below the bound, or two (split by Pollard–Brent).
7. **cycles** — partials are edges between large primes; a cycle multiplies to a full relation.
8. **matrix** — odd exponents over GF(2); singleton pruning; dense elimination for dependencies.
9. **square root** — X = ∏ y, Y = ∏ p^(Σe/2); gcd(X − Y, n).

## Small exact lessons (`square_rulerpack.lessons`)

For a class that wants to see the whole idea on numbers that fit on a blackboard, before the
real sieve. Everything is plain Python integers over a handful of samples; no threshold, no
large primes, no engine. Carried over from Lars Rocha's `square_research` package.

```python
import square_rulerpack as sq

p = sq.construct_polynomial(221, [5, 7])      # a = 35; every CRT term of b; (c, 2b, a) = (1, 32, 35)
p["crt_terms"]                                # per prime: roots of 221 mod p, chosen root, a/p, inverse, term

w = sq.siqs_window(p["k"], p["a"], p["b"], start=-2, width=5, seeds=[2, 5, 7, 11])
w.dataframe()                                 # x, root = a·x + b, Q(x), root² − 221 = a·Q(x), factors, smooth
sq.sieve_walk(w)                              # the sieve by hand: roots of Q mod p and the positions they predict

lesson = sq.factor_walkthrough(w)             # FactorWalkthrough(n=221, relations=3, rank=2, nullity=1, split=(13, 17))
lesson.matrix()                               # parity matrix: rows = relations, columns = (−1, 2, 5, 7, 11)
lesson.dataframe("elimination")               # the Gaussian elimination, row by row
lesson.check_dependency([1, 2])               # 4900 = 70², x = 138, gcd(138 − 70, 221) = 17, gcd(138 + 70, 221) = 13

sq.relation_study()                           # all four CRT root choices of a = 35, left/centre/right windows
sq.factor_walkthrough()                       # the same N = 221 lesson with no arguments
```

The bridge to the real engine: `engine_polynomials` asks this package's own selector (the same
families as the Rust engine) for a few polynomials and turns them into such windows.

```python
cap = sq.engine_polynomials(1009732533765211, base=64, half=128, families=8, per_family=4)
cap.dataframe()                               # family, index, a, b, c, a_primes, coefficients, identity_valid
windows = [cap.window(i, start=-128, width=257) for i in range(16)]
sq.factor_walkthrough(windows, original_n=cap.n).split      # (11344301, 89007911)
```

Also there: `square_window(n)` for x² − n, `polynomial_roots`, `sqrt_mod` (all roots), `prime_seeds`,
`survivor_mask`, `dependency`, and `LessonWindow.export` for JSON.

## Fast path: the Rust engine as an extension module

The teaching code above is Python on purpose. For actually factoring, the same engine is also
available compiled, the way gmpy2 wraps GMP: `native/` in this repository is a PyO3 crate that links
the published `square-rulerpack 0.2.0` and builds into the extension module `square_rulerpack_native`
(PyPI package `square-rulerpack-native`, one abi3 wheel per platform for every Python ≥ 3.10).

```sh
pip install square-rulerpack[native]
```

Once it is installed, `factor` and `factor_completely` use it by default (`engine="auto"`); the
`Walkthrough` and the lessons stay in Python.

```python
import square_rulerpack as sq
sq.HAVE_NATIVE                                                   # True
sq.factor_completely(n, workers=8)                               # Rust, 8 rayon threads, GIL released
sq.factor_completely(n, workers=8, engine="python")              # this port, for comparison
sq.factor_completely(n, workers=8, verbose=True)                 # the engine's own progress on stderr
sq.native.factor_completely(n, threads=8, progress=lambda stage, full, need, partials, polys: ...)
sq.native.is_prime(2**89 - 1)
```

`Settings` still applies: its fields are exported as the engine's `SIQS_*` variables for the
duration of the call. Python ints cross as big integers; there is no string round trip. Ctrl-C
does not stop a running sieve (the engine has no cancel hook).

Measured on an i9-11900H (8 cores), 8 threads, quiet machine:


sample 80 digit: 50726044447758565750355444206008211963939957754977146054012847169421070122136537
sample 60 digit: 100433627766186892221372630609062766858404681029709092356097

| input | Rust extension | this Python port | `sympy.factorint` |
|---|---|---|---|
| 80-digit semiprime 5072…6537 (from an IPython session) | **15.5 s** | an hour-ish | not attempted |
| 40-digit semiprime, two random 20-digit primes | 0.01 s | 0.27 s | 7.8 s |
| 50-digit semiprime, two random 25-digit primes | 0.06 s | 0.33 s | 8.8 s |
| 60-digit semiprime, two random 30-digit primes | 0.80 s | 14.8 s | 116.6 s |

using the workers=1 acheive theabove results, the 80 digit cna be down in 19  minutes with more workers.

`sympy.factorint` has Pollard rho, p−1 and ECM but no quadratic sieve, so its time grows with the
size of the *smaller factor*; the sieve's grows with the size of *n* and does not care how the
factors are split. For balanced semiprimes above 40 digits the difference is minutes versus
fractions of a second.

Building from source instead of the wheel (`pip install square-rulerpack-native --no-binary :all:`
or `maturin develop --release` in `native/`) needs a Rust toolchain; `.github/workflows/native-wheels.yml`
builds the Linux, macOS and Windows wheels.

## Speed

Measured on an i9-11900H (8 cores), numba kernel, random semiprimes:

| digits | workers | time | note |
|---|---|---|---|
| 16 | 1 | 0.3 s | the Rust crate's README example |
| 40 | 1 | 0.2 s | |
| 50 | 1 | 1.4 s | |
| 60 (197 bits) | 2 | 3.7 s | hand table row (base 3000, M = 2^16); Rust 0.26 s at 8 threads |
| 60 (199 bits) | 1 | 66 s | YAFU table row (base 5744, M = 2^15, slack 1.3) |
| 60 (199 bits) | 4 | 23 s | same row |
| 66 | 4 | 54 s | exact threshold |
| 70 | 8 | 62 s | measured on a busy machine, see `TIMINGS.md` |

Python is roughly 30–60× slower than the Rust engine per polynomial. The point of this port
is the walkthrough; for the 90- and 100-digit runs use the Rust crate (`cargo install
square-rulerpack`) or the browser page. Above about 75 digits the dense GF(2) elimination also
becomes the long pole (it is cubic; the Rust engine switches to block Lanczos at 6,000 rows).

## Files

- `square_rulerpack/arith.py` — primes, Legendre, Tonelli–Shanks, Miller–Rabin, Pollard–Brent, xorshift64.
- `square_rulerpack/params.py` — the YAFU direct-sieve table and the hand rows, `Params.for_n`.
- `square_rulerpack/engine.py` — multiplier, factor base, families, `FamilySieve` (sieve, scan,
  candidates), `Cycles`, `relations`, `square_root`, `factor_completely`.
- `square_rulerpack/linalg.py` — singleton pruning and bit-packed Gaussian elimination.
- `square_rulerpack/workers.py` — the process pool; families merged in index order, so the
  relation set is the same at any worker count.
- `square_rulerpack/walkthrough.py` — `Walkthrough` and `SieveResult`.
- `square_rulerpack/lessons.py` — `construct_polynomial`, `LessonWindow`, `factor_walkthrough`,
  `relation_study`, `engine_polynomials`.
- `square_rulerpack/native.py` — loader for the Rust extension, `Settings` → `SIQS_*` bridge, wrappers.
- `native/` — the PyO3 crate (`square-rulerpack-native`): `Cargo.toml`, `pyproject.toml` (maturin), `src/lib.rs`.
- `tests/test_factor.py`, `tests/test_lessons.py`, `tests/test_native.py` — `pytest tests` (25 tests, a few seconds; the native tests skip when the extension is absent).

## License

MIT, as the Rust crate.
