# square-rulerpack-native

Lars Rocha's [square-rulerpack](https://crates.io/crates/square-rulerpack) factoring engine
(self-initialising quadratic sieve with double large primes, AVX2 kernels, rayon threads, in Rust)
compiled as a Python extension module with PyO3. It is the fast path of the pure-Python
[square-rulerpack](https://pypi.org/project/square-rulerpack/) package, which does the teaching;
this module does the factoring.

```sh
pip install square-rulerpack[native]      # or: pip install square-rulerpack-native
```

```python
import square_rulerpack_native as srn
srn.factor_completely(50726044447758565750355444206008211963939957754977146054012847169421070122136537, threads=8)
srn.factor(n, threads=8)                  # one nontrivial factor from the sieve, or None
srn.is_prime(2**61 - 1)                   # strong probable-prime test
srn.factor_completely(n, threads=8, progress=lambda stage, full, need, partials, polys: print(stage, full, need))
```

Through the pure package it is the default engine once installed:

```python
import square_rulerpack as sq
sq.factor_completely(n, workers=8)        # engine="auto": Rust when installed, Python otherwise
```

Python ints cross as big integers (no string round trip). The GIL is released while the sieve
runs; Ctrl-C does not interrupt it. The engine's switches (`SIQS_Q2`, `SIQS_SMALL`, `SIQS_EXACT`,
`SIQS_DLP`, …) are read from the environment at call time, exactly as the Rust CLI reads them;
the pure package sets them from its `Settings`. Wheels are abi3 (one wheel per platform for
every Python ≥ 3.10); building from the sdist needs a Rust toolchain.

MIT, as the crate.
