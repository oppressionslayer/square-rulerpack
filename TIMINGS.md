# Timings of the Python port

Machine: i9-11900H (8 cores / 16 threads), Python 3.14, numpy 2.5, numba 0.68. Random semiprimes
(two primes of equal size). "quiet" = nothing else running; "busy" = measured while the browser
lab was factoring a 100-digit number on 8 workers (Chrome at 830 % CPU), so those rows are roughly
half of what a quiet machine gives.

| digits | bits | table row | threshold | workers | polynomials | polys/s | sieving | total | machine |
|---|---|---|---|---|---|---|---|---|---|
| 16 | 50 | hand (base 150, M 2^13) | plain | 1 | 2 | — | 0.3 s | 0.34 s | quiet |
| 30 | 99 | hand (base 150) | plain | 1 | 60 | — | 0.0 s | 0.05 s | quiet |
| 40 | 133 | hand (base 400, M 2^14) | plain | 1 | 608 | 3,000 | 0.2 s | 0.22 s | quiet |
| 50 | 165 | hand (base 1100, M 2^15) | plain | 1 | 2,368 | 1,700 | 1.4 s | 1.44 s | quiet |
| 60 | 197 | hand (base 3000, M 2^16) | plain | 2 | 9,600 | 3,100 | 3.1 s | 3.65 s | quiet |
| 60 | 199 | YAFU (base 5744, M 2^15, slack 1.3) | plain | 1 | 7,872 | 126 | 62.6 s | 65.6 s | quiet |
| 60 | 199 | YAFU (base 5744) | plain | 4 | 9,920 | 516 | 19.2 s | 22.8 s | quiet |
| 66 | 219 | YAFU (base 9664) | exact | 4 | 38,272 | 1,100 | 34.8 s | 54.4 s | quiet |
| 60 | 198 | YAFU (base 5472) | plain | 1 | 8,256 | 104 | 79.2 s | 83.9 s | busy |
| 66 | 217 | YAFU (base 9888) | exact | 4 | 31,616 | 469 | 67.4 s | 99.4 s | busy |
| 70 | 230 | YAFU (base 8432) | exact | 8 | 87,168 | 2,200 | 39.6 s | 61.6 s | busy |

Single-family rate on one 66-digit number, 128 polynomials, one process: 267 polys/s quiet,
140 polys/s busy. The Rust engine does about 8,000 polys/s per thread at this size, so the
Python port is 30–60× slower per polynomial.

Where the Python time goes (cProfile, 66 digits, exact threshold, one family): about 45 % in
the candidate stage (big-integer division of each survivor, Pollard–Brent for two-large-prime
cofactors), 25 % in the exact report test of the survivors (small-prime valuations), the sieve
kernel itself under 5 %. The rest of a run is the matrix: dense GF(2) elimination is cubic, and
took 20 s at 8,500 rows and 32 s at 10,000 rows (busy machine). The Rust engine switches to
block Lanczos above 6,000 rows; the Python port does not have Lanczos, so above about 75 digits
the matrix is the long pole.

The YAFU table rows (slack 1.3, which lets far more survivors through to exact division) are
tuned for the Rust engine's cheap candidate stage. In Python the 60-digit number that happened to
fall under the hand row (197 bits, base 3000, slack 0.7) was 20× faster than the one under the
YAFU row (199 bits). A Python-tuned table (smaller slack, larger M) would help; it was not done
because the point of this port is the walkthrough, not speed.

## The Rust engine as an extension (`native/`, package `square-rulerpack-native`)

Same machine, quiet, 8 threads, measured through `python -m square_rulerpack … --workers 8` with
the extension installed (`engine=auto` picks it up):

| input | Rust extension | notes |
|---|---|---|
| 80-digit semiprime 5072…6537 (265 bits, kN 269 bits, multiplier 13, base 11,856, Q1) | 15.5 s wall, 266 MB peak | factors 6245…1527 × 8122…9631, both probable prime, product checks |
| the same number through the Rust CLI (`target/release/examples/factor`, `SIQS_Q2=1 SIQS_SMALL=512 RAYON_NUM_THREADS=8`) | 15.6 s wall | same factors; the extension adds no measurable overhead |

Random balanced semiprimes, Rust extension vs this Python port vs `sympy.factorint` (sympy in a
subprocess, stopped at 300 s):

| digits | n | Rust extension, 8 threads | Python port, 8 workers | `sympy.factorint` |
|---|---|---|---|---|
| 40 | 6213204885379705444888466767095804123277 = 75583987530970793191 × 82202660753163146347 | 0.01 s | 0.27 s | 7.8 s |
| 50 | 76026501459413395163330630762096504839507238065451 = 8410324106821786346161301 × 9039663691170562716009151 | 0.06 s | 0.33 s | 8.8 s |
| 60 | 230183077131629724616750314912636607708734569350788817355159 = 324178577113326424999993911721 × 710050242003382811814015785279 | 0.80 s | 14.8 s | 116.6 s |

(The 40- and 50-digit numbers fall under the hand-tuned parameter rows, where the Python port with
8 workers is within a factor of 30 of Rust; the 60-digit one lands on the YAFU row (slack 1.3) and
is 18× slower. sympy's time at 60 digits is its ECM finding a 30-digit factor; at 70 digits it
would be hours. An earlier draft of this table had the Python column wrong because
`factor_completely(engine="python")` re-entered the default engine on the recursive split; fixed
and re-measured.)

One thing the native tests caught: `rulerpack 0.3.1`'s `big::is_prime_big` returns false for the
odd primes 5 through 97, because each of those is one of its own Miller–Rabin witnesses and a
witness ≡ 0 (mod n) fails the strong test. Composites are never called prime and primes ≥ 101
are fine. The extension routes n < 2^64 through the engine's deterministic `is_prime_u64`
instead; the crate itself should skip witnesses that divide n (or special-case n < 100).
