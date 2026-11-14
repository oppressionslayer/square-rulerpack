//! square-rulerpack for Python. Python ints cross as `BigUint` (pyo3's `num-bigint` feature), the
//! sieve runs with the GIL released on a rayon pool of the requested size, and the engine's
//! progress hook can call back into Python.
//!
//! The engine reads its switches (`SIQS_Q2`, `SIQS_SMALL`, `SIQS_EXACT`, …) from the process
//! environment at call time; the pure-Python package sets them from its `Settings` before calling.
use std::sync::Mutex;

use num_bigint::BigUint;
use num_traits::{One, Zero};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

/// The Python callable the progress hook forwards to (set for the duration of one call).
static PROGRESS: Mutex<Option<Py<PyAny>>> = Mutex::new(None);

/// The engine's hook is a plain `fn` pointer kept in a thread-local, so it fires on the thread
/// that called `factor`; that thread has released the GIL, so we re-attach to call Python.
fn py_hook(stage: u32, full: usize, need: usize, partials: usize, polys: usize) {
    let cb = PROGRESS.lock().ok().and_then(|g| g.as_ref().map(|c| c.clone_ref_detached()));
    if let Some(cb) = cb {
        Python::attach(|py| {
            if let Err(e) = cb.call1(py, (stage, full, need, partials, polys)) {
                e.print(py);
            }
        });
    }
}

trait CloneDetached {
    fn clone_ref_detached(&self) -> Py<PyAny>;
}
impl CloneDetached for Py<PyAny> {
    fn clone_ref_detached(&self) -> Py<PyAny> {
        Python::attach(|py| self.clone_ref(py))
    }
}

fn pool(threads: Option<usize>) -> PyResult<rayon::ThreadPool> {
    let t = match threads {
        Some(0) | None => std::thread::available_parallelism().map(|n| (n.get() / 2).max(1)).unwrap_or(1),
        Some(t) => t,
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(t)
        .build()
        .map_err(|e| PyRuntimeError::new_err(format!("rayon pool: {e}")))
}

fn with_progress<T>(py: Python<'_>, progress: Option<Py<PyAny>>, threads: Option<usize>,
                    f: impl FnOnce() -> T + Send) -> PyResult<T>
where
    T: Send,
{
    let pool = pool(threads)?;
    let hooked = progress.is_some();
    if let Some(cb) = progress {
        *PROGRESS.lock().map_err(|_| PyRuntimeError::new_err("progress lock"))? = Some(cb);
    }
    let out = py.detach(|| {
        pool.install(|| {
            if hooked {
                square_rulerpack::siqs::set_progress_hook(Some(py_hook));
            }
            let r = f();
            square_rulerpack::siqs::set_progress_hook(None);
            r
        })
    });
    if hooked {
        if let Ok(mut g) = PROGRESS.lock() {
            *g = None;
        }
    }
    Ok(out)
}

fn check(n: &BigUint) -> PyResult<()> {
    if n.is_zero() {
        return Err(PyValueError::new_err("n must be a positive integer"));
    }
    Ok(())
}

/// The complete factorisation of n: sorted prime factors with multiplicity (small primes by
/// trial division, perfect powers, then the self-initialising quadratic sieve on every composite
/// piece). `threads` is the rayon pool size (default half the hardware threads). `progress`, if
/// given, is called as progress(stage, full, need, partials, polys) with stage 0 preparing,
/// 1 collecting, 2 matrix, 3 square root. The GIL is released while the sieve runs; Ctrl-C does
/// not interrupt it.
#[pyfunction]
#[pyo3(signature = (n, threads=None, progress=None))]
fn factor_completely(py: Python<'_>, n: BigUint, threads: Option<usize>, progress: Option<Py<PyAny>>) -> PyResult<Vec<BigUint>> {
    check(&n)?;
    if n.is_one() {
        return Ok(vec![]);
    }
    with_progress(py, progress, threads, move || square_rulerpack::siqs::factor_completely(&n))
}

/// One nontrivial factor of n from the sieve alone, or None (prime, a prime power the sieve
/// does not handle, too small, or unlucky). Use `factor_completely` for the whole factorisation.
#[pyfunction]
#[pyo3(signature = (n, threads=None, progress=None))]
fn factor(py: Python<'_>, n: BigUint, threads: Option<usize>, progress: Option<Py<PyAny>>) -> PyResult<Option<BigUint>> {
    check(&n)?;
    with_progress(py, progress, threads, move || square_rulerpack::siqs::factor(&n))
}

/// Primality: deterministic Miller–Rabin below 2^64 (the engine's `is_prime_u64`), rulerpack's
/// strong probable-prime test above. (rulerpack 0.3.1's `is_prime_big` returns false for the
/// primes it uses as witnesses, 2..37, so small n must not go through it.)
#[pyfunction]
fn is_prime(n: BigUint) -> bool {
    match u64::try_from(&n) {
        Ok(v) => square_rulerpack::siqs::is_prime_u64(v),
        Err(_) => rulerpack::big::is_prime_big(&n),
    }
}

/// The version of this binding.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// The version of the square-rulerpack engine crate linked in.
#[pyfunction]
fn engine_version() -> &'static str {
    "0.2.0"
}

/// Whether the engine will use its AVX2 kernels on this machine.
#[pyfunction]
fn have_avx2() -> bool {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        std::env::var("SIQS_NO_AVX2").is_err() && std::arch::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        false
    }
}

#[pymodule]
fn square_rulerpack_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(factor_completely, m)?)?;
    m.add_function(wrap_pyfunction!(factor, m)?)?;
    m.add_function(wrap_pyfunction!(is_prime, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(engine_version, m)?)?;
    m.add_function(wrap_pyfunction!(have_avx2, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
