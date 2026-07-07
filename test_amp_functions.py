#!/usr/bin/env python3
"""
test_amp_functions.py
Per-function comparison of Fortran and C amplitude evaluation libraries.

For each generated (igroup, iint) pair, calls both the Fortran and C
dispatcher with identical momenta and compares the output amplitudes.

By default, uses the single reference momentum stored in
Library/amp{G}_{I}_lib.data.  Use --n-random to additionally test with
random phase-space points generated from a simple isotropic prescription
(uniform on the mass shell with E_cm=13000 GeV), which catches differences
that only appear at specific kinematic configurations.

Usage:
  python3 test_amp_functions.py                  # reference point only
  python3 test_amp_functions.py --verbose        # show all results
  python3 test_amp_functions.py --n-random 100   # also test 100 random points
  python3 test_amp_functions.py --tol 1e-12
  python3 test_amp_functions.py --data-only      # compare .data files, no library load
"""

import ctypes
import numpy as np
import re
import sys
from pathlib import Path


# Number of (pair, momentum) jobs sent to the GPU per CUDA call. Lower this to
# get more frequent partial-result output at the cost of extra host<->device
# round-trips; raise it to minimise overhead at the cost of only seeing
# results once the whole batch finishes.
CU_BATCH_SIZE = 1


# ─── helpers ────────────────────────────────────────────────────────────────

def parse_dispatcher_dims(lib_dir: Path):
    """Extract max_next and max_amps from the generated amplibc.c dispatcher."""
    text = (lib_dir / 'amplibc.c').read_text()
    m_next = re.search(r'AC_D_FP p_arr\[(\d+)\]\[4\]', text)
    m_amps = re.search(r'AC_D_CX amps_arr\[(\d+)\]', text)
    if not m_next or not m_amps:
        raise RuntimeError("Cannot parse dispatcher dimensions from Library/amplibc.c")
    return int(m_next.group(1)), int(m_amps.group(1))


def parse_cuda_dims(lib_dir: Path):
    """Extract AC_NEXT and AC_NAMP from the generated amplib.cuh CUDA header."""
    text = (lib_dir / 'amplib.cuh').read_text()
    m_next = re.search(r'#define\s+AC_NEXT\s+(\d+)', text)
    m_amps = re.search(r'#define\s+AC_NAMP\s+(\d+)', text)
    if not m_next or not m_amps:
        raise RuntimeError("Cannot parse CUDA dispatcher dimensions from Library/amplib.cuh")
    return int(m_next.group(1)), int(m_amps.group(1))


def parse_pair_dims(hpp_path: Path, G: int, I: int):
    """Extract (next, n_amps) for pair (G,I) from its generated .h header."""
    text = hpp_path.read_text()
    pat = (rf'evaluate_amp{G}_{I}\s*\('
           rf'const AC_D_FP p\[(\d+)\]\[4\],\s*'
           rf'AC_D_CX amps\[(\d+)\]\)')
    m = re.search(pat, text)
    if not m:
        raise RuntimeError(f"Cannot find signature for amp{G}_{I} in {hpp_path.name}")
    return int(m.group(1)), int(m.group(2))


def read_fort_data(path: Path, next_val: int):
    """
    Read Fortran stream-binary data file.
    Layout: real(kind=8) p(0:3, next)  →  4*next float64 values
            complex(kind=8) amps(n)    →  n complex128 values
    Fortran column-major p(comp, part) has the same flat layout as C
    row-major double[part][comp], so no reordering is needed.
    """
    raw = np.fromfile(path, dtype=np.float64)
    n_p = 4 * next_val
    p_flat = raw[:n_p].copy()
    amps = raw[n_p:].view(np.complex128).copy()
    return p_flat, amps


def find_pairs(lib_dir: Path):
    pairs = []
    for f in sorted(lib_dir.glob('amp*_lib.data')):
        m = re.match(r'amp(\d+)_(\d+)_lib$', f.stem)
        if m:
            pairs.append((int(m.group(1)), int(m.group(2)), f))
    return pairs


def format_complex(z, spec: str = '') -> str:
    """Format a complex number as '(real, imag j)', with an optional
    per-component format spec (e.g. '.4e')."""
    return f"({format(z.real, spec)}, {format(z.imag, spec)} j)"


def rel_err(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    nan_mask = ~np.isfinite(a) | ~np.isfinite(b)
    denom = np.abs(a) + np.abs(b)
    result = np.where(denom > 0, np.abs(a - b) / denom, 0.0)
    result[nan_mask] = np.inf
    return result


# ─── data-file-only mode ────────────────────────────────────────────────────

def run_data_comparison(lib_dir: Path, pairs, tol: float, verbose: bool):
    """Compare Library/amp{G}_{I}_lib.data vs amp{G}_{I}_libc.data directly."""
    print("Mode: direct data-file comparison (no library loading)\n")
    n_pass = n_fail = n_skip = 0
    max_err_global = 0.0
    failures = []

    for G, I, fort_path in pairs:
        c_path = lib_dir / f'amp{G}_{I}_libc.data'
        label = f'amp{G}_{I}'

        hpp_path = lib_dir / f'amp{G}_{I}_lib.h'
        if not hpp_path.exists():
            n_skip += 1
            continue
        try:
            next_val, n_amps = parse_pair_dims(hpp_path, G, I)
        except RuntimeError:
            n_skip += 1
            continue

        if not c_path.exists():
            if verbose:
                print(f"  SKIP  {label:<12}  (no _libc.data file)")
            n_skip += 1
            continue

        p_fort, amps_fort = read_fort_data(fort_path, next_val)
        p_c,  amps_c  = read_fort_data(c_path,  next_val)

        # Check that the reference momenta match (both files should use the same p)
        if not np.allclose(p_fort, p_c, rtol=1e-15, atol=0):
            print(f"  WARN  {label:<12}  momenta differ between _lib.data and _libc.data")

        err = rel_err(amps_fort[:n_amps], amps_c[:n_amps])
        max_err = float(np.max(err))
        max_err_global = max(max_err_global, max_err)
        passed = max_err < tol

        if passed:
            n_pass += 1
            if verbose:
                print(f"  PASS  {label:<12}  max_rel_err={max_err:.2e}")
        else:
            n_fail += 1
            worst = int(np.argmax(err))
            failures.append((G, I, max_err, worst, amps_fort[worst], amps_c[worst]))
            print(f"  FAIL  {label:<12}  max_rel_err={max_err:.2e}  worst_amp={worst+1}")
            if verbose:
                print(f"         Fortran: {format_complex(amps_fort[worst])}")
                print(f"         C:       {format_complex(amps_c[worst])}")

    return n_pass, n_fail, n_skip, max_err_global, failures


# ─── live-dispatch mode ──────────────────────────────────────────────────────

def run_dispatch_comparison(lib_dir: Path, pairs, tol: float, verbose: bool,
                             random_momenta=None):
    """Load both .so files and call evaluate_amp for each pair with identical momenta."""
    max_next, max_amps = parse_dispatcher_dims(lib_dir)
    print(f"Dispatcher: max_next={max_next}, max_amps={max_amps}")

    try:
        fort_so = ctypes.CDLL('./libamplicolmadspace_f.so')
        c_so  = ctypes.CDLL('./libamplicolmadspace_c.so')
    except OSError as e:
        sys.exit(f"ERROR: cannot load libraries: {e}\n"
                 f"Try --data-only to compare saved data files instead.")

    # Fortran: __amp_lib_MOD_evaluate_amp(ichan*, iint*, p(0:3,max_next)*, amps(*)*)
    fort_eval = getattr(fort_so, '__amp_lib_MOD_evaluate_amp')
    fort_eval.restype = None

    c_eval = c_so.evaluate_amp
    c_eval.restype = None

    print(f"Loaded Fortran and C dispatcher libraries\n")

    if random_momenta is None:
        random_momenta = []

    n_pass = n_fail = n_skip = 0
    max_err_global = 0.0
    failures = []

    for G, I, data_path in pairs:
        hpp_path = lib_dir / f'amp{G}_{I}_lib.h'
        label = f'amp{G}_{I}'

        if not hpp_path.exists():
            n_skip += 1
            continue
        try:
            next_val, n_amps = parse_pair_dims(hpp_path, G, I)
        except RuntimeError as e:
            print(f"  SKIP  {label}: {e}")
            n_skip += 1
            continue

        # Collect all momenta to test: the stored reference + random points
        p_ref, _ = read_fort_data(data_path, next_val)
        momenta_list = [p_ref]
        for rng_p in random_momenta:
            momenta_list.append(rng_p[:4 * next_val])

        pair_max_err = 0.0
        pair_worst   = (0, None, None)
        ichan = ctypes.c_int(G)
        iint  = ctypes.c_int(I)

        for p_flat in momenta_list:
            p_padded = np.zeros(4 * max_next, dtype=np.float64)
            p_padded[:4 * next_val] = p_flat

            amps_fort = np.zeros(max_amps, dtype=np.complex128)
            amps_c  = np.zeros(max_amps, dtype=np.complex128)
            p_ptr = p_padded.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

            fort_eval(
                ctypes.byref(ichan), ctypes.byref(iint),
                p_ptr,
                amps_fort.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )
            c_eval(
                ctypes.byref(ichan), ctypes.byref(iint),
                p_ptr,
                amps_c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )

            af = amps_fort[:n_amps]
            ac = amps_c[:n_amps]
            err = rel_err(af, ac)
            max_err = float(np.max(err))
            if max_err >= pair_max_err:
                pair_max_err = max_err
                worst_idx = int(np.argmax(err))
                pair_worst = (worst_idx, af[worst_idx], ac[worst_idx])

        max_err_global = max(max_err_global, pair_max_err)
        passed = pair_max_err <= tol

        if passed:
            n_pass += 1
            if verbose:
                print(f"  PASS  {label:<12}  max_rel_err={pair_max_err:.2e}")
        else:
            n_fail += 1
            failures.append((G, I, pair_max_err, *pair_worst))
            idx, af_w, ac_w = pair_worst
            vals_str = (f"  fort={format_complex(af_w, '.4e')}  c={format_complex(ac_w, '.4e')}"
                        if af_w is not None else "")
            print(f"  FAIL  {label:<12}  max_rel_err={pair_max_err:.2e}  worst_amp={idx+1}"
                  + vals_str)

    return n_pass, n_fail, n_skip, max_err_global, failures


# ─── CUDA live-dispatch mode ─────────────────────────────────────────────────

def run_cu_dispatch_comparison(lib_dir: Path, pairs, tol: float, verbose: bool,
                                random_momenta=None):
    """
    Load libamplicolmadspace_cu.so (if it exists) and compare against the
    Fortran dispatcher for each pair.  Returns None if the CUDA library is
    unavailable or fails to load.

    All (pair, momentum) combinations are split into sequential batches of
    CU_BATCH_SIZE jobs, each dispatched as one CUDA call, trading runtime
    for more frequent partial-result output.  Fortran evaluations still run
    sequentially since the Fortran dispatcher has no batch interface.
    """
    cu_so_path = Path('./libamplicolmadspace_cu.so')
    if not cu_so_path.exists():
        return None

    try:
        cu_so = ctypes.CDLL(str(cu_so_path))
    except OSError as e:
        print(f"Note: CUDA library found but failed to load: {e}")
        return None

    try:
        max_next_cu, max_amps_cu = parse_cuda_dims(lib_dir)
    except RuntimeError as e:
        print(f"Note: {e}")
        return None

    try:
        fort_so = ctypes.CDLL('./libamplicolmadspace_f.so')
    except OSError as e:
        sys.exit(f"ERROR: cannot load Fortran library: {e}")

    max_next_f, max_amps_f = parse_dispatcher_dims(lib_dir)

    fort_eval = getattr(fort_so, '__amp_lib_MOD_evaluate_amp')
    fort_eval.restype = None

    cu_eval = cu_so.evaluate_amp
    cu_eval.restype = None

    print(f"CUDA dispatcher: AC_NEXT={max_next_cu}, AC_NAMP={max_amps_cu}")
    print(f"Loaded CUDA library — comparing against Fortran\n")

    if random_momenta is None:
        random_momenta = []

    # ── Phase 1: enumerate all valid (pair, momentum) jobs ───────────────────
    # pair_meta[i] is None for skipped pairs, else (G, I, n_amps, job_start, n_jobs)
    pair_meta = []
    # all_jobs[k] = (G, I, next_val, n_amps, p_flat)
    all_jobs = []

    for G, I, data_path in pairs:
        hpp_path = lib_dir / f'amp{G}_{I}_lib.h'
        if not hpp_path.exists():
            pair_meta.append(None)
            continue
        try:
            next_val, n_amps = parse_pair_dims(hpp_path, G, I)
        except RuntimeError:
            pair_meta.append(None)
            continue

        p_ref, _ = read_fort_data(data_path, next_val)
        momenta_list = [p_ref] + [rng_p[:4 * next_val] for rng_p in random_momenta]

        job_start = len(all_jobs)
        for p_flat in momenta_list:
            all_jobs.append((G, I, next_val, n_amps, p_flat))
        pair_meta.append((G, I, n_amps, job_start, len(momenta_list)))

    total = len(all_jobs)

    # ── Phase 2: Fortran — sequential (no batch interface) ───────────────────
    fort_results = []   # fort_results[k] = complex128 array of length n_amps
    ichan_c = ctypes.c_int(0)
    iint_c  = ctypes.c_int(0)
    for G, I, next_val, n_amps, p_flat in all_jobs:
        p_f = np.zeros(4 * max_next_f, dtype=np.float64)
        p_f[:4 * next_val] = p_flat
        amps_f = np.zeros(max_amps_f, dtype=np.complex128)
        ichan_c.value = G
        iint_c.value  = I
        fort_eval(
            ctypes.byref(ichan_c), ctypes.byref(iint_c),
            p_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            amps_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
        fort_results.append(amps_f[:n_amps].copy())

    # ── Phase 3: CUDA — sequential batches of CU_BATCH_SIZE jobs ─────────────
    amps_batch = np.zeros(total * max_amps_cu, dtype=np.complex128)
    n_batches = (total + CU_BATCH_SIZE - 1) // CU_BATCH_SIZE

    for b in range(n_batches):
        start = b * CU_BATCH_SIZE
        end = min(start + CU_BATCH_SIZE, total)
        batch_jobs = all_jobs[start:end]
        n_batch = end - start

        ichans  = np.array([j[0] for j in batch_jobs], dtype=np.int32)
        iints   = np.array([j[1] for j in batch_jobs], dtype=np.int32)
        p_batch = np.zeros(n_batch * 4 * max_next_cu, dtype=np.float64)
        for k, (_, _, next_val, _, p_flat) in enumerate(batch_jobs):
            p_batch[k * 4 * max_next_cu : k * 4 * max_next_cu + 4 * next_val] = p_flat
        amps_out = np.zeros(n_batch * max_amps_cu, dtype=np.complex128)

        cu_eval(
            ichans.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            iints.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            p_batch.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            amps_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int(n_batch),
        )

        amps_batch[start * max_amps_cu : end * max_amps_cu] = amps_out

        if verbose:
            print(f"  CUDA batch {b + 1}/{n_batches}: jobs {start}-{end - 1} "
                  f"({n_batch} events)")

    # ── Phase 4: compare per pair ─────────────────────────────────────────────
    n_pass = n_fail = n_skip = 0
    max_err_global = 0.0
    failures = []

    for pair_idx, (G, I, *_) in enumerate(pairs):
        label = f'amp{G}_{I}'
        meta = pair_meta[pair_idx]
        if meta is None:
            n_skip += 1
            continue
        _, _, n_amps, job_start, n_jobs = meta

        pair_max_err = 0.0
        pair_worst   = (0, None, None)

        for k in range(n_jobs):
            job_idx = job_start + k
            af  = fort_results[job_idx]
            acu = amps_batch[job_idx * max_amps_cu : job_idx * max_amps_cu + n_amps]
            err = rel_err(af, acu)
            max_err = float(np.max(err))
            if max_err >= pair_max_err:
                pair_max_err = max_err
                worst_idx = int(np.argmax(err))
                pair_worst = (worst_idx, af[worst_idx], acu[worst_idx])

        max_err_global = max(max_err_global, pair_max_err)
        passed = pair_max_err <= tol

        if passed:
            n_pass += 1
            if verbose:
                print(f"  PASS  {label:<12}  max_rel_err={pair_max_err:.2e}")
        else:
            n_fail += 1
            failures.append((G, I, pair_max_err, *pair_worst))
            idx, af_w, acu_w = pair_worst
            vals_str = (f"  fort={format_complex(af_w, '.4e')}  cu={format_complex(acu_w, '.4e')}"
                        if af_w is not None else "")
            print(f"  FAIL  {label:<12}  max_rel_err={pair_max_err:.2e}  worst_amp={idx+1}"
                  + vals_str)

    return n_pass, n_fail, n_skip, max_err_global, failures

def make_random_momenta(n_points: int, max_next: int, ecm: float = 13000.0, seed: int = 42):
    """
    Generate n_points flat momentum arrays of shape (4*max_next,) for massless
    particles with overall energy-momentum conservation.  Uses a simple
    back-to-back + isotropic decay prescription; the exact kinematics don't
    matter as long as they are physical.

    Returns list of 1-D float64 arrays, each of length 4*max_next.
    p_flat layout: p(comp=0..3, part=1..next) column-major = particle-major.
    """
    rng = np.random.default_rng(seed)
    result = []
    for _ in range(n_points):
        n = max_next
        # Random massless momenta summing to (ecm, 0, 0, 0)
        # Use RAMBO-like flat generation: exponential energies, isotropic angles
        u = rng.random(n)
        cos_t = 2.0 * rng.random(n) - 1.0
        phi   = 2.0 * np.pi * rng.random(n)
        sin_t = np.sqrt(1.0 - cos_t**2)
        q = np.column_stack([
            -np.log(u * rng.random(n)),   # energies (unnormalised)
            sin_t * np.cos(phi),
            sin_t * np.sin(phi),
            cos_t,
        ])
        q[:, 1:] *= q[:, 0:1]            # 3-momenta = E * direction

        # Boost to sum to (ecm/2, 0, 0, 0) for the first two (incoming),
        # and distribute remainder among outgoing.  For simplicity, rescale
        # all energies so that sum(E) = ecm with incoming along z.
        scale = ecm / q[:, 0].sum()
        q *= scale

        # Pack as Fortran p(0:3, next) column-major
        p_flat = q.T.ravel().astype(np.float64)   # shape (4*next,)
        result.append(p_flat)
    return result


# ─── main ────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--verbose', '-v', action='store_true')
    ap.add_argument('--tol', type=float, default=0.0,
                    help='Relative error tolerance (default 0.0 = exact match)')
    ap.add_argument('--data-only', action='store_true',
                    help='Compare saved .data files; do not load shared libraries')
    ap.add_argument('--n-random', type=int, default=0,
                    help='Number of additional random phase-space points to test (default 0)')
    ap.add_argument('--seed', type=int, default=42,
                    help='RNG seed for random momenta (default 42)')
    ap.add_argument('--libdir', default='Library')
    args = ap.parse_args()

    lib_dir = Path(args.libdir)
    pairs = find_pairs(lib_dir)
    print(f"Found {len(pairs)} (igroup, iint) pairs in {lib_dir}/\n")

    n_fail_total = 0

    if args.data_only:
        n_pass, n_fail, n_skip, max_err, failures = run_data_comparison(
            lib_dir, pairs, args.tol, args.verbose)
        n_fail_total = n_fail
        print(f"\n{'='*70}")
        print(f"Results: {n_pass} passed, {n_fail} failed, {n_skip} skipped  "
              f"(tol={args.tol:.0e})")
        print(f"Global max relative error: {max_err:.2e}")
        if failures:
            print(f"\nFailing pairs (sorted by error magnitude):")
            for G, I, err, idx, af_w, ac_w in sorted(failures, key=lambda x: -x[2]):
                print(f"  amp{G}_{I:<4}  err={err:.2e}  amp_index={idx+1}")
                print(f"             Fortran: {format_complex(af_w)}")
                print(f"             C:      {format_complex(ac_w)}")
    else:
        max_next, _ = parse_dispatcher_dims(lib_dir)
        random_momenta = (make_random_momenta(args.n_random, max_next, seed=args.seed)
                          if args.n_random > 0 else [])
        if random_momenta:
            print(f"Testing {1 + len(random_momenta)} momenta per pair "
                  f"(1 reference + {len(random_momenta)} random)\n")

        # ─── Fortran vs C ──────────────────────────────────────────────────
        n_pass, n_fail, n_skip, max_err, failures = run_dispatch_comparison(
            lib_dir, pairs, args.tol, args.verbose, random_momenta)
        n_fail_total += n_fail
        print(f"\n{'='*70}")
        print(f"Fortran vs C   {n_pass} passed, {n_fail} failed, {n_skip} skipped  "
              f"(tol={args.tol:.0e})   max_rel_err={max_err:.2e}")
        if failures:
            print(f"\nFailing pairs — Fortran vs C (sorted by error magnitude):")
            for G, I, err, idx, af_w, ac_w in sorted(failures, key=lambda x: -x[2]):
                print(f"  amp{G}_{I:<4}  err={err:.2e}  amp_index={idx+1}")
                print(f"             Fortran: {format_complex(af_w)}")
                print(f"             C:       {format_complex(ac_w)}")

        # ─── Fortran vs CUDA (skipped if library not present) ────────────────
        cu_result = run_cu_dispatch_comparison(
            lib_dir, pairs, args.tol, args.verbose, random_momenta)
        if cu_result is not None:
            cu_pass, cu_fail, cu_skip, cu_max_err, cu_failures = cu_result
            n_fail_total += cu_fail
            print(f"\n{'='*70}")
            print(f"Fortran vs CUDA  {cu_pass} passed, {cu_fail} failed, {cu_skip} skipped  "
                  f"(tol={args.tol:.0e})   max_rel_err={cu_max_err:.2e}")
            if cu_failures:
                print(f"\nFailing pairs — Fortran vs CUDA (sorted by error magnitude):")
                for G, I, err, idx, af_w, acu_w in sorted(cu_failures, key=lambda x: -x[2]):
                    print(f"  amp{G}_{I:<4}  err={err:.2e}  amp_index={idx+1}")
                    print(f"             Fortran: {format_complex(af_w)}")
                    print(f"             CUDA:    {format_complex(acu_w)}")

    return 1 if n_fail_total > 0 else 0


if __name__ == '__main__':
    sys.exit(main())
