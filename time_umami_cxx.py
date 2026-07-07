import madspace as ms
import numpy as np
import time

np.random.seed(12345)

me = ms.default_context().load_matrix_element("libamplicolmadspace_f.so", "")
me2 = ms.default_context().load_matrix_element("libamplicolmadspace_c.so", "")

m = 10
n = 10000000
tol = 1e-3

psmap = ms.PhaseSpaceMapping([0., 0., 0., 0.], 13000.)

func = ms.MatrixElement(
    me,
    [
        ms.MatrixElement.momenta_in,
        ms.MatrixElement.alpha_s_in,
        ms.MatrixElement.flavor_in,
        ms.MatrixElement.channel_in,
        ms.MatrixElement.random_helicity_in
    ],
    [ms.MatrixElement.matrix_element_out, ms.MatrixElement.helicity_index_out]
)

func2 = ms.MatrixElement(
    me2,
    [
        ms.MatrixElement.momenta_in,
        ms.MatrixElement.alpha_s_in,
        ms.MatrixElement.flavor_in,
        ms.MatrixElement.channel_in,
        ms.MatrixElement.random_helicity_in
    ],
    [ms.MatrixElement.matrix_element_out, ms.MatrixElement.helicity_index_out]
)

print(f"Benchmarking Fortran and C implementations over m = {m} iterations, "
      f"each with n = {n} random phase space points.")

fortran_times = np.zeros(m)
cxx_times = np.zeros(m)


print("Running warm-up iterations to avoid cold-start effects...")

p_ext, x1, x2, det = psmap.map_forward([np.random.rand(n, psmap.random_dim())])
alpha_s = 0.118 * np.ones(n)
flavors = np.zeros(n, dtype=np.int32)
channels = np.zeros(n, dtype=np.int32)
rnd_hel = np.random.rand(n)
amp2, hel = func(p_ext, alpha_s, flavors, channels, rnd_hel)
amp22, hel2 = func2(p_ext, alpha_s, flavors, channels, rnd_hel)

print("Warm-up iterations completed. Starting benchmark...")

for i in range(m):
    p_ext, x1, x2, det = psmap.map_forward([np.random.rand(n, psmap.random_dim())])
    alpha_s = 0.118 * np.ones(n)
    flavors = np.zeros(n, dtype=np.int32)
    channels = np.zeros(n, dtype=np.int32)
    rnd_hel = np.random.rand(n)

    t0 = time.perf_counter()
    amp2, hel = func(p_ext, alpha_s, flavors, channels, rnd_hel)
    t1 = time.perf_counter()
    fortran_times[i] = t1 - t0

    t0 = time.perf_counter()
    amp22, hel2 = func2(p_ext, alpha_s, flavors, channels, rnd_hel)
    t1 = time.perf_counter()
    cxx_times[i] = t1 - t0

    sums = amp2 + amp22
    errors = np.abs(amp2 - amp22)
    rel_errors = np.where(sums != 0, errors / sums, 0.0)
    max_error = np.max(rel_errors)
    if max_error > tol:
        print(f"  Iteration {i}: WARNING - numerical consistency check failed "
              f"(max relative error {max_error:.3e} > {tol:.1e}).")
    else:
        print(f"  Iteration {i}: consistency OK "
              f"(max relative error {max_error:.3e}).")

print()
print(f"Fortran runtime: mean = {np.mean(fortran_times):.6f} s, "
      f"std = {np.std(fortran_times):.6f} s")
print(f"C runtime:       mean = {np.mean(cxx_times):.6f} s, "
      f"std = {np.std(cxx_times):.6f} s")
