import argparse
import madspace as ms
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("filename", nargs="?", default="events_output.txt")
args = ap.parse_args()

np.random.seed(12345)

print("Loading C++...")
me = ms.default_context().load_matrix_element("libamplicolmadspace_c.so", "")
n_part = me.particle_count()
n = 1000000

psmap = ms.PhaseSpaceMapping(np.zeros(n_part), 13000., leptonic=True)
res = psmap.map_forward([np.random.rand(n, psmap.random_dim())])
p_ext = res.momenta
alpha_s = 0.118 * np.ones(n)
flavors = np.zeros(n, dtype=np.int32)
channels = np.zeros(n, dtype=np.int32)
rnd_hel = np.random.rand(n)

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
amp2, hel = func(p_ext, alpha_s, flavors, channels, rnd_hel)

with open(args.filename, "w") as f:
    for i in range(n):
        f.write(f"Event #{i + 1}\n")
        f.write("-" * 79 + "\n")
        f.write(f"{' n':>3}{'E':>26}{'px':>26}{'py':>26}{'pz':>26}\n")
        ev = p_ext[i]
        for j in range(n_part):
            f.write(f"{j + 1:>3}{ev[j, 0]:>26.16e}{ev[j, 1]:>26.16e}{ev[j, 2]:>26.16e}{ev[j, 3]:>26.16e}\n")
        f.write("-" * 79 + "\n")
        f.write(f" Matrix element = {amp2[i]:.16e} GeV^-4\n")
        f.write("-" * 79 + "\n")

print(f"Wrote {n} events to {args.filename}")
