import json

import madspace as ms
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(12345)

me  = ms.default_context().load_matrix_element("libamplicolmadspace_f.so", "")

n_part = me.particle_count()

me2 = ms.default_context().load_matrix_element("libamplicolmadspace_c.so", "")
if me2.particle_count() != n_part:
    raise RuntimeError("C++ and Fortran implementations have different particle counts, check library consistency.")

n = 1000000

with open("madspace.json") as f:
    _madspace_data = json.load(f)
_channel1 = next(c for c in _madspace_data["channels"] if c["channel"] == 1)
_proc0 = _channel1["processes"][0]
_me_flavor0 = next(m for m in _proc0["matrix_elements"] if m["label"] == 1)
pids = _me_flavor0["pdg_ids"]
if len(pids) != n_part:
    raise RuntimeError("Not same particle cound in madspace.json and n_par")

jet_pids = ms.Observable.jet_pids
lepton_pids = ms.Observable.lepton_pids
Opt = ms.Observable.ObservableOption

def _obs(observable, *select_groups):
    return ms.Observable(pids, observable, list(select_groups))

cuts = ms.Cuts([
    ms.CutItem(_obs(Opt.obs_pt, jet_pids), min=20.0),
    ms.CutItem(_obs(Opt.obs_eta_abs, jet_pids), max=5.0),
    ms.CutItem(_obs(Opt.obs_delta_r, jet_pids), min=0.4),                # jet-jet

    ms.CutItem(_obs(Opt.obs_pt, lepton_pids), min=10.0),
    ms.CutItem(_obs(Opt.obs_eta_abs, lepton_pids), max=2.5),
    ms.CutItem(_obs(Opt.obs_delta_r, lepton_pids), min=0.4),             # lepton-lepton

    ms.CutItem(_obs(Opt.obs_delta_r, jet_pids, lepton_pids), min=0.4),   # jet-lepton

    ms.CutItem(ms.Observable(pids, Opt.obs_sqrt_s, []), min=100.0),
])

psmap = ms.PhaseSpaceMapping(
    np.zeros(n_part), 13000.,
    invariant_power=0.8,
    cuts=cuts,
    leptonic=True,
)

p_ext, x1, x2, det = psmap.map_forward([np.random.rand(n, psmap.random_dim())])
alpha_s = 0.118 * np.ones(n)
flavors = np.zeros(n, dtype=np.int32)
channels = np.zeros(n, dtype=np.int32)
rnd_hel = np.random.rand(n)

# print("Fortran:")
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
# print(amp2)
# print(hel)
# print("C++:")
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
amp22, hel2 = func2(p_ext, alpha_s, flavors, channels, rnd_hel)

print(f"Testing Fortran and C++ implementations for consistency, using n = {n} random phase space points.")

sums = amp2 + amp22
errors = np.abs(amp2 - amp22)
rel_errors = np.where(sums != 0, errors / sums, 0.0)
tot_errors = np.sum(rel_errors)
print(f"Total error between Fortran and C++ implementations: {tot_errors:.15e}")
mean_error = np.mean(rel_errors)
print(f"Mean error between Fortran and C++ implementations: {mean_error:.15e}")
median_error = np.median(rel_errors)
print(f"Median error between Fortran and C++ implementations: {median_error:.15e}")
n_pts_median = np.sum(rel_errors <= median_error)
print(f"Fraction of points with error less than or equal to the median: {n_pts_median / n:.3e}")
max_error = np.max(rel_errors)
print(f"Maximum error between Fortran and C++ implementations: {max_error:.15e}")
max_error_ind = np.argmax(rel_errors)
print(f"Max error for event {max_error_ind}")

print(f"Phase space point for that specific event")
print("-" * 79 + "\n")
print(f"{' n':>3}{'E':>26}{'px':>26}{'py':>26}{'pz':>26}\n")
ev = p_ext[max_error_ind]
for j in range(n_part):
    print(f"{j + 1:>3}{ev[j, 0]:>26.16e}{ev[j, 1]:>26.16e}{ev[j, 2]:>26.16e}{ev[j, 3]:>26.16e}\n")
print("-" * 79 + "\n")
print(f" Matrix element F = {amp2[max_error_ind]:.16e} GeV^-4\n")
print("-" * 79 + "\n")
print(f" Matrix element C = {amp22[max_error_ind]:.16e} GeV^-4\n")
print("-" * 79 + "\n")
# print(amp22)
# print(hel2)

eps = 1e-17
digits = -np.log10(np.maximum(rel_errors, eps))

digits_int = np.floor(digits).astype(int)

fig,ax = plt.subplots(figsize=(8, 5))

counts, bins, patches = ax.hist(
    digits_int,
    bins=np.arange(digits_int.min(), digits_int.max() + 2) - 0.5,
    edgecolor="black",
)


ax.set_xlabel(f"Number of agreeing decimal digits         SUm={len(digits)}")
ax.set_ylabel("Number of phase-space points")
ax.set_title("Agreement between Fortran and C++ implementations")
ax.set_xticks(np.arange(digits_int.min(), digits_int.max() + 1))
ax.grid(True, alpha=0.3)

# Add counts above each bar
for count, patch in zip(counts, patches):
    if count > 0:
        ax.text(
            patch.get_x() + patch.get_width() / 2,
            count,
            f"{int(count)}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

plt.tight_layout()
plt.show()
plt.savefig("FortranVSC.pdf")
