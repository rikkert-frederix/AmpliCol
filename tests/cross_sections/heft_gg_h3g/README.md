# HEFT `g g > h + gluons` cross-section comparison

This comparison uses complete MadGraph event generation, not its survey-only
estimate.  With the matched inputs below, the required 100,000-event AmpliCol
sample and the 10,000-event MadGraph sample agree within one combined standard
deviation.

## Matched inputs

- proton-proton collisions at 14 TeV;
- the explicit partonic process `g g > h g g g`;
- `NNPDF23_nlo_as_0119_qed`, member 0, LHAPDF ID 244800;
- `mu_R = mu_F = H_T/2`, where `H_T` includes the Higgs transverse mass and
  all three gluon transverse momenta;
- for resolved final-state gluons, `pT(g) > 30 GeV`, `|eta(g)| < 6`, and
  `DeltaR(g,g) > 0.4`, with no dijet mass cut;
- `mH = 125 GeV`, zero external-Higgs width, `kappa = 1`, and the derived
  `sin(theta_W) = 0.4714302554840721` and `v = 246.2184581014718 GeV`;
- strict HEFT in MadGraph's bundled `hgg_plugin`: `MT = 1e9 GeV` removes its
  finite-top expansion, and `Gf = 1.16639000000326e-5` gives the same vev.

Both programs update the Wilson coefficient event by event.  AmpliCol sets
`gH = alpha_s(mu_R)/(3*pi*v)` from LHAPDF before evaluating each amplitude and
stores that `alpha_s` in LHE `AQCDUP`; the full-colour reweighter reads the
event's `AQCDUP`.  MadGraph updates `G` at the dynamical renormalization scale
and then recomputes the model's `GH`, which is proportional to `G**2/v`.

## Reproduction

Generate and reweight the required AmpliCol sample from a scratch directory:

```sh
mkdir -p amplicol_gg_h3g/Outputs
cd amplicol_gg_h3g
/path/to/heft/process_list.py --heft --serial 'g g > h g g g'
/path/to/heft/amplicol_generate \
  --process=processes.txt \
  --input=/path/to/heft/tests/cross_sections/heft_gg_h3g/amplicol_run_card.dat \
  --nevents=100000 --seed=260822 --tag=gg_h3g_100k_seed260822
/path/to/heft/amplicol_reweight \
  Outputs/gg_h3g_100k_seed260822_events.lhe \
  --input=/path/to/heft/tests/cross_sections/heft_gg_h3g/amplicol_run_card.dat
```

Generate the MadGraph output, apply every matched card setting, and generate
10,000 events:

```sh
cd /path/to/heft/tests/cross_sections/heft_gg_h3g
/path/to/MG5_aMC/bin/mg5_aMC madgraph_process.mg5
./configure_madgraph.py MG5_gg_h3g
MG5_gg_h3g/bin/madevent generate_events matched_h3g_events \
  -f --multicore --nb_core=4
```

The MadGraph card helper is intentionally run before `generate_events` and
sets `nevents = 10000`.  A standalone `madevent survey` is insufficient for
this comparison: in the recorded runs its central values were 0.60826 pb with
strategy 1 and 0.60955 pb with strategy 2, while refinement during event
generation produced the stable result below.

Summarize and statistically compare the two generated samples with:

```sh
./compare_rates.py \
  /path/to/amplicol_gg_h3g/Outputs/gg_h3g_100k_seed260822_events.lhe \
  /path/to/amplicol_gg_h3g/Outputs/gg_h3g_100k_seed260822_events.lhe.rwgt \
  MG5_gg_h3g/SubProcesses/results.dat \
  --madgraph-banner \
  MG5_gg_h3g/Events/matched_h3g_events/matched_h3g_events_tag_1_banner.txt
```

The one-body check uses the dedicated rapidity map, with
`x1*x2 = mH^2/S`, and therefore requires PDFs.  Reproduce it with:

```sh
mkdir -p amplicol_gg_h/Outputs
cd amplicol_gg_h
/path/to/heft/process_list.py --heft --serial 'g g > h'
/path/to/heft/amplicol_generate \
  --process=processes.txt \
  --input=/path/to/heft/tests/cross_sections/heft_gg_h3g/amplicol_run_card.dat \
  --nevents=100000 --seed=260824 --tag=gg_h_100k
/path/to/heft/amplicol_reweight \
  Outputs/gg_h_100k_events.lhe \
  --input=/path/to/heft/tests/cross_sections/heft_gg_h3g/amplicol_run_card.dat

cd /path/to/heft/tests/cross_sections/heft_gg_h3g
/path/to/MG5_aMC/bin/mg5_aMC madgraph_onebody_process.mg5
./configure_madgraph.py MG5_gg_h \
  --settings madgraph_onebody_settings.txt
MG5_gg_h/bin/madevent generate_events matched_h_events \
  -f --multicore --nb_core=4
```

The separate one-body settings file omits jet-cut keys because MadGraph does
not put them in a `2->1` run card; all applicable beam, PDF, scale, mass,
coupling, width, helicity, integration, and seed settings are still updated
before event generation.

## Recorded results

| Process | AmpliCol events | AmpliCol full colour [pb] | MadGraph events | MadGraph [pb] | Pull |
|---|---:|---:|---:|---:|---:|
| `g g > h` | 100,000 | 17.528460 +/- 0.004722 | 10,000 | 17.5200 +/- 0.014932 | 0.54 |
| `g g > h g` | 100,000 | 7.1688290 +/- 0.0037256 | 10,000 | 7.1521592 +/- 0.013948 | 1.15 |
| `g g > h g g` | 100,000 | 2.2226501 +/- 0.0017198 | 10,000 | 2.2220480 +/- 0.010004 | 0.06 |
| `g g > h g g g` | 100,000 | 0.6164148 +/- 0.0005207 | 10,000 | 0.6145901 +/- 0.0019338 | 0.91 |

The primary AmpliCol leading-colour result is
`0.69334809 +/- 0.00058566 pb`; full-colour reweighting changes every event by
`8/9` for this pure five-gluon colour structure.  An independent second
100,000-event seed gave `0.6167573 +/- 0.0006849 pb`.  Their inverse-variance
combination is `0.6165402 +/- 0.0004145 pb`, also consistent with MadGraph
(`0.99 sigma`).  Large LHE files are retained as run artifacts rather than
committed to the repository.  The samples were generated with the earlier
independent value `v = 246.21965 GeV`; the recorded rates have been normalized
to the now-derived vev by the exact common factor
`(246.21965/246.2184581014718)^2 = 1.00000968166`.
