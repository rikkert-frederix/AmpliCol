# HEFT MadGraph reference points

The constants in `heft_regression.f03` were generated independently with
MadGraph5_aMC@NLO 3.7.2 and its bundled `hgg_plugin`:

```text
import model sm
add model hgg_plugin
generate g g > g h
add process g g > g g h
add process g g > h g g g
add process u u~ > h g g g g
add process t t~ > g h
output standalone <directory> --force
```

The standalone `matrix(p, hel, ic)` routines were called at fixed momenta and
fixed helicities, without initial-state spin or colour averaging. The resulting
full-colour squared matrix elements are committed directly in the regression,
so running the AmpliCol tests does not require a MadGraph installation.

The common inputs were

```text
g_s = 1.21771578477671971
g_H = 5.24865999342440800e-5
m_H = 125 GeV
m_t = 173 GeV
Gamma_t = 1.4915 GeV
alpha_EW = 1/132.507
```

The `g g > h g g g` cross-section event point instead uses the actual first
event's running couplings,

```text
mu_R = 265.81348 GeV
alpha_s(mu_R) = 0.10293058
g_s = 1.13730550681465803
g_H = 4.43558169905496733e-5
```

For this point the plugin was put in the strict heavy-top limit with
`MT = 1e9 GeV`. The fixed reference used `Gf = 1.16637870752389e-5`, giving
the plugin `v = 246.21965 GeV`; AmpliCol is passed the resulting `g_H`
explicitly, so this Feynman-rule comparison is independent of the runtime vev.

Here positive `g_H` denotes the coefficient in AmpliCol's convention
`L_HEFT = -g_H h G^a_(mu nu) G^(a,mu nu)/4`. At the inputs above, the bundled
plugin's internal parameter is `GH = -5.24865999342440800e-5`; its `Hgg`,
`Hggg`, and `Hgggg` UFO couplings use that signed parameter. The comparison
therefore passes `g_H = -GH` to AmpliCol. Since every pure-HEFT amplitude is
linear in `g_H`, this normalization test is independent of the runtime choice
`g_H(mu_R) = heft_kappa alpha_s(mu_R)/(3 pi heft_vev)`.

The four comparisons cover:

- `g g > g h`: Hgg and Hggg rules, four helicity configurations;
- `g g > g g h`: Hgggg auxiliary-tensor factorization and the complete colour
  interference matrix, four helicity configurations;
- `g g > h g g g`: all 24 five-gluon trace coefficients and their complete
  colour interference matrix, both the all-minus helicity and the sum over all
  32 helicity configurations at an event from the 100k-event rate sample;
- `u u~ > h g g g g`: all 24 open-quark-line colour coefficients and their
  complete interference matrix, eight nonzero helicity configurations plus an
  equal-quark-helicity zero, at a non-planar 2-to-5 phase-space point;
- `t t~ > g h`: simultaneous SM top-Yukawa and HEFT sectors, all eight fermion
  and gluon helicity configurations, fixing their relative phase.

Upstream references: [MG5_aMC v3.7.2](https://github.com/mg5amcnlo/mg5amcnlo/releases/tag/v3.7.2)
and the bundled [hgg_plugin](https://github.com/mg5amcnlo/mg5amcnlo/tree/3.x/models/hgg_plugin).
