from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import dataclass, field

import torch
import torch.nn as nn

import madspace as ms


# =============================================================================
# 1.  Config / cuts
# =============================================================================
@dataclass
class Cuts:
    ptj_min: float = 30.0;  etaj_max: float = 6.0
    pta_min: float = 30.0;  etaa_max: float = 6.0
    ptl_min: float = -1.0;  etal_max: float = -1.0
    drjj_min: float = 0.4;  draa_min: float = 0.4
    drll_min: float = -1.0


_MASS_TABLE = {6: 173.2, 24: 80.4, 23: 91.19, 25: 125.0}   # |pdg| -> mass; else massless


def _pdg_mass(pdg):
    return _MASS_TABLE.get(abs(int(pdg)), 0.0)


@dataclass
class Config:
    sqrts: float = 14000.0
    pdfset: str = "NNPDF23_nlo_as_0119_qed"
    scale_ref: float | None = None         # None -> dynamic H_T/2 (Fortran scale_choice=2);
                                            # a float -> FIXED ren/fact scale (scale_choice=0, e.g. M_Z=91.188)
    t_invariant_power: float = 0.3
    cuts: Cuts = field(default_factory=Cuts)


# =============================================================================
# 2.  Channel loader  (converted format OR original)
# =============================================================================
@dataclass
class Channel:
    index: int
    color_order: list[int]
    init_states: list[tuple[int, int]]
    pdg_final: list[int]
    # one entry per PROCESS (= unique amplitude = the ME `flavor_in`/iint index);
    # each process is a list of its flavor configs: (factor, pdg_a, pdg_b, final).
    processes: list[list[tuple[float, int, int, tuple]]] = field(default_factory=list)


def load_channels(path):
    data = json.load(open(path))
    converted = "color_orders" in data
    channels, n_out = [], None
    for ci, ch in enumerate(data["channels"]):
        pso = ch["phasespace_order"]
        color_order = list(pso) if converted else [c - 1 for c in pso]
        n_out = len(color_order) - 2
        inits, pdg_final, procs = set(), [], []
        for proc in ch["processes"]:
            cfgs = []
            for me in proc["matrix_elements"]:
                pdg = list(data["pdg_ids"][me["pdg_ids"]] if converted else me["pdg_ids"])
                if len(pdg) != len(color_order) or not all(isinstance(v, int) for v in pdg):
                    continue
                factor = float(me.get("factor", 1.0))
                cfgs.append((factor, pdg[0], pdg[1], tuple(pdg[2:])))
                inits.add((pdg[0], pdg[1]))
                if not pdg_final:
                    pdg_final = pdg[2:]
            if cfgs:
                procs.append(cfgs)
        channels.append(Channel(ci, color_order, sorted(inits) or [(21, 21)],
                                pdg_final or [21] * n_out,
                                procs or [[(1.0, 21, 21, tuple([21] * n_out))]]))
    return channels, n_out


# =============================================================================
# 3.  Cuts -> ms.Cuts (applied INSIDE PhaseSpaceMapping)
# =============================================================================
def build_ms_cuts(cfg: Config, ch: Channel):
    c = cfg.cuts
    a, b = ch.init_states[0]
    pids = [a, b, *ch.pdg_final]
    O = ms.Observable
    items = []

    def add(obs, sel, **kw):
        items.append(ms.CutItem(ms.Observable(pids, obs, [sel]), **kw))

    if c.ptj_min > 0:  add(O.obs_pt,      O.jet_pids,    min=c.ptj_min)
    if c.etaj_max > 0: add(O.obs_eta_abs, O.jet_pids,    max=c.etaj_max)
    if c.drjj_min > 0: add(O.obs_delta_r, O.jet_pids,    min=c.drjj_min)
    if c.pta_min > 0:  add(O.obs_pt,      O.photon_pids, min=c.pta_min)
    if c.etaa_max > 0: add(O.obs_eta_abs, O.photon_pids, max=c.etaa_max)
    if c.draa_min > 0: add(O.obs_delta_r, O.photon_pids, min=c.draa_min)
    if c.ptl_min > 0:  add(O.obs_pt,      O.lepton_pids, min=c.ptl_min)
    if c.etal_max > 0: add(O.obs_eta_abs, O.lepton_pids, max=c.etal_max)
    return ms.Cuts(items) if items else None


# =============================================================================
# 4.  PDFs & alpha_s + matrix element  (MadSpace primitives, train.py style)
# =============================================================================
def load_pdf_and_coupling(cfg: Config, ctx):
    """Return (pdf_grid, running_coupling), initialised on `ctx`. lhapdf is used
    only to locate the data files."""
    import lhapdf
    lhapdf.setVerbosity(0)
    base = lhapdf.paths()[0]
    s = cfg.pdfset
    pdf_grid = ms.PdfGrid(os.path.join(base, s, f"{s}_0000.dat"))
    pdf_grid.initialize_globals(ctx)
    ag = ms.AlphaSGrid(os.path.join(base, s, f"{s}.info"))
    run_coupl = ms.RunningCoupling(ag)
    ag.initialize_globals(ctx)
    return pdf_grid, run_coupl


def load_matrix_element(ctx, library, card=""):
    """AmpliCol UMAMI ME with the inputs DifferentialCrossSection expects
    (test_umami.py pattern). alpha_s_in is supplied by dcs internally."""
    api = ctx.load_matrix_element(library, card)
    return ms.MatrixElement(
        api,
        [ms.MatrixElement.momenta_in, ms.MatrixElement.alpha_s_in,
         ms.MatrixElement.flavor_in, ms.MatrixElement.channel_in,
         ms.MatrixElement.random_helicity_in],
        [ms.MatrixElement.matrix_element_out, ms.MatrixElement.helicity_index_out],
    )


# =============================================================================
# 5.  One color-ordered channel:  randoms -> weight
# =============================================================================
class ColorOrderedChannel:
    """Wraps ONE PhaseSpaceMapping (color_ordered) for this colour ordering, plus a
    DifferentialCrossSection turning (momenta, flavor, channel, x1, x2) into
    |A_sigma(flavor)|^2 f_a f_b / (2 s_hat).  If `me_func` is None (self-test) the
    weight is just the cut phase-space volume det.  Always called with torch tensors."""

    def __init__(self, ch: Channel, n_out, cfg: Config,
                 me_func=None, run_coupl=None, pdf_grid=None):
        self.ch = ch
        self.cfg = cfg
        self.n_out = n_out
        masses = [0.0, 0.0] + [_pdg_mass(p) for p in ch.pdg_final]
        self.psmap = ms.PhaseSpaceMapping(
            masses, cfg.sqrts,
            leptonic=False,
            invariant_power=cfg.t_invariant_power,
            mode=ms.PhaseSpaceMapping.color_ordered,
            color_order=ch.color_order,
            cuts=build_ms_cuts(cfg, ch),
        )
        self.random_dim = self.psmap.random_dim()
        self.discrete_dim = self.psmap.discrete_dim()
        self.dim = self.random_dim                    # MadNIS flow dim (continuous only)

        # flatten the channel's flavor configs; each is a dcs pid_option (for f_a f_b)
        # carrying its process index (the ME flavor_in) and its idenCOandMAPfactor.
        self.cfg_proc, self.cfg_factor, pid_options = [], [], []
        for iproc, cfgs in enumerate(ch.processes):
            for (fac, a, b, final) in cfgs:
                self.cfg_proc.append(iproc)
                self.cfg_factor.append(fac)
                pid_options.append([a, b, *final])
        self.n_config = len(pid_options)

        self.dcs = None
        if me_func is not None:
            # scale_choice=2 in common.f03 is H_T/2 == EnergyScale.half_transverse_mass
            # (massless partons: m_T = p_T, so half the transverse mass == H_T/2).
            if cfg.scale_ref is None:
                escale = ms.EnergyScale(2 + n_out, ms.EnergyScale.half_transverse_mass)
            else:
                escale = ms.EnergyScale(2 + n_out, float(cfg.scale_ref))
            self.dcs = ms.DifferentialCrossSection(
                me_func, cfg.sqrts, run_coupl,
                escale,
                pid_options,
                pdf_grid is not None, pdf_grid is not None,
                pdf_grid, pdf_grid,
                False, True,
            )

    def weight(self, r, disc=None):
        """r: continuous randoms (n, random_dim).  disc: optional discrete choices
        (n, discrete_dim).  madspace returns the per-branch Jacobian only (the 2^dd
        solution multiplicity is NOT baked into det), so this function owns it:
          - disc None (we sample uniformly): multiply by 2^dd, which compensates the
            uniform 1/2^dd averaging and recovers the discrete sum;
          - disc supplied by MadNIS (its MixedFlow carries the learned q_disc): no
            factor, since the integrator already owns the discrete measure."""
        dev, dt = r.device, r.dtype
        mult = 1.0
        if self.discrete_dim:
            if disc is None:                          # uniform: apply 2^dd ourselves
                disc = torch.randint(0, 2, (r.shape[0], self.discrete_dim),
                                     device=dev, dtype=torch.int32)
                mult = float(2 ** self.discrete_dim)
            else:                                     # MadNIS owns q_disc: no factor
                disc = disc.to(torch.int32)
            inputs = [r, disc]
        else:
            inputs = [r]
        p, x1, x2, det = self.psmap.map_forward(inputs)
        good = torch.isfinite(det) & (det > 0)        # cut-surviving, non-degenerate
        out = torch.zeros(r.shape[0], dtype=det.dtype, device=dev)

        if self.dcs is None:                          # self-test: cut phase-space volume
            return torch.where(good, det * mult, out)
        if not bool(good.any()):
            return out

        # evaluate the ME ONLY on surviving events (AmpliCol STOPs on a zero-energy parton)
        pg, x1g, x2g, detg = p[good], x1[good], x2[good], det[good]
        ng = x1g.shape[0]
        chan_id = torch.full((ng,), self.ch.index, device=dev, dtype=torch.int32)
        total = torch.zeros(ng, dtype=detg.dtype, device=dev)
        for c in range(self.n_config):
            flavor = torch.full((ng,), self.cfg_proc[c], device=dev, dtype=torch.int32)
            pdf_id = torch.full((ng,), c, device=dev, dtype=torch.int32)
            rnd_hel = torch.rand(ng, device=dev, dtype=dt)   # helicity is summed exactly
            res = self.dcs(pg.contiguous(), flavor, chan_id, rnd_hel, x1g, x2g, pdf_id)
            # dcs replaces matrix_element_out (slot 0) with the differential cross
            # section; take that, drop the (informational) helicity index.
            dxs = res[0] if isinstance(res, (tuple, list)) else res
            total = total + dxs * float(self.cfg_factor[c])

        out[good] = total * detg * mult #torch.nan_to_num(total * detg * mult, nan=0.0, posinf=0.0, neginf=0.0)
        return out


# =============================================================================
# 6.  Per-ordering integration: one independent integrator per colour ordering
# =============================================================================
def integrate_per_channel(channels, n_out, cfg, ctx, me_library=None, param_card="",
                          use_pdf=True, batch_size=512, epochs=2000,
                          n_points=1_000_000, seed=0, vegas_only=False,
                          learn_discrete=False):
    """Build an INDEPENDENT single-channel integrator for each colour ordering and
    sum the results.  The colour-ordered cross section is Sum_sigma integral|A_sigma|^2 --
    genuinely independent integrals, NOT a multichannel (no 1/Sum g_i) -- so each
    ordering gets its own flow (or VEGAS grid), its own VEGAS pre-training, and its
    own error; the total is the sum with errors added in quadrature.

    learn_discrete=True hands the 2->3 discrete solution choices to MadNIS's MixedFlow
    (a learned DiscreteMADE) instead of sampling them uniformly: input_dim grows by
    discrete_dim, the discrete columns ride at the END of x, and weight() drops the
    2^dd multiplicity since the flow now carries q_disc.  Only continuous-dim orderings
    fall back to the plain single-channel flow.

    NOTE: each ordering is integrated with `n_points`, so the TOTAL cost is
    len(channels) * n_points.  To match a combined run of N points, pass n_points=N/nc."""
    from madnis.integrator import (Integrator, Integrand, VegasPreTraining,
                                   stratified_variance_softclip)

    torch.set_default_dtype(torch.float64)
    me_func = None if me_library is None else load_matrix_element(ctx, me_library, param_card)
    pdf_grid = run_coupl = None
    if use_pdf and me_func is not None:
        pdf_grid, run_coupl = load_pdf_and_coupling(cfg, ctx)

    results = []                       # (index, color_order, I_sigma, err_sigma, rsd)
    total, var = 0.0, 0.0
    for c in channels:
        ch = ColorOrderedChannel(c, n_out, cfg, me_func, run_coupl, pdf_grid)
        dd = ch.discrete_dim

        if learn_discrete and dd:
            def fn(x, ch=ch, nd=dd):
                x = x.double()
                disc = x[:, :nd]
                cont = x[:, nd:].contiguous()
                return ch.weight(cont, disc=disc).double()
            integrand = Integrand(fn, input_dim=ch.random_dim + dd, discrete_dims=[2] * dd)
        else:
            def fn(x, ch=ch):          # continuous only; weight() samples disc uniformly
                return ch.weight(x.double()).double()
            integrand = Integrand(fn, input_dim=ch.dim)

        def vegas_callback(status):
            print(f"[VEGAS] Iteration {status.step + 1}: loss={status.variance:.6f}")

        def callback(status):
            if (status.step + 1) % (epochs // 10) == 0:
                print(f"[MadNIS] Iteration {status.step + 1:4d}: loss={status.loss:.6f}")

        flow_dict = dict(layers=3, units=64, bins=10, activation=nn.LeakyReLU)
        integrator = Integrator(
            integrand,
            flow_kwargs=flow_dict,
            batch_size=batch_size, dtype=torch.float64,
            optimizer=lambda p: torch.optim.Adam(p, lr=1e-3),
            scheduler=lambda o: torch.optim.lr_scheduler.CosineAnnealingLR(o, T_max=epochs),
            loss=stratified_variance_softclip)
        vegas = VegasPreTraining(integrator, bins=64, damping=0.7)

        print(f"Running VEGAS {'pre-' if not vegas_only else ''}training")
        vegas.train([5000, 10000, 20000, 40000, 80000, 100_000], callback=vegas_callback)
        if vegas_only:
            Ii, ei = vegas.integrate(n_points)
        else:
            print(f"Initialize integrator for channel {c.index} (color_order={c.color_order})")
            vegas.initialize_integrator()
            print(f"Running MadNIS training")
            integrator.train(epochs, callback=callback)
            Ii, ei = integrator.integrate(n_points)

        # RSD = sigma_w / mu_w = (err/I) * sqrt(N): per-sample relative spread,
        # independent of how many points were drawn -- the sampler-quality metric.
        rsd = (ei / Ii) * (n_points ** 0.5) if Ii else float("nan")
        results.append((c.index, c.color_order, Ii, ei, rsd))
        total += Ii
        var += ei * ei
        rel = f"  ({ei / Ii * 100:.2f}%)" if Ii else ""
        print(f"  [channel {c.index}] color_order={c.color_order}: "
              f"{Ii:.6e} +/- {ei:.2e}{rel}  RSD={rsd:.3f}")
    return total, var ** 0.5, results


# =============================================================================
# 7.  CLI
# =============================================================================
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json", nargs="?", default="madspace_converted.json")
    ap.add_argument("--library", default="libamplicolmadspace.so")
    ap.add_argument("--param-card", default="")
    ap.add_argument("--no-pdf", action="store_true", help="partonic (no PDF) AmpliCol run")
    ap.add_argument("--sqrts", type=float, default=14000.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n", type=int, default=1_000_000,
                    help="integration points per colour ordering (total = nc * n)")
    ap.add_argument("--vegas-only", action="store_true",
                    help="VEGAS grid only, no normalizing flow")
    ap.add_argument("--learn-discrete", action="store_true",
                    help="hand the 2->3 discrete choices to MadNIS (learned q_disc) "
                         "instead of uniform sampling")
    ap.add_argument("--selftest", action="store_true",
                    help="no ME/PDF: cross-check the colour-ordered phase-space sum")
    args = ap.parse_args()

    cfg = Config(sqrts=args.sqrts)
    channels, n_out = load_channels(args.json)
    print(f"Loaded {len(channels)} channels, n_out = {n_out}")
    for c in channels:
        print(f"  channel {c.index}: color_order={c.color_order}  init={c.init_states[:3]}")

    ctx = ms.default_context()
    me_library = None if args.selftest else args.library
    use_pdf = not (args.selftest or args.no_pdf)

    t0 = time.time()
    sigma, err, per_channel = integrate_per_channel(
        channels, n_out, cfg, ctx, me_library, args.param_card, use_pdf,
        seed=args.seed, n_points=args.n, vegas_only=args.vegas_only,
        learn_discrete=args.learn_discrete)
    n_total = len(channels) * args.n
    tag = "per-channel VEGAS-only integral" if args.vegas_only else "per-channel MadNIS integral"

    unit = "" if args.selftest else " pb"
    print("\n==============================================")
    print(f"  {tag}: {sigma:.6e} +/- {err:.2e}{unit}"
          + (f"   ({err / sigma * 100:.2f}%)" if sigma else ""))
    if sigma:
        print(f"   RSD={(err / sigma) * (n_total ** 0.5):.3f}  (N={n_total:,})")
    print(f"  ({time.time() - t0:.1f}s)")
    print("==============================================")

    # per-channel error budget: which ordering dominates the variance?
    tot_var = sum(e * e for *_, e, _ in per_channel)
    print("  per-channel error budget (RSD is sample-count independent):")
    for idx, order, I, e, rsd in sorted(per_channel, key=lambda r: -r[3]):
        share = (e * e / tot_var * 100) if tot_var else 0.0
        print(f"    channel {idx}: I={I:.4e}  err={e:.2e}  RSD={rsd:7.3f}  "
              f"({share:5.1f}% of variance)  color_order={order}")
    if args.selftest:
        print("[selftest] |A|^2=1, no PDF: each channel integrates the cut "
              "phase-space volume; all orderings must agree.")


if __name__ == "__main__":
    main()
