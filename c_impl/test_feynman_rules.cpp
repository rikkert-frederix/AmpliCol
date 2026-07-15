// test_feynman_rules.cpp
//
// Two-mode test suite for the FeynmanRules C library.
//
// MODE 1 – Self-consistency (always runs):
//   Tests algebraic properties of the library that must hold
//   independently of the Fortran: symmetry relations, propagator
//   inverses, variant equivalences (real/complex/coupl), massless
//   limits, etc.
//
// MODE 2 – Cross-language diff (runs when reference file is present):
//   Reads a binary reference file produced by the companion Fortran
//   driver (test_feynman_rules_reference.f03) and compares every
//   output value to the C result.  Invoke the Fortran driver first:
//
//     gfortran -O2 test_feynman_rules_reference.f03 feynmanrules.f03
//              -o fortran_ref && ./fortran_ref
//   then run this binary – it will find "fortran_reference.dat" automatically.
//
// Build:
//   g++ -std=c++17 -O2 -Wall test_feynman_rules.cpp FeynmanRules.o -o test_feynman_rules
// Run:
//   ./test_feynman_rules

#include "FeynmanRules.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// ─── imaginary unit constant for test expressions ──────────────────────────
static const AC_CX cImag(0.0, 1.0);

// ─── tiny test framework ───────────────────────────────────────────────────

static int g_pass = 0, g_fail = 0;
static std::string g_suite;

void suite(const std::string& name) {
    g_suite = name;
    std::cout << "\n=== " << name << " ===\n";
}

void check(bool ok, const std::string& label, AC_FP got = 0, AC_FP tol = 0) {
    if (ok) {
        ++g_pass;
        std::cout << "  PASS  " << label << "\n";
    } else {
        ++g_fail;
        std::cout << "  FAIL  " << label;
        if (tol > 0) std::cout << "  (|err|=" << std::scientific << got << ", tol=" << tol << ")";
        std::cout << "\n";
    }
}

// Check |a - b| / (|b| + 1) < tol  (relative+absolute mixed norm)
bool near(AC_CX a, AC_CX b, AC_FP tol = 1e-12) {
    AC_FP err   = std::abs(a - b);
    AC_FP scale = std::abs(b) + 1.0;
    return err / scale < tol;
}
bool near(AC_FP a, AC_FP b, AC_FP tol = 1e-12) {
    return std::abs(a - b) / (std::abs(b) + 1.0) < tol;
}

// Max mixed-norm error over n complex elements
static AC_FP maxErr_cx(const AC_CX* a, const AC_CX* b, int n) {
    AC_FP e = 0;
    for (int i = 0; i < n; ++i)
        e = std::max(e, (AC_FP)(std::abs(a[i] - b[i]) / (std::abs(b[i]) + 1.0)));
    return e;
}
// Max mixed-norm error over n real elements
static AC_FP maxErr_fp(const AC_FP* a, const AC_FP* b, int n) {
    AC_FP e = 0;
    for (int i = 0; i < n; ++i)
        e = std::max(e, (AC_FP)(std::abs(a[i] - b[i]) / (std::abs(b[i]) + 1.0)));
    return e;
}

// ─── reference momenta ─────────────────────────────────────────────────────

static const AC_D_FP pGeneric[]   = {10.0, 3.0, 4.0, 5.0};    // massive, off-shell
static const AC_D_FP pMassless[]  = {10.0, 6.0, 0.0, 8.0};    // p^2 = 0
static const AC_D_FP pAlongZ[]    = { 5.0, 0.0, 0.0, 5.0};    // pt=0, p[3]>0
static const AC_D_FP pAlongNZ[]   = { 5.0, 0.0, 0.0,-5.0};    // pt=0, p[3]<0
static const AC_D_FP pIncoming[]  = {-std::sqrt(50.0), 3.0, 4.0, 5.0};    // p[0]<0, on-shell massless (E^2=px^2+py^2+pz^2=50)
static const AC_D_FP pIncomingZ[] = {-5.0, 0.0, 0.0, 5.0};    // p[0]<0, pt=0
static const AC_D_FP pMass2[]     = { 7.0, 1.0, 2.0, 3.0};    // second momentum

// ─── helpers for symmetry checks ───────────────────────────────────────────

static AC_CX minkDot(const AC_CX* a, const AC_CX* b) {
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
}

// ============================================================================
// SECTION 1: external wavefunctions
// ============================================================================

void test_ext_gluon_cmplx_transverse() {
    suite("ext_gluon_cmplx – transversality & masslessness");
    for (int hel : {1, -1}) {
        AC_CX wf[6] = {};
        ext_gluon_cmplx(pMassless, hel, wf);
        AC_CX dot = wf[0]*(AC_FP)pMassless[0] - wf[1]*(AC_FP)pMassless[1]
                  - wf[2]*(AC_FP)pMassless[2] - wf[3]*(AC_FP)pMassless[3];
        std::string lab = "transversality hel=" + std::to_string(hel);
        AC_FP err = std::abs(dot) / (std::abs(wf[0]) + 1.0);
        check(err < 1e-12, lab, err, 1e-12);
    }
    // ε(+) · ε(-) = +1  (HELAS normalisation with metric (+,-,-,-))
    AC_CX wfP[6] = {}, wfM[6] = {};
    ext_gluon_cmplx(pMassless, +1, wfP);
    ext_gluon_cmplx(pMassless, -1, wfM);
    AC_CX dot = minkDot(wfP, wfM);
    check(near(dot, AC_CX(1.0, 0.0), 1e-12), "ε(+)·ε(-) = +1",
          std::abs(dot - AC_CX(1.0, 0.0)), 1e-12);
}

void test_ext_gluon_real_vs_cmplx() {
    suite("ext_gluon_real – equivalence to linear combination of ±1 helicities");
    const AC_FP sqh = std::sqrt(0.5);
    AC_CX wfP[6] = {}, wfM[6] = {};
    ext_gluon_cmplx(pMassless, +1, wfP);
    ext_gluon_cmplx(pMassless, -1, wfM);

    AC_FP wfR[4] = {}, wfRm[4] = {};
    ext_gluon_real(pMassless, +1, wfR);
    ext_gluon_real(pMassless, -1, wfRm);

    for (int i = 0; i < 4; ++i) {
        AC_FP expect_p  = (cImag * (wfP[i] + wfM[i])).real() * sqh;
        AC_FP expect_m  = -(wfP[i] - wfM[i]).real() * sqh;
        check(near(wfR[i],  expect_p, 1e-14), "real ihel=+1 component " + std::to_string(i));
        check(near(wfRm[i], expect_m, 1e-14), "real ihel=-1 component " + std::to_string(i));
    }
}

void test_ext_vector_mass_massless_limit() {
    suite("ext_vector_mass – massless limit matches ext_gluon_cmplx");
    for (int hel : {1, -1}) {
        AC_CX wfMass[6] = {}, wfCmplx[6] = {};
        ext_vector_mass(pMassless, hel, +1, wfMass, 0.0);
        ext_gluon_cmplx(pMassless, hel, wfCmplx);
        AC_FP err = maxErr_cx(wfMass, wfCmplx, 6);
        check(err < 1e-12,
              "massless limit hel=" + std::to_string(hel) + " matches cmplx", err, 1e-12);
    }
}

void test_ext_vector_mass_transverse() {
    suite("ext_vector_mass – transversality for massive vector");
    const AC_D_FP pmassive[] = {10.0, 3.0, 4.0, 0.0};
    AC_FP vmass = 5.0;
    for (int hel : {1, -1, 0}) {
        AC_CX wf[6] = {};
        ext_vector_mass(pmassive, hel, hel, wf, vmass);
        AC_FP norm = 0;
        for (int i = 0; i < 4; ++i) norm += std::norm(wf[i]);
        check(norm > 1e-20, "massive wf non-zero hel=" + std::to_string(hel));
    }
}

void test_ext_quark_antiquark_consistency() {
    suite("ext_quark / ext_antiquark – massless helicity spinors");
    for (AC_FP fmass : {0.0, 1.5}) {
        AC_CX wfQ1[6] = {}, wfQ2[6] = {};
        ext_quark(pMassless, +1, wfQ1, fmass);
        ext_quark(pMassless, -1, wfQ2, fmass);
        if (fmass == 0.0) {
            AC_CX ip(0);
            for (int i = 0; i < 4; ++i) ip += std::conj(wfQ1[i]) * wfQ2[i];
            check(std::abs(ip) < 1e-12, "massless quark ortho hel+1 vs -1",
                  std::abs(ip), 1e-12);
        }
        AC_FP n1 = 0, n2 = 0;
        for (int i = 0; i < 4; ++i) { n1 += std::norm(wfQ1[i]); n2 += std::norm(wfQ2[i]); }
        check(n1 > 1e-20, "quark wf non-zero hel=+1 fmass=" + std::to_string((int)fmass));
        check(n2 > 1e-20, "quark wf non-zero hel=-1 fmass=" + std::to_string((int)fmass));
    }
}

void test_ext_quark_incoming_vs_outgoing() {
    suite("ext_quark – incoming (p[0]<0) branch executes without NaN");
    AC_CX wf[6] = {};
    ext_quark(pIncoming, +1, wf, 0.0);
    bool ok = true;
    for (int i = 0; i < 4; ++i) ok = ok && std::isfinite(wf[i].real()) && std::isfinite(wf[i].imag());
    check(ok, "no NaN/Inf for pIncoming massless");

    ext_quark(pIncoming, -1, wf, 2.0);
    ok = true;
    for (int i = 0; i < 4; ++i) ok = ok && std::isfinite(wf[i].real()) && std::isfinite(wf[i].imag());
    check(ok, "no NaN/Inf for pIncoming massive");

    ext_quark(pIncomingZ, +1, wf, 0.0);
    ok = true;
    for (int i = 0; i < 4; ++i) ok = ok && std::isfinite(wf[i].real()) && std::isfinite(wf[i].imag());
    check(ok, "no NaN/Inf for pIncomingZ (pt=0) massless");
}

void test_ext_gluon_along_z_axis() {
    suite("ext_gluon_cmplx – along ±z axis (pt=0 branch)");
    for (const AC_D_FP *p : {pAlongZ, pAlongNZ}) {
        for (int hel : {1, -1}) {
            AC_CX wf[6] = {};
            ext_gluon_cmplx(p, hel, wf);
            bool ok = true;
            for (int i = 0; i < 4; ++i)
                ok = ok && std::isfinite(wf[i].real()) && std::isfinite(wf[i].imag());
            check(ok, std::string("finite pt=0 p[3]=") +
                      std::to_string((int)p[3]) + " hel=" + std::to_string(hel));
        }
    }
}

void test_ext_scalar() {
    suite("ext_scalar");
    AC_CX wf[1] = {};
    ext_scalar(pGeneric, wf);
    check(near(wf[0], AC_CX(1.0, 0.0), 1e-15), "wf[0] == 1+0i");
}

// ============================================================================
// SECTION 2: Vertices – internal consistency
// ============================================================================

void test_three_gluon_antisymmetry() {
    suite("ThreeGluon – antisymmetry under 1↔2 swap");
    AC_CX wf1[6] = {}, wf2[6] = {}, wfOut1[6] = {}, wfOut2[6] = {};
    ext_gluon_cmplx(pGeneric, +1, wf1);
    ext_gluon_cmplx(pMass2,   -1, wf2);
    ThreeGluon(wf1, pGeneric, wf2, pMass2,   wfOut1);
    ThreeGluon(wf2, pMass2,   wf1, pGeneric, wfOut2);
    for (int i = 0; i < 4; ++i) {
        AC_FP err = std::abs(wfOut1[i] + wfOut2[i]) / (std::abs(wfOut1[i]) + 1.0);
        check(err < 1e-12, "antisymmetry component " + std::to_string(i), err, 1e-12);
    }
}

void test_three_gluon_real_vs_cmplx_real_input() {
    suite("ThreeGluon_real – real gluon wf consistency with ThreeGluon");
    AC_FP wfR1[4] = {}, wfR2[4] = {}, wfROut[4] = {};
    ext_gluon_real(pMassless, +1, wfR1);
    ext_gluon_real(pMassless, -1, wfR2);
    ThreeGluon_real(wfR1, pMassless, wfR2, pAlongZ, wfROut);

    AC_CX wfC1[6] = {}, wfC2[6] = {}, wfCOut[6] = {};
    for (int i = 0; i < 4; ++i) { wfC1[i] = wfR1[i]; wfC2[i] = wfR2[i]; }
    ThreeGluon(wfC1, pMassless, wfC2, pAlongZ, wfCOut);
    // wfCOut = i/sqrt(2)*(...),  wfROut = 1/sqrt(2)*(...) → wfCOut = i * wfROut
    for (int i = 0; i < 4; ++i) {
        AC_FP err = std::abs(wfCOut[i] - cImag * AC_CX(wfROut[i])) / (std::abs(wfCOut[i]) + 1.0);
        check(err < 1e-12, "real vs cmplx ThreeGluon comp " + std::to_string(i), err, 1e-12);
    }
}

void test_three_gluon_coupl_unit_coupling() {
    suite("ThreeGluon_coupl – unit coupling reproduces ThreeGluon");
    AC_CX wf1[6] = {}, wf2[6] = {}, wfNoC[6] = {}, wfC[6] = {};
    ext_gluon_cmplx(pGeneric,  +1, wf1);
    ext_gluon_cmplx(pMassless, -1, wf2);
    ThreeGluon(wf1, pGeneric, wf2, pMassless, wfNoC);
    const AC_FP coupl[] = {1.0, 0.0};
    ThreeGluon_coupl(wf1, pGeneric, wf2, pMassless, wfC, coupl);
    AC_FP err = maxErr_cx(wfNoC, wfC, 6);
    check(err < 1e-12, "unit coupl matches no-coupl variant", err, 1e-12);
}

void test_four_gluon_bose_symmetry() {
    suite("FourGluon – Bose symmetry 1↔3 gives same amplitude");
    AC_CX wf1[6] = {}, wf2[6] = {}, wf3[6] = {}, wfOut13[6] = {}, wfOut31[6] = {};
    ext_gluon_cmplx(pGeneric,  +1, wf1);
    ext_gluon_cmplx(pMassless, -1, wf2);
    ext_gluon_cmplx(pMass2,    +1, wf3);
    FourGluon(wf1, wf2, wf3, wfOut13);
    FourGluon(wf3, wf2, wf1, wfOut31);
    AC_FP err = maxErr_cx(wfOut13, wfOut31, 6);
    check(err < 1e-12, "1↔3 symmetry", err, 1e-12);
}

void test_tensor_roundtrip() {
    suite("TwoGluontoTensor / TensorGluontoGluon – roundtrip");
    AC_CX wfg1[6] = {}, wfg2[6] = {}, wfg3[6] = {}, wfT[6] = {}, result[6] = {};
    ext_gluon_cmplx(pGeneric,  +1, wfg1);
    ext_gluon_cmplx(pMassless, -1, wfg2);
    ext_gluon_cmplx(pMass2,    +1, wfg3);
    TwoGluontoTensor(wfg1, wfg2, wfT);
    TensorGluontoGluon(wfT, wfg3, result);
    AC_FP norm = 0;
    bool finite = true;
    for (int i = 0; i < 4; ++i) {
        norm += std::norm(result[i]);
        finite = finite && std::isfinite(result[i].real()) && std::isfinite(result[i].imag());
    }
    check(finite, "TwoGluon→Tensor→Gluon finite");
    check(norm > 1e-20, "TwoGluon→Tensor→Gluon non-zero");
}

void test_tensor_antisymmetry() {
    suite("TwoGluontoTensor – antisymmetry wfg1 ↔ wfg2");
    AC_CX wfg1[6] = {}, wfg2[6] = {}, wfT12[6] = {}, wfT21[6] = {};
    ext_gluon_cmplx(pGeneric,  +1, wfg1);
    ext_gluon_cmplx(pMassless, -1, wfg2);
    TwoGluontoTensor(wfg1, wfg2, wfT12);
    TwoGluontoTensor(wfg2, wfg1, wfT21);
    for (int i = 0; i < 6; ++i) {
        AC_FP err = std::abs(wfT12[i] + wfT21[i]) / (std::abs(wfT12[i]) + 1.0);
        check(err < 1e-12, "antisymmetry tensor component " + std::to_string(i), err, 1e-12);
    }
}

void test_tensor_coupl_vs_uncoupl() {
    suite("TwoGluontoTensor_coupl – coupl[0]=1 matches uncoupled");
    AC_CX wfg1[6] = {}, wfg2[6] = {}, wfT[6] = {}, wfTC[6] = {};
    ext_gluon_cmplx(pGeneric,  +1, wfg1);
    ext_gluon_cmplx(pMassless, -1, wfg2);
    TwoGluontoTensor(wfg1, wfg2, wfT);
    const AC_FP coupl[] = {1.0, 0.0};
    TwoGluontoTensor_coupl(wfg1, wfg2, wfTC, coupl);
    AC_FP err = maxErr_cx(wfT, wfTC, 6);
    check(err < 1e-12, "TwoGluontoTensor coupl=1 matches uncoupled", err, 1e-12);
}

void test_tensor_real_vs_cmplx() {
    suite("TwoGluontoTensor_real – real inputs match complex variant");
    AC_FP wfg1r[4] = {}, wfg2r[4] = {}, wfTr[6] = {};
    ext_gluon_real(pMassless, +1, wfg1r);
    ext_gluon_real(pMassless, -1, wfg2r);
    TwoGluontoTensor_real(wfg1r, wfg2r, wfTr);

    AC_CX wfg1c[6] = {}, wfg2c[6] = {}, wfTc[6] = {};
    for (int i = 0; i < 4; ++i) { wfg1c[i] = wfg1r[i]; wfg2c[i] = wfg2r[i]; }
    TwoGluontoTensor(wfg1c, wfg2c, wfTc);
    for (int i = 0; i < 6; ++i) {
        AC_FP err = std::abs(wfTc[i] - AC_CX(wfTr[i])) / (std::abs(wfTc[i]) + 1.0);
        check(err < 1e-12, "real vs cmplx tensor comp " + std::to_string(i), err, 1e-12);
    }
}

// ─── Quark/gluon vertex families ───────────────────────────────────────────

void test_quark_gluon_coupl_vs_uncoupl() {
    suite("QuarkGluontoQuark – coupl variants vs. uncoupled (symmetric coupling)");
    AC_CX wfq[6] = {}, wfg[6] = {}, wfOut[6] = {}, wfOutC[6] = {};
    ext_quark(pGeneric,  +1, wfq, 0.0);
    ext_gluon_cmplx(pMassless, -1, wfg);
    QuarkGluontoQuark(wfq, wfg, wfOut);
    const AC_FP coupl[] = {1.0, 1.0};
    QuarkGluontoQuark_coupl(wfq, wfg, wfOutC, coupl);
    AC_FP err = maxErr_cx(wfOut, wfOutC, 6);
    check(err < 1e-12, "coupl={1,1} matches uncoupled", err, 1e-12);
}

void test_gluon_quark_symmetry() {
    suite("GluonQuarktoQuark – consistent with QuarkGluontoQuark (Bose)");
    AC_CX wfq[6] = {}, wfg[6] = {}, wfOut1[6] = {}, wfOut2[6] = {};
    ext_quark(pGeneric,  +1, wfq, 0.0);
    ext_gluon_cmplx(pMassless, -1, wfg);
    QuarkGluontoQuark(wfq, wfg, wfOut1);
    GluonQuarktoQuark(wfg, wfq, wfOut2);
    AC_FP err = maxErr_cx(wfOut1, wfOut2, 6);
    check(err < 1e-12, "QGQ == GQQ Lorentz structure", err, 1e-12);
}

void test_gluon_quark_coupl_vs_uncoupl() {
    suite("GluonQuarktoQuark_coupl – coupl={1,1} matches uncoupled");
    AC_CX wfq[6] = {}, wfg[6] = {}, wfOut[6] = {}, wfOutC[6] = {};
    ext_quark(pGeneric,  -1, wfq, 0.0);
    ext_gluon_cmplx(pMassless, +1, wfg);
    GluonQuarktoQuark(wfg, wfq, wfOut);
    const AC_FP coupl[] = {1.0, 1.0};
    GluonQuarktoQuark_coupl(wfg, wfq, wfOutC, coupl);
    AC_FP err = maxErr_cx(wfOut, wfOutC, 6);
    check(err < 1e-12, "coupl={1,1} matches uncoupled", err, 1e-12);
}

void test_aquark_gluon_coupl_vs_uncoupl() {
    suite("AQuarkGluontoAQuark – coupl={1,1} matches uncoupled");
    AC_CX wfq[6] = {}, wfg[6] = {}, wfOut[6] = {}, wfOutC[6] = {};
    ext_antiquark(pGeneric, +1, wfq, 0.0);
    ext_gluon_cmplx(pMassless, -1, wfg);
    AQuarkGluontoAQuark(wfq, wfg, wfOut);
    const AC_FP coupl[] = {1.0, 1.0};
    AQuarkGluontoAQuark_coupl(wfq, wfg, wfOutC, coupl);
    AC_FP err = maxErr_cx(wfOut, wfOutC, 6);
    check(err < 1e-12, "coupl={1,1} matches uncoupled", err, 1e-12);
}

void test_gluon_aquark_symmetry() {
    suite("GluonAQuarktoAQuark – consistent with AQuarkGluontoAQuark");
    AC_CX wfq[6] = {}, wfg[6] = {}, wfOut1[6] = {}, wfOut2[6] = {};
    ext_antiquark(pGeneric, +1, wfq, 0.0);
    ext_gluon_cmplx(pMassless, -1, wfg);
    AQuarkGluontoAQuark(wfq, wfg, wfOut1);
    GluonAQuarktoAQuark(wfg, wfq, wfOut2);
    AC_FP err = maxErr_cx(wfOut1, wfOut2, 6);
    check(err < 1e-12, "AQG == GAQ Lorentz structure", err, 1e-12);
}

void test_quark_gluon_real_vs_cmplx() {
    suite("QuarkGluontoQuark_real – real gluon gives consistent result");
    AC_FP wfgR[4] = {};
    ext_gluon_real(pMassless, +1, wfgR);
    AC_CX wfq[6] = {}, wfOutR[6] = {}, wfOutC[6] = {}, wfgC[6] = {};
    for (int i = 0; i < 4; ++i) wfgC[i] = wfgR[i];
    ext_quark(pGeneric, +1, wfq, 0.0);
    QuarkGluontoQuark_real(wfq, wfgR, wfOutR);
    QuarkGluontoQuark(wfq, wfgC, wfOutC);
    AC_FP err = maxErr_cx(wfOutR, wfOutC, 6);
    check(err < 1e-12, "real gluon variant matches complex", err, 1e-12);
}

void test_gluon_quark_real_vs_cmplx() {
    suite("GluonQuarktoQuark_real – real gluon gives consistent result");
    AC_FP wfgR[4] = {};
    ext_gluon_real(pMassless, -1, wfgR);
    AC_CX wfq[6] = {}, wfOutR[6] = {}, wfOutC[6] = {}, wfgC[6] = {};
    for (int i = 0; i < 4; ++i) wfgC[i] = wfgR[i];
    ext_quark(pGeneric, -1, wfq, 0.0);
    GluonQuarktoQuark_real(wfgR, wfq, wfOutR);
    GluonQuarktoQuark(wfgC, wfq, wfOutC);
    AC_FP err = maxErr_cx(wfOutR, wfOutC, 6);
    check(err < 1e-12, "real gluon variant matches complex", err, 1e-12);
}

void test_scalar_vertices() {
    suite("Scalar vertices – GluonGluon→Scalar roundtrip with ScalarGluon→Gluon");
    AC_CX wfg1[6] = {}, wfg2[6] = {}, wfs[1] = {};
    const AC_FP coupl[] = {1.0, 0.0};
    ext_gluon_cmplx(pGeneric,  +1, wfg1);
    ext_gluon_cmplx(pMassless, -1, wfg2);
    GluonGluontoScalar(wfg1, wfg2, wfs, coupl);
    bool finite = std::isfinite(wfs[0].real()) && std::isfinite(wfs[0].imag());
    check(finite, "GGS finite");
    check(std::abs(wfs[0]) > 1e-20, "GGS non-zero");

    AC_CX wfgOut[6] = {};
    ScalarGluontoGluon(wfs, wfg2, wfgOut, coupl);
    AC_FP norm = 0;
    for (int i = 0; i < 4; ++i) norm += std::norm(wfgOut[i]);
    check(norm > 1e-20, "SGG non-zero");

    AC_CX wfgOut2[6] = {};
    GluonScalartoGluon(wfg2, wfs, wfgOut2, coupl);
    AC_FP err = maxErr_cx(wfgOut, wfgOut2, 6);
    check(err < 1e-12, "SGG == GSG symmetry", err, 1e-12);
}

void test_scalar_scalar_scalar() {
    suite("ScalarScalartoScalar – coupl(2) flag");
    AC_CX wfs1[] = {AC_CX(1.5, 0.3)};
    AC_CX wfs2[] = {AC_CX(0.7, -0.2)};
    AC_CX wfsOut1[1] = {}, wfsOut2[1] = {};
    const AC_FP coupl_real[] = {2.0, 0.0};
    const AC_FP coupl_imag[] = {2.0, -10.0};  // triggers TMP = i
    ScalarScalartoScalar(wfs1, wfs2, wfsOut1, coupl_real);
    ScalarScalartoScalar(wfs1, wfs2, wfsOut2, coupl_imag);
    // coupl(2) = -10 → result2 = i * result1
    AC_FP err = std::abs(wfsOut2[0] - cImag * wfsOut1[0]) /
                (std::abs(wfsOut1[0]) + 1.0);
    check(err < 1e-14, "coupl(2)=-10 gives extra factor i", err, 1e-14);
}

void test_lepton_vertices_symmetry() {
    suite("LeptonAleptontoGluon / AleptonLeptontoGluon – argument-swap antisymmetry");
    AC_CX wfq1[6] = {}, wfq2[6] = {}, wfg1[6] = {}, wfg2[6] = {};
    const AC_FP coupl[] = {1.0, 0.5};
    ext_quark(pGeneric,  +1, wfq1, 0.0);
    ext_antiquark(pMass2, -1, wfq2, 0.0);
    LeptonALeptontoGluon(wfq1, wfq2, wfg1, coupl);
    ALeptonLeptontoGluon(wfq2, wfq1, wfg2, coupl);
    AC_FP err = maxErr_cx(wfg1, wfg2, 6);
    check(err < 1e-12, "LAG == ALG with swapped args", err, 1e-12);
}

void test_quark_antiquark_to_gluon() {
    suite("QuarkAQuarktoGluon / AQuarkQuarktoGluon – finite & non-zero");
    AC_CX wfq[6] = {}, wfaq[6] = {}, wfg1[6] = {}, wfg2[6] = {};
    const AC_FP coupl[] = {1.0, 0.0};
    ext_quark(pGeneric,   +1, wfq,  0.0);
    ext_antiquark(pMass2, -1, wfaq, 0.0);
    QuarkAQuarktoGluon(wfq, wfaq, wfg1, coupl);
    AQuarkQuarktoGluon(wfaq, wfq, wfg2);
    bool f1 = true, f2 = true;
    for (int i = 0; i < 4; ++i) {
        f1 = f1 && std::isfinite(wfg1[i].real()) && std::isfinite(wfg1[i].imag());
        f2 = f2 && std::isfinite(wfg2[i].real()) && std::isfinite(wfg2[i].imag());
    }
    check(f1, "QAQ→G finite");
    check(f2, "AQQ→G finite");
}

// ============================================================================
// SECTION 3: Propagators
// ============================================================================

void test_gluon_propagator_massless() {
    suite("GluonPropagator (massless) – multiplies by -i/p²");
    AC_CX wfg_orig[6] = {}, wfg[6] = {};
    ext_gluon_cmplx(pMassless, +1, wfg);
    std::copy(wfg, wfg+6, wfg_orig);

    AC_CX wfg2_orig[6] = {}, wfg2[6] = {};
    ext_gluon_cmplx(pGeneric, +1, wfg2);
    std::copy(wfg2, wfg2+6, wfg2_orig);
    GluonPropagator(wfg2, pGeneric);
    AC_FP p2 = pGeneric[0]*pGeneric[0] - pGeneric[1]*pGeneric[1]
             - pGeneric[2]*pGeneric[2] - pGeneric[3]*pGeneric[3];
    AC_CX factor = -cImag / p2;
    AC_FP err = maxErr_cx(wfg2, wfg2_orig, 6);
    check(err > 1e-14, "propagator modifies wf");
    for (int i = 0; i < 4; ++i) {
        AC_FP e = std::abs(wfg2[i] - factor * wfg2_orig[i]) / (std::abs(wfg2[i]) + 1.0);
        check(e < 1e-13, "propagator factor comp " + std::to_string(i), e, 1e-13);
    }
}

void test_gluon_propagator_real() {
    suite("GluonPropagator_real – multiplies by 1/p²");
    AC_FP wfg2[]      = {1.0, 0.5, -0.3, 0.8};
    AC_FP wfg2_orig[] = {1.0, 0.5, -0.3, 0.8};
    GluonPropagator_real(wfg2, pGeneric);
    AC_FP p2 = pGeneric[0]*pGeneric[0] - pGeneric[1]*pGeneric[1]
             - pGeneric[2]*pGeneric[2] - pGeneric[3]*pGeneric[3];
    AC_FP factor = 1.0 / p2;
    AC_FP err = maxErr_fp(wfg2, wfg2_orig, 4);
    check(err > 1e-14, "real propagator modifies wf");
    for (int i = 0; i < 4; ++i) {
        AC_FP e = std::abs(wfg2[i] - factor * wfg2_orig[i]) / (std::abs(wfg2[i]) + 1.0);
        check(e < 1e-14, "real prop factor comp " + std::to_string(i), e, 1e-14);
    }
}

void test_gluon_propagator_mass() {
    suite("GluonPropagator_mass – zero width reduces to simple form");
    AC_CX wfg_orig[] = {AC_CX(1,0), AC_CX(0.5,0.3), AC_CX(-0.2,0.1), AC_CX(0.4,0)};
    AC_CX wfg[]      = {AC_CX(1,0), AC_CX(0.5,0.3), AC_CX(-0.2,0.1), AC_CX(0.4,0)};
    GluonPropagator_mass(wfg, pGeneric, 91.2, 0.0);
    bool finite = true, changed = false;
    for (int i = 0; i < 4; ++i) {
        finite  = finite  && std::isfinite(wfg[i].real()) && std::isfinite(wfg[i].imag());
        changed = changed || (std::abs(wfg[i] - wfg_orig[i]) > 1e-20);
    }
    check(finite,  "massive gluon prop finite");
    check(changed, "massive gluon prop modifies wf");
}

void test_quark_propagator() {
    suite("QuarkPropagator – finite & modifies wf");
    AC_CX wfq[6] = {}, wfq_orig[6] = {};
    ext_quark(pGeneric, +1, wfq, 0.0);
    std::copy(wfq, wfq+4, wfq_orig);
    QuarkPropagator(wfq, pMass2, 1.5, 0.0);
    bool finite = true, changed = false;
    for (int i = 0; i < 4; ++i) {
        finite  = finite  && std::isfinite(wfq[i].real()) && std::isfinite(wfq[i].imag());
        changed = changed || (std::abs(wfq[i] - wfq_orig[i]) > 1e-20);
    }
    check(finite,  "quark propagator finite");
    check(changed, "quark propagator modifies wf");
}

void test_aquark_propagator() {
    suite("AQuarkPropagator – finite & modifies wf");
    AC_CX wfq[6] = {}, wfq_orig[6] = {};
    ext_antiquark(pGeneric, +1, wfq, 0.0);
    std::copy(wfq, wfq+4, wfq_orig);
    AQuarkPropagator(wfq, pMass2, 1.5, 0.0);
    bool finite = true, changed = false;
    for (int i = 0; i < 4; ++i) {
        finite  = finite  && std::isfinite(wfq[i].real()) && std::isfinite(wfq[i].imag());
        changed = changed || (std::abs(wfq[i] - wfq_orig[i]) > 1e-20);
    }
    check(finite,  "antiquark propagator finite");
    check(changed, "antiquark propagator modifies wf");
}

void test_quark_propagator_scalar_factor() {
    suite("QuarkPropagator – off-diagonal structure");
    AC_CX wfq[6] = {}, wfq_orig[6] = {};
    for (int i = 0; i < 4; ++i) wfq[i] = wfq_orig[i] = AC_CX(1, 0);
    QuarkPropagator(wfq, pGeneric, 0.0, 0.0);
    QuarkPropagator(wfq, pGeneric, 0.0, 0.0);
    bool finite = true;
    for (int i = 0; i < 4; ++i)
        finite = finite && std::isfinite(wfq[i].real()) && std::isfinite(wfq[i].imag());
    check(finite, "AC_FP quark propagator finite");
}

void test_scalar_propagator() {
    suite("ScalarPropagator – scalar factor i/(p²-m²+imΓ)");
    AC_CX wfs[] = {AC_CX(1.0, 0.0)};
    AC_FP sm = 125.0, sw = 0.004;
    ScalarPropagator(wfs, pGeneric, sm, sw);
    AC_FP p2 = pGeneric[0]*pGeneric[0] - pGeneric[1]*pGeneric[1]
             - pGeneric[2]*pGeneric[2] - pGeneric[3]*pGeneric[3];
    AC_CX expected = cImag / AC_CX(p2 - sm*sm, sm*sw);
    AC_FP err = std::abs(wfs[0] - expected) / (std::abs(expected) + 1.0);
    check(err < 1e-14, "scalar propagator factor", err, 1e-14);
}

// ============================================================================
// SECTION 4: Cross-language diff against Fortran reference dump
// ============================================================================

struct FortranRecord {
    char name[64];
    std::vector<AC_D_FP> values;  // reference file always stores real(kind=8) == double
};

static bool read_fortran_reference(const std::string& path,
                                   std::vector<FortranRecord>& records)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    while (f.good()) {
        FortranRecord rec;
        if (!f.read(rec.name, 64)) break;
        rec.name[63] = '\0';
        for (int i = 62; i >= 0 && rec.name[i] == ' '; --i) rec.name[i] = '\0';
        int32_t n = 0;
        if (!f.read(reinterpret_cast<char*>(&n), 4)) break;
        rec.values.resize(n);
        if (!f.read(reinterpret_cast<char*>(rec.values.data()), n * 8)) break;
        records.push_back(rec);
    }
    return !records.empty();
}

static AC_FP compareComplex(const AC_CX* cpp, int n, const AC_D_FP* fortran) {
    AC_FP err = 0;
    for (int i = 0; i < n; ++i) {
        AC_CX ref(fortran[2*i], fortran[2*i+1]);
        err = std::max(err, (AC_FP)(std::abs(cpp[i] - ref) / (std::abs(ref) + 1.0)));
    }
    return err;
}
static AC_FP compareReal(const AC_FP* cpp, int n, const AC_D_FP* fortran) {
    AC_FP err = 0;
    for (int i = 0; i < n; ++i)
        err = std::max(err, (AC_FP)(std::abs(cpp[i] - fortran[i]) / (std::abs(fortran[i]) + 1.0)));
    return err;
}

// ─── bitwise / ULP diagnostics ──────────────────────────────────────────────
//
// A relative-error tolerance of 1e-12 hides single-ULP drift (~1e-16 rel.
// error) that nonetheless accumulates across many chained vertex calls in
// a full amplitude evaluation. These helpers measure exact ULP distance so
// such drift can be localized to the specific FR routine that introduces it.

static int64_t ulpOrdered(double d) {
    int64_t i;
    std::memcpy(&i, &d, sizeof(i));
    return (i < 0) ? (int64_t)(INT64_MIN - i) : i;
}
static int64_t ulpDist(double a, double b) {
    if (a == b) return 0;
    if (!std::isfinite(a) || !std::isfinite(b)) return -1;  // sentinel: incomparable
    int64_t oa = ulpOrdered(a), ob = ulpOrdered(b);
    return std::abs(oa - ob);
}
static int64_t maxUlpComplex(const AC_CX* cpp, int n, const AC_D_FP* fortran) {
    int64_t m = 0;
    for (int i = 0; i < n; ++i) {
        m = std::max(m, ulpDist(cpp[i].real(), fortran[2*i]));
        m = std::max(m, ulpDist(cpp[i].imag(), fortran[2*i+1]));
    }
    return m;
}
static int64_t maxUlpReal(const AC_FP* cpp, int n, const AC_D_FP* fortran) {
    int64_t m = 0;
    for (int i = 0; i < n; ++i)
        m = std::max(m, ulpDist((double)cpp[i], fortran[i]));
    return m;
}

static int g_bitwise_exact = 0, g_bitwise_drift = 0;

void run_cross_language_tests(const std::vector<FortranRecord>& recs) {
    suite("Cross-language diff (Fortran reference)");
    const AC_FP TOL = 1e-12;

    std::unordered_map<std::string, const FortranRecord*> map;
    for (auto& r : recs) map[r.name] = &r;

    auto report_ulp = [&](const std::string& tag, int64_t ulp) {
        if (ulp == 0) { ++g_bitwise_exact; return; }
        ++g_bitwise_drift;
        std::cout << "        [NOT BITWISE IDENTICAL]  " << tag
                   << "  max|Δulp|=" << ulp << "\n";
    };

    auto check_rec = [&](const std::string& tag, const AC_FP* cpp, int n) {
        auto it = map.find(tag);
        if (it == map.end()) { check(false, tag + " [NOT FOUND IN REFERENCE]"); return; }
        if ((int)it->second->values.size() != n) { check(false, tag + " [WRONG SIZE in reference]"); return; }
        AC_FP err = compareReal(cpp, n, it->second->values.data());
        check(err < TOL, tag, err, TOL);
        report_ulp(tag, maxUlpReal(cpp, n, it->second->values.data()));
    };
    auto check_rec_c = [&](const std::string& tag, const AC_CX* cpp, int n) {
        auto it = map.find(tag);
        if (it == map.end()) { check(false, tag + " [NOT FOUND IN REFERENCE]"); return; }
        if ((int)it->second->values.size() != 2*n) { check(false, tag + " [WRONG SIZE in reference]"); return; }
        AC_FP err = compareComplex(cpp, n, it->second->values.data());
        check(err < TOL, tag, err, TOL);
        report_ulp(tag, maxUlpComplex(cpp, n, it->second->values.data()));
    };

    // ext_gluon_cmplx
    for (auto& [p, pname] : std::vector<std::pair<const AC_D_FP*, std::string>>{
            {pMassless, "massless"}, {pGeneric, "generic"},
            {pAlongZ,   "alongZ"},  {pAlongNZ, "alongNZ"}}) {
        for (int hel : {1, -1}) {
            AC_CX wf[6] = {};
            ext_gluon_cmplx(p, hel, wf);
            check_rec_c("ext_gluon_cmplx_" + pname + "_" + std::to_string(hel), wf, 6);
        }
    }

    // ext_gluon_real
    for (int hel : {1, -1}) {
        AC_FP wf[4] = {};
        ext_gluon_real(pMassless, hel, wf);
        check_rec("ext_gluon_real_" + std::to_string(hel), wf, 4);
    }

    // ext_vector_mass
    {
        AC_CX wf[6] = {};
        ext_vector_mass(pGeneric, 1, 1, wf, 5.0);
        check_rec_c("ext_gluon_mass_1", wf, 6);
        ext_vector_mass(pGeneric, -1, -1, wf, 5.0);
        check_rec_c("ext_gluon_mass_m1", wf, 6);
        ext_vector_mass(pGeneric, 0, 0, wf, 5.0);
        check_rec_c("ext_gluon_mass_0", wf, 6);
    }

    // ext_quark  (outgoing & incoming, massless & massive)
    for (int hel : {1, -1}) {
        for (auto& [fm, fmtag] : std::vector<std::pair<AC_FP, std::string>>{{0.0,"0"},{1.5,"15"}}) {
            AC_CX wf[6] = {};
            ext_quark(pGeneric, hel, wf, fm);
            check_rec_c("ext_quark_" + std::to_string(hel) + "_" + fmtag, wf, 6);
            ext_quark(pIncoming, hel, wf, fm);
            check_rec_c("ext_quark_inc_" + std::to_string(hel) + "_" + fmtag, wf, 6);
        }
    }

    // ext_antiquark
    for (int hel : {1, -1}) {
        for (auto& [fm, fmtag] : std::vector<std::pair<AC_FP, std::string>>{{0.0,"0"},{1.5,"15"}}) {
            AC_CX wf[6] = {};
            ext_antiquark(pGeneric, hel, wf, fm);
            check_rec_c("ext_antiquark_" + std::to_string(hel) + "_" + fmtag, wf, 6);
        }
    }

    // Vertex functions
    {
        AC_CX wf1[6] = {}, wf2[6] = {};
        ext_gluon_cmplx(pGeneric,  +1, wf1);
        ext_gluon_cmplx(pMassless, -1, wf2);

        AC_CX wfOut[6] = {};
        ThreeGluon(wf1, pGeneric, wf2, pMassless, wfOut);
        check_rec_c("ThreeGluon", wfOut, 6);

        FourGluon(wf1, wf2, wf1, wfOut);
        check_rec_c("FourGluon", wfOut, 6);

        AC_CX wfT[6] = {};
        TwoGluontoTensor(wf1, wf2, wfT);
        check_rec_c("TwoGluontoTensor", wfT, 6);

        AC_CX wfg3[6] = {};
        ext_gluon_cmplx(pMass2, +1, wfg3);
        TensorGluontoGluon(wfT, wfg3, wfOut);
        check_rec_c("TensorGluontoGluon", wfOut, 6);

        GluonTensortoGluon(wfg3, wfT, wfOut);
        check_rec_c("GluonTensortoGluon", wfOut, 6);

        AC_CX wfq[6] = {}, wfaq[6] = {};
        ext_quark(pGeneric, +1, wfq, 0.0);
        QuarkGluontoQuark(wfq, wf2, wfOut);
        check_rec_c("QuarkGluontoQuark", wfOut, 6);
        GluonQuarktoQuark(wf2, wfq, wfOut);
        check_rec_c("GluonQuarktoQuark", wfOut, 6);

        ext_antiquark(pMass2, -1, wfaq, 0.0);
        AQuarkGluontoAQuark(wfaq, wf2, wfOut);
        check_rec_c("AQuarkGluontoAQuark", wfOut, 6);
        GluonAQuarktoAQuark(wf2, wfaq, wfOut);
        check_rec_c("GluonAQuarktoAQuark", wfOut, 6);

        const AC_FP coupl[] = {0.7, 0.3};
        QuarkAQuarktoGluon(wfq, wfaq, wfOut, coupl);
        check_rec_c("QuarkAQuarktoGluon", wfOut, 6);

        AQuarkQuarktoGluon(wfaq, wfq, wfOut);
        check_rec_c("AQuarkQuarktoGluon", wfOut, 6);

        AC_CX wfs[1] = {};
        GluonGluontoScalar(wf1, wf2, wfs, coupl);
        check_rec_c("GluonGluontoScalar", wfs, 1);

        ScalarGluontoGluon(wfs, wf2, wfOut, coupl);
        check_rec_c("ScalarGluontoGluon", wfOut, 6);

        for (int i = 0; i < 6; ++i) wfOut[i] = {};
        LeptonALeptontoGluon(wfq, wfaq, wfOut, coupl);
        check_rec_c("LeptonALeptontoGluon", wfOut, 6);

        for (int i = 0; i < 6; ++i) wfOut[i] = {};
        ALeptonLeptontoGluon(wfaq, wfq, wfOut, coupl);
        check_rec_c("ALeptonLeptontoGluon", wfOut, 6);

        AC_CX wfsIn[1] = {};
        ext_scalar(pGeneric, wfsIn);

        for (int i = 0; i < 6; ++i) wfOut[i] = {};
        QuarkScalartoQuark(wfq, wfsIn, wfOut, coupl);
        check_rec_c("QuarkScalartoQuark", wfOut, 6);

        for (int i = 0; i < 6; ++i) wfOut[i] = {};
        GluonScalartoGluon(wf1, wfsIn, wfOut, coupl);
        check_rec_c("GluonScalartoGluon", wfOut, 6);

        AC_CX wfsOut[1] = {};
        ScalarScalartoScalar(wfsIn, wfsIn, wfsOut, coupl);
        check_rec_c("ScalarScalartoScalar", wfsOut, 1);

        const AC_FP couplFlag[] = {0.7, -10.0};
        ScalarScalartoScalar(wfsIn, wfsIn, wfsOut, couplFlag);
        check_rec_c("ScalarScalartoScalar_flag", wfsOut, 1);
    }

    // Propagators
    {
        AC_CX wfg[6] = {};
        ext_gluon_cmplx(pGeneric, +1, wfg);
        GluonPropagator(wfg, pGeneric);
        check_rec_c("GluonPropagator", wfg, 6);
    }
    {
        AC_FP wfg[] = {1.0, 0.5, -0.3, 0.8};
        GluonPropagator_real(wfg, pGeneric);
        check_rec("GluonPropagator_real", wfg, 4);
    }
    {
        AC_CX wfg[] = {AC_CX(1,0), AC_CX(0.5,0.3), AC_CX(-0.2,0.1), AC_CX(0.4,0)};
        GluonPropagator_mass(wfg, pGeneric, 91.2, 2.5);
        check_rec_c("GluonPropagator_mass", wfg, 6);
    }
    {
        AC_CX wfq[6] = {};
        ext_quark(pGeneric, +1, wfq, 0.0);
        QuarkPropagator(wfq, pMass2, 1.5, 0.1);
        check_rec_c("QuarkPropagator", wfq, 6);
    }
    {
        AC_CX wfq[6] = {};
        ext_antiquark(pGeneric, +1, wfq, 0.0);
        AQuarkPropagator(wfq, pMass2, 1.5, 0.1);
        check_rec_c("AQuarkPropagator", wfq, 6);
    }
    {
        AC_CX wfs[] = {AC_CX(1.0, 0.0)};
        ScalarPropagator(wfs, pGeneric, 125.0, 0.004);
        check_rec_c("ScalarPropagator", wfs, 1);
    }
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char** argv) {
    std::cout << "FeynmanRules test suite\n";
    std::cout << std::string(60, '=') << "\n";

    test_ext_gluon_cmplx_transverse();
    test_ext_gluon_real_vs_cmplx();
    test_ext_vector_mass_massless_limit();
    test_ext_vector_mass_transverse();
    test_ext_quark_antiquark_consistency();
    test_ext_quark_incoming_vs_outgoing();
    test_ext_gluon_along_z_axis();
    test_ext_scalar();

    test_three_gluon_antisymmetry();
    test_three_gluon_real_vs_cmplx_real_input();
    test_three_gluon_coupl_unit_coupling();
    test_four_gluon_bose_symmetry();
    test_tensor_roundtrip();
    test_tensor_antisymmetry();
    test_tensor_coupl_vs_uncoupl();
    test_tensor_real_vs_cmplx();

    test_quark_gluon_coupl_vs_uncoupl();
    test_gluon_quark_symmetry();
    test_gluon_quark_coupl_vs_uncoupl();
    test_aquark_gluon_coupl_vs_uncoupl();
    test_gluon_aquark_symmetry();
    test_quark_gluon_real_vs_cmplx();
    test_gluon_quark_real_vs_cmplx();
    test_scalar_vertices();
    test_scalar_scalar_scalar();
    test_lepton_vertices_symmetry();
    test_quark_antiquark_to_gluon();

    test_gluon_propagator_massless();
    test_gluon_propagator_real();
    test_gluon_propagator_mass();
    test_quark_propagator();
    test_aquark_propagator();
    test_quark_propagator_scalar_factor();
    test_scalar_propagator();

    std::string refpath = "fortran_reference.dat";
    if (argc > 1) refpath = argv[1];
    std::vector<FortranRecord> recs;
    if (read_fortran_reference(refpath, recs)) {
        std::cout << "\n[Found Fortran reference: " << refpath
                  << " (" << recs.size() << " records)]\n";
        run_cross_language_tests(recs);
    } else {
        std::cout << "\n[No Fortran reference file found at '" << refpath
                  << "' – skipping cross-language diff.]\n"
                  << "[Build the reference with:]\n"
                  << "  gfortran -O2 test_feynman_rules_reference.f03 feynmanrules.f03 "
                  << "-o fortran_ref && ./fortran_ref\n";
    }

    std::cout << "\n" << std::string(60, '=') << "\n";
    int total = g_pass + g_fail;
    std::cout << "Results: " << g_pass << "/" << total << " passed";
    if (g_fail == 0) std::cout << "  ✓ ALL PASSED\n";
    else             std::cout << "  ✗ " << g_fail << " FAILED\n";

    if (g_bitwise_exact + g_bitwise_drift > 0) {
        int bwtotal = g_bitwise_exact + g_bitwise_drift;
        std::cout << "Bitwise:  " << g_bitwise_exact << "/" << bwtotal
                   << " exact (0 ULP)";
        if (g_bitwise_drift == 0) std::cout << "  ✓ ALL BITWISE IDENTICAL\n";
        else                      std::cout << "  ✗ " << g_bitwise_drift << " DRIFTED\n";
    }
    return (g_fail == 0) ? 0 : 1;
}
