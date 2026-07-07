#pragma once
// FeynmanRules.h
// C translation of feynmanrules.f03 / FeynmanRules.hpp.
// Momentum 4-vectors:  p[0..3]  (p[0] = energy)
// Wavefunction arrays: wf[0..3] (size 4), wfT[0..5] (size 6), wfs[0] (size 1)
// Coupling arrays:     coupl[0..1] (size 2)
//
// Dual C/C++ header.  C++ callers get extern "C" linkage automatically;
// AC_CX resolves to std::complex<double> (same ABI as C's double complex).

#include "AmpliColTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

// ─── external wavefunctions ────────────────────────────────────────────────

void ext_gluon_real   (const AC_D_FP p_D[4], int ihel, AC_FP wf[4]);
void ext_gluon_cmplx  (const AC_D_FP p_D[4], int ihel, AC_CX wf[6]);
void ext_vector_mass  (const AC_D_FP p_D[4], int nhel, int nsv,
                       AC_CX wf[6], AC_FP vmass);
void ext_quark        (const AC_D_FP p_D[4], int nhel, AC_CX wf[6], AC_FP fmass);
void ext_antiquark    (const AC_D_FP p_D[4], int nhel, AC_CX wf[6], AC_FP fmass);
void ext_scalar       (const AC_D_FP p_D[4], AC_CX wf[1]);

// ─── vertices ──────────────────────────────────────────────────────────────

void ThreeGluon              (const AC_CX wf1[6], const AC_D_FP pwf1[4],
                              const AC_CX wf2[6], const AC_D_FP pwf2[4],
                              AC_CX wf[6]);
void ThreeGluon_real         (const AC_FP wf1[4], const AC_D_FP pwf1[4],
                              const AC_FP wf2[4], const AC_D_FP pwf2[4],
                              AC_FP wf[4]);
void ThreeGluon_coupl        (const AC_CX wf1[6], const AC_D_FP pwf1[4],
                              const AC_CX wf2[6], const AC_D_FP pwf2[4],
                              AC_CX wf[6], const AC_FP coupl[2]);
void FourGluon               (const AC_CX wf1[6], const AC_CX wf2[6],
                              const AC_CX wf3[6], AC_CX wf[6]);
void TwoGluontoTensor        (const AC_CX wfg1[6], const AC_CX wfg2[6],
                              AC_CX wfT[6]);
void TwoGluontoTensor_real   (const AC_FP wfg1[4], const AC_FP wfg2[4],
                              AC_FP wfT[6]);
void TwoGluontoTensor_coupl  (const AC_CX wfg1[6], const AC_CX wfg2[6],
                              AC_CX wfT[6], const AC_FP coupl[2]);
void TensorGluontoGluon      (const AC_CX wfT1[6], const AC_CX wfg2[6],
                              AC_CX wfg[6]);
void TensorGluontoGluon_real (const AC_FP wfT1[6], const AC_FP wfg2[4],
                              AC_FP wfg[4]);
void TensorGluontoGluon_coupl(const AC_CX wfT1[6], const AC_CX wfg2[6],
                              AC_CX wfg[6], const AC_FP coupl[2]);
void GluonTensortoGluon      (const AC_CX wfg1[6], const AC_CX wfT2[6],
                              AC_CX wfg[6]);
void GluonTensortoGluon_real (const AC_FP wfg1[4], const AC_FP wfT2[6],
                              AC_FP wfg[4]);
void GluonTensortoGluon_coupl(const AC_CX wfg1[6], const AC_CX wfT2[6],
                              AC_CX wfg[6], const AC_FP coupl[2]);

void QuarkGluontoQuark         (const AC_CX wfq1[6], const AC_CX wfg2[6],
                                AC_CX wfq[6]);
void QuarkGluontoQuark_real    (const AC_CX wfq1[6], const AC_FP wfg2[4],
                                AC_CX wfq[6]);
void QuarkGluontoQuark_coupl   (const AC_CX wfq1[6], const AC_CX wfg2[6],
                                AC_CX wfq[6], const AC_FP coupl[2]);
void GluonQuarktoQuark         (const AC_CX wfg1[6], const AC_CX wfq2[6],
                                AC_CX wfq[6]);
void GluonQuarktoQuark_real    (const AC_FP wfg1[4], const AC_CX wfq2[6],
                                AC_CX wfq[6]);
void GluonQuarktoQuark_coupl   (const AC_CX wfg1[6], const AC_CX wfq2[6],
                                AC_CX wfq[6], const AC_FP coupl[2]);
void AQuarkGluontoAQuark       (const AC_CX wfq1[6], const AC_CX wfg2[6],
                                AC_CX wfq[6]);
void AQuarkGluontoAQuark_coupl (const AC_CX wfq1[6], const AC_CX wfg2[6],
                                AC_CX wfq[6], const AC_FP coupl[2]);
void GluonAQuarktoAQuark       (const AC_CX wfg1[6], const AC_CX wfq2[6],
                                AC_CX wfq[6]);
void GluonAQuarktoAQuark_coupl (const AC_CX wfg1[6], const AC_CX wfq2[6],
                                AC_CX wfq[6], const AC_FP coupl[2]);
void QuarkAQuarktoGluon        (const AC_CX wfq1[6], const AC_CX wfq2[6],
                                AC_CX wfg[6], const AC_FP coupl[2]);
void AQuarkQuarktoGluon        (const AC_CX wfq1[6], const AC_CX wfq2[6],
                                AC_CX wfg[6]);
void LeptonALeptontoGluon      (const AC_CX wfq1[6], const AC_CX wfq2[6],
                                AC_CX wfg[6], const AC_FP coupl[2]);
void ALeptonLeptontoGluon      (const AC_CX wfq1[6], const AC_CX wfq2[6],
                                AC_CX wfg[6], const AC_FP coupl[2]);
void QuarkScalartoQuark        (const AC_CX wfq1[6], const AC_CX wfs2[1],
                                AC_CX wfq[6], const AC_FP coupl[2]);
void GluonGluontoScalar        (const AC_CX wfg1[6], const AC_CX wfg2[6],
                                AC_CX wfs[1], const AC_FP coupl[2]);
void ScalarGluontoGluon        (const AC_CX wfs1[1], const AC_CX wfg2[6],
                                AC_CX wfg[6], const AC_FP coupl[2]);
void GluonScalartoGluon        (const AC_CX wfg1[6], const AC_CX wfs2[1],
                                AC_CX wfg[6], const AC_FP coupl[2]);
void ScalarScalartoScalar      (const AC_CX wfs1[1], const AC_CX wfs2[1],
                                AC_CX wfs[1], const AC_FP coupl[2]);

// ─── propagators ───────────────────────────────────────────────────────────

void GluonPropagator      (AC_CX wfg[6],  const AC_D_FP p[4]);
void GluonPropagator_real (AC_FP wfg[4],  const AC_D_FP p[4]);
void GluonPropagator_mass (AC_CX wfg[6],  const AC_D_FP p[4], AC_D_FP vm, AC_D_FP vw);
void QuarkPropagator      (AC_CX wfq[6],  const AC_D_FP p[4], AC_D_FP fm, AC_D_FP fw);
void AQuarkPropagator     (AC_CX wfq[6],  const AC_D_FP p[4], AC_D_FP fm, AC_D_FP fw);
void ScalarPropagator     (AC_CX wfs[1],  const AC_D_FP p[4], AC_D_FP sm, AC_D_FP sw);

#ifdef __cplusplus
} // extern "C"
#endif
