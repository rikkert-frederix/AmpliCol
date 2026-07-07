// FeynmanRules_device.cuh
// Device-side declarations of all Feynman rule functions.
//
// Implementation lives in FeynmanRules_device.cu. These are __device__
// functions with external linkage (not static/__forceinline__), so any
// translation unit that calls them must be device-linked against
// FeynmanRules_device.o (nvcc -rdc=true ... -dlink).
//
// CUDA C++ only (nvcc). AC_CX = std::complex<double> throughout.
// Mirrors FeynmanRules.h / FeynmanRules.c, with:
//   • CMPLX(r,i)   → AC_CX(r,i)       (C++ std::complex constructor)
//   • creal(z)     → (z).real()
//   • fprintf/exit → __trap()          (non-recoverable bad input)

#pragma once
#include "AmpliColTypes.h"   // AC_FP, AC_CX

// ─── external wavefunctions ────────────────────────────────────────────────

__device__ void ext_gluon_real   (const AC_D_FP p_D[4], int ihel, AC_FP wf[4]);
__device__ void ext_gluon_cmplx  (const AC_D_FP p_D[4], int ihel, AC_CX wf[6]);
__device__ void ext_vector_mass  (const AC_D_FP p_D[4], int nhel, int nsv,
                                  AC_CX wf[6], AC_FP vmass);
__device__ void ext_quark        (const AC_D_FP p_D[4], int nhel, AC_CX wf[6], AC_FP fmass);
__device__ void ext_antiquark    (const AC_D_FP p_D[4], int nhel, AC_CX wf[6], AC_FP fmass);
__device__ void ext_scalar       (const AC_D_FP p_D[4], AC_CX wf[1]);

// ─── vertices ──────────────────────────────────────────────────────────────

__device__ void ThreeGluon              (const AC_CX wf1[6], const AC_D_FP pwf1[4],
                                         const AC_CX wf2[6], const AC_D_FP pwf2[4],
                                         AC_CX wf[6]);
__device__ void ThreeGluon_real         (const AC_FP wf1[4], const AC_D_FP pwf1[4],
                                         const AC_FP wf2[4], const AC_D_FP pwf2[4],
                                         AC_FP wf[4]);
__device__ void ThreeGluon_coupl        (const AC_CX wf1[6], const AC_D_FP pwf1[4],
                                         const AC_CX wf2[6], const AC_D_FP pwf2[4],
                                         AC_CX wf[6], const AC_FP coupl[2]);
__device__ void FourGluon               (const AC_CX wf1[6], const AC_CX wf2[6],
                                         const AC_CX wf3[6], AC_CX wf[6]);
__device__ void TwoGluontoTensor        (const AC_CX wfg1[6], const AC_CX wfg2[6],
                                         AC_CX wfT[6]);
__device__ void TwoGluontoTensor_real   (const AC_FP wfg1[4], const AC_FP wfg2[4],
                                         AC_FP wfT[6]);
__device__ void TwoGluontoTensor_coupl  (const AC_CX wfg1[6], const AC_CX wfg2[6],
                                         AC_CX wfT[6], const AC_FP coupl[2]);
__device__ void TensorGluontoGluon      (const AC_CX wfT1[6], const AC_CX wfg2[6],
                                         AC_CX wfg[6]);
__device__ void TensorGluontoGluon_real (const AC_FP wfT1[6], const AC_FP wfg2[4],
                                         AC_FP wfg[4]);
__device__ void TensorGluontoGluon_coupl(const AC_CX wfT1[6], const AC_CX wfg2[6],
                                         AC_CX wfg[6], const AC_FP coupl[2]);
__device__ void GluonTensortoGluon      (const AC_CX wfg1[6], const AC_CX wfT2[6],
                                         AC_CX wfg[6]);
__device__ void GluonTensortoGluon_real (const AC_FP wfg1[4], const AC_FP wfT2[6],
                                         AC_FP wfg[4]);
__device__ void GluonTensortoGluon_coupl(const AC_CX wfg1[6], const AC_CX wfT2[6],
                                         AC_CX wfg[6], const AC_FP coupl[2]);

__device__ void QuarkGluontoQuark         (const AC_CX wfq1[6], const AC_CX wfg2[6],
                                           AC_CX wfq[6]);
__device__ void QuarkGluontoQuark_real    (const AC_CX wfq1[6], const AC_FP wfg2[4],
                                           AC_CX wfq[6]);
__device__ void QuarkGluontoQuark_coupl   (const AC_CX wfq1[6], const AC_CX wfg2[6],
                                           AC_CX wfq[6], const AC_FP coupl[2]);
__device__ void GluonQuarktoQuark         (const AC_CX wfg1[6], const AC_CX wfq2[6],
                                           AC_CX wfq[6]);
__device__ void GluonQuarktoQuark_real    (const AC_FP wfg1[4], const AC_CX wfq2[6],
                                           AC_CX wfq[6]);
__device__ void GluonQuarktoQuark_coupl   (const AC_CX wfg1[6], const AC_CX wfq2[6],
                                           AC_CX wfq[6], const AC_FP coupl[2]);
__device__ void AQuarkGluontoAQuark       (const AC_CX wfq1[6], const AC_CX wfg2[6],
                                           AC_CX wfq[6]);
__device__ void AQuarkGluontoAQuark_coupl (const AC_CX wfq1[6], const AC_CX wfg2[6],
                                           AC_CX wfq[6], const AC_FP coupl[2]);
__device__ void GluonAQuarktoAQuark       (const AC_CX wfg1[6], const AC_CX wfq2[6],
                                           AC_CX wfq[6]);
__device__ void GluonAQuarktoAQuark_coupl (const AC_CX wfg1[6], const AC_CX wfq2[6],
                                           AC_CX wfq[6], const AC_FP coupl[2]);
__device__ void QuarkAQuarktoGluon        (const AC_CX wfq1[6], const AC_CX wfq2[6],
                                           AC_CX wfg[6], const AC_FP coupl[2]);
__device__ void AQuarkQuarktoGluon        (const AC_CX wfq1[6], const AC_CX wfq2[6],
                                           AC_CX wfg[6]);
__device__ void LeptonALeptontoGluon      (const AC_CX wfq1[6], const AC_CX wfq2[6],
                                           AC_CX wfg[6], const AC_FP coupl[2]);
__device__ void ALeptonLeptontoGluon      (const AC_CX wfq1[6], const AC_CX wfq2[6],
                                           AC_CX wfg[6], const AC_FP coupl[2]);
__device__ void QuarkScalartoQuark        (const AC_CX wfq1[6], const AC_CX wfs2[1],
                                           AC_CX wfq[6], const AC_FP coupl[2]);
__device__ void GluonGluontoScalar        (const AC_CX wfg1[6], const AC_CX wfg2[6],
                                           AC_CX wfs[1], const AC_FP coupl[2]);
__device__ void ScalarGluontoGluon        (const AC_CX wfs1[1], const AC_CX wfg2[6],
                                           AC_CX wfg[6], const AC_FP coupl[2]);
__device__ void GluonScalartoGluon        (const AC_CX wfg1[6], const AC_CX wfs2[1],
                                           AC_CX wfg[6], const AC_FP coupl[2]);
__device__ void ScalarScalartoScalar      (const AC_CX wfs1[1], const AC_CX wfs2[1],
                                           AC_CX wfs[1], const AC_FP coupl[2]);

// ─── propagators ───────────────────────────────────────────────────────────

__device__ void GluonPropagator      (AC_CX wfg[6],  const AC_D_FP p[4]);
__device__ void GluonPropagator_real (AC_FP wfg[4],  const AC_D_FP p[4]);
__device__ void GluonPropagator_mass (AC_CX wfg[6],  const AC_D_FP p[4], AC_D_FP vm, AC_D_FP vw);
__device__ void QuarkPropagator      (AC_CX wfq[6],  const AC_D_FP p[4], AC_D_FP fm, AC_D_FP fw);
__device__ void AQuarkPropagator     (AC_CX wfq[6],  const AC_D_FP p[4], AC_D_FP fm, AC_D_FP fw);
__device__ void ScalarPropagator     (AC_CX wfs[1],  const AC_D_FP p[4], AC_D_FP sm, AC_D_FP sw);
