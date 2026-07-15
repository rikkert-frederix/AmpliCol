// FeynmanRules.c
// C translation of feynmanrules.f03.
//
// All Fortran real(kind=8)    → AC_FP  (double)
// All Fortran complex(kind=8) → AC_CX  (double complex)

#include "FeynmanRules.h"
#include "AmpliColTypes.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ─── helpers ───────────────────────────────────────────────────────────────

static const AC_FP rZero  = 0.0;
static const AC_FP rHalf  = 0.5;
static const AC_FP rOne   = 1.0;
static const AC_FP rTwo   = 2.0;
static const AC_FP sqh    = 0.707106781186547524401;
static const AC_FP isqTwo = 1.0 / 1.414213562373095048801688724209;  // kept for bitwise-identical Fortran results
static const AC_FP tiny   = 1e-8;

// Complex constants via CMPLX (C11 constant-expression macro)
#define cImag   CMPLX(0.0, 1.0)
#define cZero   CMPLX(0.0, 0.0)
#define imHalf  CMPLX(0.0, 0.5)
// imSqh uses sqh (= dsqrt(0.5d0) in Fortran, 0x3fe6a09e667f3bcd)
#define imSqh   CMPLX(0.0, 0.707106781186547524401)
#define imIsqTwo CMPLX(0.0, 1.0 / 1.414213562373095048801688724209)

#define FR_min(a,b) ((a) < (b) ? (a) : (b))
#define FR_max(a,b) ((a) > (b) ? (a) : (b))

static inline AC_FP fsign(AC_FP a, AC_FP b) {
    return (b >= 0.0) ? fabs(a) : -fabs(a);
}

// Minkowski dot products (+,-,-,-)
static inline AC_CX dot4cx(const AC_CX a[6], const AC_CX b[6]) {
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
}
static inline AC_CX dot4cx_fp(const AC_CX a[6], const AC_D_FP b[4]) {
    return a[0]*(AC_FP)b[0] - a[1]*(AC_FP)b[1] - a[2]*(AC_FP)b[2] - a[3]*(AC_FP)b[3];
}
static inline AC_FP dot4fp(const AC_FP a[4], const AC_FP b[4]) {
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
}

// ─────────────────────────────────────────────────────────────────────────────
// External wavefunctions
// ─────────────────────────────────────────────────────────────────────────────

void ext_gluon_real(const AC_D_FP p_D[4], int ihel, AC_FP wf[4])
{
    AC_CX wf1[6], wf0[6];
    ext_gluon_cmplx(p_D,  1, wf1);
    ext_gluon_cmplx(p_D, -1, wf0);

    if (ihel == 1) {
        for (int i = 0; i < 4; ++i)
            wf[i] = creal(cImag * (wf1[i] + wf0[i])) * sqh;
    } else if (ihel == -1) {
        for (int i = 0; i < 4; ++i)
            wf[i] = -creal(wf1[i] - wf0[i]) * sqh;
    }
}

void ext_gluon_cmplx(const AC_D_FP p_D[4], int ihel, AC_CX wf[6])
{
    AC_FP p[4] = { (AC_FP)p_D[0], (AC_FP)p_D[1], (AC_FP)p_D[2], (AC_FP)p_D[3] };
    if (p[0] == 0.0) {
        fprintf(stderr, "Cannot generate external gluon with zero energy\n");
        fprintf(stderr, "%g %g %g %g\n", p[0], p[1], p[2], p[3]);
        exit(1);
    } else if (p[0] > 0.0) {
        AC_FP hel  = (AC_FP)ihel;
        AC_FP pp   = p[0];
        AC_FP pt   = sqrt(p[1]*p[1] + p[2]*p[2]);
        wf[0] = cZero;
        wf[3] = CMPLX(hel * pt / pp * sqh, 0.0);
        if (pt != rZero) {
            AC_FP pzpt = p[3] / (pp * pt) * sqh * hel;
            wf[1] = CMPLX(-p[1] * pzpt, -p[2] / pt * sqh);
            wf[2] = CMPLX(-p[2] * pzpt,  p[1] / pt * sqh);
        } else {
            wf[1] = CMPLX(-hel * sqh, 0.0);
            wf[2] = CMPLX(rZero, fsign(sqh, p[3]));
        }
    } else {
        AC_FP hel  = (AC_FP)(-ihel);
        AC_FP pp   = -p[0];
        AC_FP pt   = sqrt(p[1]*p[1] + p[2]*p[2]);
        wf[0] = cZero;
        wf[3] = CMPLX(hel * pt / pp * sqh, 0.0);
        if (pt != rZero) {
            AC_FP pzpt = -p[3] / (pp * pt) * sqh * hel;
            wf[1] = CMPLX( p[1] * pzpt,  p[2] / pt * sqh);
            wf[2] = CMPLX( p[2] * pzpt, -p[1] / pt * sqh);
        } else {
            wf[1] = CMPLX(-hel * sqh, 0.0);
            wf[2] = CMPLX(rZero, -fsign(sqh, p[3]));
        }
    }
}

void ext_vector_mass(const AC_D_FP p_D[4], int nhel, int nsv,
                     AC_CX wf[6], AC_FP vmass)
{
    AC_FP p[4] = { (AC_FP)p_D[0], (AC_FP)p_D[1], (AC_FP)p_D[2], (AC_FP)p_D[3] };
    AC_FP hel    = (AC_FP)nhel;
    int   nsvahl = nsv * abs(nhel);
    AC_FP pt2    = p[1]*p[1] + p[2]*p[2];
    AC_FP pp     = FR_min(p[0], sqrt(pt2 + p[3]*p[3]));
    AC_FP pt     = FR_min(pp,   sqrt(pt2));

    if (vmass != rZero) {
        AC_FP hel0 = rOne - fabs(hel);
        if (pp == rZero) {
            wf[0] = CMPLX(rZero, 0.0);
            wf[1] = CMPLX(-hel * sqh, 0.0);
            wf[2] = CMPLX(rZero, nsvahl * sqh);
            wf[3] = CMPLX(hel0, 0.0);
        } else {
            AC_FP emp = p[0] / (vmass * pp);
            wf[0] = CMPLX(hel0 * pp / vmass, 0.0);
            wf[3] = CMPLX(hel0 * p[3] * emp + hel * pt / pp * sqh, 0.0);
            if (pt != rZero) {
                AC_FP pzpt = p[3] / (pp * pt) * sqh * hel;
                wf[1] = CMPLX(hel0*p[1]*emp - p[1]*pzpt, -nsvahl*p[2]/pt*sqh);
                wf[2] = CMPLX(hel0*p[2]*emp - p[2]*pzpt,  nsvahl*p[1]/pt*sqh);
            } else {
                wf[1] = CMPLX(-hel * sqh, 0.0);
                wf[2] = CMPLX(rZero, nsvahl * fsign(sqh, p[3]));
            }
        }
    } else {
        pp = p[0];
        pt = sqrt(p[1]*p[1] + p[2]*p[2]);
        wf[0] = CMPLX(rZero, 0.0);
        wf[3] = CMPLX(hel * pt / pp * sqh, 0.0);
        if (pt != rZero) {
            AC_FP pzpt = p[3] / (pp * pt) * sqh * hel;
            wf[1] = CMPLX(-p[1]*pzpt, -nsv*p[2]/pt*sqh);
            wf[2] = CMPLX(-p[2]*pzpt,  nsv*p[1]/pt*sqh);
        } else {
            wf[1] = CMPLX(-hel * sqh, 0.0);
            wf[2] = CMPLX(rZero, nsv * fsign(sqh, p[3]));
        }
    }
}

void ext_quark(const AC_D_FP p_D[4], int nhel, AC_CX wf[6], AC_FP fmass)
{
    AC_FP p[4] = { (AC_FP)p_D[0], (AC_FP)p_D[1], (AC_FP)p_D[2], (AC_FP)p_D[3] };
    AC_CX chi[2];

    if (p[0] > 0.0) {
        if (fabs(fmass) < tiny) {
            AC_FP sqp0p3;
            if (p[1] == 0.0 && p[2] == 0.0 && p[3] < 0.0)
                sqp0p3 = 0.0;
            else if (p[3] < 0.0)
                sqp0p3 = sqrt(FR_max((p[1]*p[1] + p[2]*p[2]) / (p[0] - p[3]), rZero));
            else
                sqp0p3 = sqrt(FR_max(p[0] + p[3], rZero));
            chi[0] = CMPLX(sqp0p3, 0.0);
            if (sqp0p3 == rZero)
                chi[1] = (double complex)(-nhel) * sqrt(rTwo * p[0]);
            else
                chi[1] = CMPLX((double)(nhel) * p[1], -p[2]) / sqp0p3;
            if (nhel == 1) {
                wf[0] = chi[0]; wf[1] = chi[1];
                wf[2] = cZero;  wf[3] = cZero;
            } else {
                wf[0] = cZero;  wf[1] = cZero;
                wf[2] = chi[1]; wf[3] = chi[0];
            }
        } else {
            int   nsf  = +1;
            int   nh   = nsf * nhel;
            AC_FP pp   = fabs(sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]));
            AC_FP sf[2], omega[2], sfomeg[2];
            sf[0] = (AC_FP)(1 + nsf + (1 - nsf)*nh) * 0.5;
            sf[1] = (AC_FP)(1 + nsf - (1 - nsf)*nh) * 0.5;
            omega[0] = sqrt(p[0] + pp);
            omega[1] = fmass / omega[0];
            int ip = (3 + nh) / 2 - 1;
            int im = (3 - nh) / 2 - 1;
            sfomeg[0] = sf[0] * omega[ip];
            sfomeg[1] = sf[1] * omega[im];
            AC_FP pp3;
            if (p[3] < 0.0)
                pp3 = FR_max((p[1]*p[1] + p[2]*p[2]) / (pp - p[3]), rZero);
            else
                pp3 = FR_max(pp + p[3], rZero);
            chi[0] = CMPLX(sqrt(pp3 * 0.5 / pp), 0.0);
            if (pp3 == rZero)
                chi[1] = (double complex)(-nh);
            else
                chi[1] = CMPLX((double)nh * p[1], -p[2]) / sqrt(rTwo * pp * pp3);
            wf[0] = sfomeg[1] * chi[im];
            wf[1] = sfomeg[1] * chi[ip];
            wf[2] = sfomeg[0] * chi[im];
            wf[3] = sfomeg[0] * chi[ip];
        }
    } else {
        if (fabs(fmass) < tiny) {
            AC_FP sqp0p3;
            if (p[1] == 0.0 && p[2] == 0.0 && p[3] > 0.0)
                sqp0p3 = 0.0;
            else if (p[3] > 0.0)
                sqp0p3 = -sqrt(FR_max((p[1]*p[1] + p[2]*p[2]) / (p[3] - p[0]), rZero));
            else
                sqp0p3 = -sqrt(FR_max(-(p[0] + p[3]), rZero));
            chi[0] = CMPLX(sqp0p3, 0.0);
            if (sqp0p3 == rZero)
                chi[1] = (double complex)(-nhel) * sqrt(rTwo * fabs(p[0]));
            else
                chi[1] = CMPLX((double)(-nhel) * (-p[1]), -(-p[2])) / sqp0p3;
            if (-nhel == 1) {
                wf[0] = chi[0]; wf[1] = chi[1];
                wf[2] = cZero;  wf[3] = cZero;
            } else {
                wf[0] = cZero;  wf[1] = cZero;
                wf[2] = chi[1]; wf[3] = chi[0];
            }
        } else {
            int   nsf  = -1;
            int   nh   = nsf * nhel;
            AC_FP pp   = fabs(sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]));
            AC_FP sf[2], omega[2], sfomeg[2];
            sf[0] = (AC_FP)(1 + nsf + (1 - nsf)*nh) * 0.5;
            sf[1] = (AC_FP)(1 + nsf - (1 - nsf)*nh) * 0.5;
            omega[0] = sqrt(fabs(p[0]) + pp);
            omega[1] = fmass / omega[0];
            int ip = (3 + nh) / 2 - 1;
            int im = (3 - nh) / 2 - 1;
            sfomeg[0] = sf[0] * omega[ip];
            sfomeg[1] = sf[1] * omega[im];
            AC_FP pp3;
            if (p[3] > 0.0)
                pp3 = FR_max((p[1]*p[1] + p[2]*p[2]) / (pp + p[3]), rZero);
            else
                pp3 = FR_max(pp - p[3], rZero);
            chi[0] = CMPLX(sqrt(pp3 * 0.5 / pp), 0.0);
            if (pp3 == rZero)
                chi[1] = (double complex)(-nh);
            else
                chi[1] = CMPLX((double)nh * (-p[1]), -(-p[2])) / sqrt(rTwo * pp * pp3);
            wf[0] = sfomeg[1] * chi[im];
            wf[1] = sfomeg[1] * chi[ip];
            wf[2] = sfomeg[0] * chi[im];
            wf[3] = sfomeg[0] * chi[ip];
        }
    }
}

void ext_antiquark(const AC_D_FP p_D[4], int nhel, AC_CX wf[6], AC_FP fmass)
{
    AC_FP p[4] = { (AC_FP)p_D[0], (AC_FP)p_D[1], (AC_FP)p_D[2], (AC_FP)p_D[3] };
    AC_CX chi[2];

    if (p[0] > 0.0) {
        if (fabs(fmass) < tiny) {
            AC_FP sqp0p3;
            if (p[1] == 0.0 && p[2] == 0.0 && p[3] < 0.0)
                sqp0p3 = 0.0;
            else if (p[3] < 0.0)
                sqp0p3 = -sqrt(FR_max((p[1]*p[1] + p[2]*p[2]) / (p[0] - p[3]), rZero));
            else
                sqp0p3 = -sqrt(FR_max(p[0] + p[3], rZero));
            chi[0] = CMPLX(sqp0p3, 0.0);
            if (sqp0p3 == rZero)
                chi[1] = (double complex)(-nhel) * sqrt(rTwo * p[0]);
            else
                chi[1] = CMPLX((double)(-nhel) * p[1], p[2]) / sqp0p3;
            if (-nhel == 1) {
                wf[0] = cZero;  wf[1] = cZero;
                wf[2] = chi[0]; wf[3] = chi[1];
            } else {
                wf[0] = chi[1]; wf[1] = chi[0];
                wf[2] = cZero;  wf[3] = cZero;
            }
        } else {
            int   nsf  = -1;
            int   nh   = nsf * nhel;
            AC_FP pp   = fabs(sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]));
            AC_FP sf[2], omega[2], sfomeg[2];
            sf[0] = (AC_FP)(1 + nsf + (1 - nsf)*nh) * 0.5;
            sf[1] = (AC_FP)(1 + nsf - (1 - nsf)*nh) * 0.5;
            omega[0] = sqrt(p[0] + pp);
            omega[1] = fmass / omega[0];
            int ip = (3 + nh) / 2 - 1;
            int im = (3 - nh) / 2 - 1;
            sfomeg[0] = sf[0] * omega[ip];
            sfomeg[1] = sf[1] * omega[im];
            AC_FP pp3;
            if (p[3] < 0.0)
                pp3 = FR_max((p[1]*p[1] + p[2]*p[2]) / (pp - p[3]), rZero);
            else
                pp3 = FR_max(pp + p[3], rZero);
            chi[0] = CMPLX(sqrt(pp3 * 0.5 / pp), 0.0);
            if (pp3 == rZero)
                chi[1] = (double complex)(-nh);
            else
                chi[1] = CMPLX((double)nh * p[1], p[2]) / sqrt(rTwo * pp * pp3);
            wf[0] = sfomeg[0] * chi[im];
            wf[1] = sfomeg[0] * chi[ip];
            wf[2] = sfomeg[1] * chi[im];
            wf[3] = sfomeg[1] * chi[ip];
        }
    } else {
        if (fabs(fmass) < tiny) {
            AC_FP sqp0p3;
            if (p[1] == 0.0 && p[2] == 0.0 && p[3] > 0.0)
                sqp0p3 = 0.0;
            else if (p[3] > 0.0)
                sqp0p3 = sqrt(FR_max((p[1]*p[1] + p[2]*p[2]) / (p[3] - p[0]), rZero));
            else
                sqp0p3 = sqrt(FR_max(-(p[0] + p[3]), rZero));
            chi[0] = CMPLX(sqp0p3, 0.0);
            if (sqp0p3 == rZero)
                chi[1] = (double complex)(-nhel) * sqrt(rTwo * fabs(p[0]));
            else
                chi[1] = CMPLX((double)nhel * (-p[1]), (-p[2])) / sqp0p3;
            if (nhel == 1) {
                wf[0] = cZero;  wf[1] = cZero;
                wf[2] = chi[0]; wf[3] = chi[1];
            } else {
                wf[0] = chi[1]; wf[1] = chi[0];
                wf[2] = cZero;  wf[3] = cZero;
            }
        } else {
            int   nsf  = +1;
            int   nh   = nsf * nhel;
            AC_FP pp   = fabs(sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]));
            AC_FP sf[2], omega[2], sfomeg[2];
            sf[0] = (AC_FP)(1 + nsf + (1 - nsf)*nh) * 0.5;
            sf[1] = (AC_FP)(1 + nsf - (1 - nsf)*nh) * 0.5;
            omega[0] = sqrt(fabs(p[0]) + pp);
            omega[1] = fmass / omega[0];
            int ip = (3 + nh) / 2 - 1;
            int im = (3 - nh) / 2 - 1;
            sfomeg[0] = sf[0] * omega[ip];
            sfomeg[1] = sf[1] * omega[im];
            AC_FP pp3;
            if (p[3] > 0.0)
                pp3 = FR_max((p[1]*p[1] + p[2]*p[2]) / (pp + p[3]), rZero);
            else
                pp3 = FR_max(pp - p[3], rZero);
            chi[0] = CMPLX(sqrt(pp3 * 0.5 / pp), 0.0);
            if (pp3 == rZero)
                chi[1] = (double complex)(-nh);
            else
                chi[1] = CMPLX((double)nh * (-p[1]), (-p[2])) / sqrt(rTwo * pp * pp3);
            wf[0] = sfomeg[0] * chi[im];
            wf[1] = sfomeg[0] * chi[ip];
            wf[2] = sfomeg[1] * chi[im];
            wf[3] = sfomeg[1] * chi[ip];
        }
    }
}

void ext_scalar(const AC_D_FP p_D[4], AC_CX wf[1])
{
    (void)p_D;
    wf[0] = CMPLX(1.0, 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Vertices
// ─────────────────────────────────────────────────────────────────────────────

void ThreeGluon(const AC_CX wf1[6], const AC_D_FP pwf1[4],
                const AC_CX wf2[6], const AC_D_FP pwf2[4],
                AC_CX wf[6])
{
    AC_CX TMP1 = dot4cx(wf1, wf2);
    AC_CX TMP2 = dot4cx_fp(wf1, pwf2);
    AC_CX TMP3 = dot4cx_fp(wf2, pwf1);
    for (int i = 0; i < 4; ++i) {
        wf[i] = CMPLX(0.0, isqTwo) * (TMP1 * CMPLX(pwf1[i] - pwf2[i], 0.0) +
                        2.0 * (TMP2 * wf2[i] - TMP3 * wf1[i]));
    }
}

void ThreeGluon_real(const AC_FP wf1[4], const AC_D_FP pwf1[4],
                     const AC_FP wf2[4], const AC_D_FP pwf2[4],
                     AC_FP wf[4])
{
    AC_FP TMP1 = dot4fp(wf1, wf2);
    AC_FP TMP2 = wf1[0]*pwf2[0] - wf1[1]*pwf2[1] - wf1[2]*pwf2[2] - wf1[3]*pwf2[3];
    AC_FP TMP3 = wf2[0]*pwf1[0] - wf2[1]*pwf1[1] - wf2[2]*pwf1[2] - wf2[3]*pwf1[3];
    for (int i = 0; i < 4; ++i)
        wf[i] = isqTwo * (TMP1 * (pwf1[i] - pwf2[i]) +
                           2.0 * (TMP2 * wf2[i] - TMP3 * wf1[i]));
}

void ThreeGluon_coupl(const AC_CX wf1[6], const AC_D_FP pwf1[4],
                      const AC_CX wf2[6], const AC_D_FP pwf2[4],
                      AC_CX wf[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = dot4cx(wf1, wf2);
    AC_CX TMP2 = cZero, TMP3 = cZero;
    for (int i = 0; i < 4; ++i) {
        AC_FP sign = (i == 0) ? 1.0 : -1.0;
        TMP2 += sign * wf1[i] * (2.0*pwf2[i] + pwf1[i]);
        TMP3 += sign * wf2[i] * (-2.0*pwf1[i] - pwf2[i]);
    }
    AC_CX TMP4 = imSqh * coupl[0];
    for (int i = 0; i < 4; ++i)
        wf[i] = TMP4 * (TMP1 * CMPLX(pwf1[i] - pwf2[i], 0.0) +
                         TMP2 * wf2[i] + TMP3 * wf1[i]);
}

void FourGluon(const AC_CX wf1[6], const AC_CX wf2[6],
               const AC_CX wf3[6], AC_CX wf[6])
{
    AC_CX TMP1 = dot4cx(wf1, wf2);
    AC_CX TMP2 = dot4cx(wf1, wf3);
    AC_CX TMP3 = dot4cx(wf2, wf3);
    for (int i = 0; i < 4; ++i)
        wf[i] = imHalf * (2.0*wf2[i]*TMP2 - wf1[i]*TMP3 - wf3[i]*TMP1);
}

void TwoGluontoTensor(const AC_CX wfg1[6], const AC_CX wfg2[6], AC_CX wfT[6])
{
    wfT[0] = wfg1[0]*wfg2[1] - wfg1[1]*wfg2[0];
    wfT[1] = wfg1[0]*wfg2[2] - wfg1[2]*wfg2[0];
    wfT[2] = wfg1[0]*wfg2[3] - wfg1[3]*wfg2[0];
    wfT[3] = wfg1[1]*wfg2[2] - wfg1[2]*wfg2[1];
    wfT[4] = wfg1[1]*wfg2[3] - wfg1[3]*wfg2[1];
    wfT[5] = wfg1[2]*wfg2[3] - wfg1[3]*wfg2[2];
}

void TwoGluontoTensor_real(const AC_FP wfg1[4], const AC_FP wfg2[4], AC_FP wfT[6])
{
    wfT[0] = wfg1[0]*wfg2[1] - wfg1[1]*wfg2[0];
    wfT[1] = wfg1[0]*wfg2[2] - wfg1[2]*wfg2[0];
    wfT[2] = wfg1[0]*wfg2[3] - wfg1[3]*wfg2[0];
    wfT[3] = wfg1[1]*wfg2[2] - wfg1[2]*wfg2[1];
    wfT[4] = wfg1[1]*wfg2[3] - wfg1[3]*wfg2[1];
    wfT[5] = wfg1[2]*wfg2[3] - wfg1[3]*wfg2[2];
}

void TwoGluontoTensor_coupl(const AC_CX wfg1[6], const AC_CX wfg2[6],
                            AC_CX wfT[6], const AC_FP coupl[2])
{
    wfT[0] = (wfg1[0]*wfg2[1] - wfg1[1]*wfg2[0]) * coupl[0];
    wfT[1] = (wfg1[0]*wfg2[2] - wfg1[2]*wfg2[0]) * coupl[0];
    wfT[2] = (wfg1[0]*wfg2[3] - wfg1[3]*wfg2[0]) * coupl[0];
    wfT[3] = (wfg1[1]*wfg2[2] - wfg1[2]*wfg2[1]) * coupl[0];
    wfT[4] = (wfg1[1]*wfg2[3] - wfg1[3]*wfg2[1]) * coupl[0];
    wfT[5] = (wfg1[2]*wfg2[3] - wfg1[3]*wfg2[2]) * coupl[0];
}

void TensorGluontoGluon(const AC_CX wfT1[6], const AC_CX wfg2[6], AC_CX wfg[6])
{
    wfg[0] = (wfT1[0]*wfg2[1] + wfT1[1]*wfg2[2] + wfT1[2]*wfg2[3]) * imHalf;
    wfg[1] = (wfT1[0]*wfg2[0] + wfT1[3]*wfg2[2] + wfT1[4]*wfg2[3]) * imHalf;
    wfg[2] = (wfT1[1]*wfg2[0] - wfT1[3]*wfg2[1] + wfT1[5]*wfg2[3]) * imHalf;
    wfg[3] = (wfT1[2]*wfg2[0] - wfT1[4]*wfg2[1] - wfT1[5]*wfg2[2]) * imHalf;
}

void TensorGluontoGluon_real(const AC_FP wfT1[6], const AC_FP wfg2[4], AC_FP wfg[4])
{
    wfg[0] = (wfT1[0]*wfg2[1] + wfT1[1]*wfg2[2] + wfT1[2]*wfg2[3]) * rHalf;
    wfg[1] = (wfT1[0]*wfg2[0] + wfT1[3]*wfg2[2] + wfT1[4]*wfg2[3]) * rHalf;
    wfg[2] = (wfT1[1]*wfg2[0] - wfT1[3]*wfg2[1] + wfT1[5]*wfg2[3]) * rHalf;
    wfg[3] = (wfT1[2]*wfg2[0] - wfT1[4]*wfg2[1] - wfT1[5]*wfg2[2]) * rHalf;
}

void TensorGluontoGluon_coupl(const AC_CX wfT1[6], const AC_CX wfg2[6],
                              AC_CX wfg[6], const AC_FP coupl[2])
{
    wfg[0] = (wfT1[0]*wfg2[1] + wfT1[1]*wfg2[2] + wfT1[2]*wfg2[3]) * imHalf * coupl[0];
    wfg[1] = (wfT1[0]*wfg2[0] + wfT1[3]*wfg2[2] + wfT1[4]*wfg2[3]) * imHalf * coupl[0];
    wfg[2] = (wfT1[1]*wfg2[0] - wfT1[3]*wfg2[1] + wfT1[5]*wfg2[3]) * imHalf * coupl[0];
    wfg[3] = (wfT1[2]*wfg2[0] - wfT1[4]*wfg2[1] - wfT1[5]*wfg2[2]) * imHalf * coupl[0];
}

void GluonTensortoGluon(const AC_CX wfg1[6], const AC_CX wfT2[6], AC_CX wfg[6])
{
    wfg[0] = (-wfg1[1]*wfT2[0] - wfg1[2]*wfT2[1] - wfg1[3]*wfT2[2]) * imHalf;
    wfg[1] = (-wfg1[0]*wfT2[0] - wfg1[2]*wfT2[3] - wfg1[3]*wfT2[4]) * imHalf;
    wfg[2] = (-wfg1[0]*wfT2[1] + wfg1[1]*wfT2[3] - wfg1[3]*wfT2[5]) * imHalf;
    wfg[3] = (-wfg1[0]*wfT2[2] + wfg1[1]*wfT2[4] + wfg1[2]*wfT2[5]) * imHalf;
}

void GluonTensortoGluon_real(const AC_FP wfg1[4], const AC_FP wfT2[6], AC_FP wfg[4])
{
    wfg[0] = (-wfg1[1]*wfT2[0] - wfg1[2]*wfT2[1] - wfg1[3]*wfT2[2]) * rHalf;
    wfg[1] = (-wfg1[0]*wfT2[0] - wfg1[2]*wfT2[3] - wfg1[3]*wfT2[4]) * rHalf;
    wfg[2] = (-wfg1[0]*wfT2[1] + wfg1[1]*wfT2[3] - wfg1[3]*wfT2[5]) * rHalf;
    wfg[3] = (-wfg1[0]*wfT2[2] + wfg1[1]*wfT2[4] + wfg1[2]*wfT2[5]) * rHalf;
}

void GluonTensortoGluon_coupl(const AC_CX wfg1[6], const AC_CX wfT2[6],
                              AC_CX wfg[6], const AC_FP coupl[2])
{
    wfg[0] = (-wfg1[1]*wfT2[0] - wfg1[2]*wfT2[1] - wfg1[3]*wfT2[2]) * imHalf * coupl[0];
    wfg[1] = (-wfg1[0]*wfT2[0] - wfg1[2]*wfT2[3] - wfg1[3]*wfT2[4]) * imHalf * coupl[0];
    wfg[2] = (-wfg1[0]*wfT2[1] + wfg1[1]*wfT2[3] - wfg1[3]*wfT2[5]) * imHalf * coupl[0];
    wfg[3] = (-wfg1[0]*wfT2[2] + wfg1[1]*wfT2[4] + wfg1[2]*wfT2[5]) * imHalf * coupl[0];
}

// ─── Quark/gluon interactions ────────────────────────────────────────────────

void QuarkGluontoQuark(const AC_CX wfq1[6], const AC_CX wfg2[6], AC_CX wfq[6])
{
    AC_CX TMP1 = wfg2[0] + wfg2[3];
    AC_CX TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = wfg2[1] + cImag * wfg2[2];
    AC_CX TMP4 = wfg2[1] - cImag * wfg2[2];
    wfq[0] = imSqh * (TMP1*wfq1[2] + TMP3*wfq1[3]);
    wfq[1] = imSqh * (TMP2*wfq1[3] + TMP4*wfq1[2]);
    wfq[2] = imSqh * (TMP2*wfq1[0] - TMP3*wfq1[1]);
    wfq[3] = imSqh * (TMP1*wfq1[1] - TMP4*wfq1[0]);
}

void QuarkGluontoQuark_real(const AC_CX wfq1[6], const AC_FP wfg2[4], AC_CX wfq[6])
{
    AC_FP  TMP1 = wfg2[0] + wfg2[3];
    AC_FP  TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = CMPLX(wfg2[1],  wfg2[2]);
    AC_CX TMP4 = CMPLX(wfg2[1], -wfg2[2]);
    wfq[0] = imSqh * (TMP1*wfq1[2] + TMP3*wfq1[3]);
    wfq[1] = imSqh * (TMP2*wfq1[3] + TMP4*wfq1[2]);
    wfq[2] = imSqh * (TMP2*wfq1[0] - TMP3*wfq1[1]);
    wfq[3] = imSqh * (TMP1*wfq1[1] - TMP4*wfq1[0]);
}

void QuarkGluontoQuark_coupl(const AC_CX wfq1[6], const AC_CX wfg2[6],
                             AC_CX wfq[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfg2[0] + wfg2[3];
    AC_CX TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = wfg2[1] + cImag * wfg2[2];
    AC_CX TMP4 = wfg2[1] - cImag * wfg2[2];
    AC_CX TMP5 = imSqh * coupl[0];   // L
    wfq[0] = TMP5 * (TMP1*wfq1[2] + TMP3*wfq1[3]);
    wfq[1] = TMP5 * (TMP2*wfq1[3] + TMP4*wfq1[2]);
    TMP5   = imSqh * coupl[1];        // R
    wfq[2] = TMP5 * (TMP2*wfq1[0] - TMP3*wfq1[1]);
    wfq[3] = TMP5 * (TMP1*wfq1[1] - TMP4*wfq1[0]);
}

void GluonQuarktoQuark(const AC_CX wfg1[6], const AC_CX wfq2[6], AC_CX wfq[6])
{
    AC_CX TMP1 = wfg1[0] + wfg1[3];
    AC_CX TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = wfg1[1] + cImag * wfg1[2];
    AC_CX TMP4 = wfg1[1] - cImag * wfg1[2];
    wfq[0] = imSqh * (TMP1*wfq2[2] + TMP3*wfq2[3]);
    wfq[1] = imSqh * (TMP2*wfq2[3] + TMP4*wfq2[2]);
    wfq[2] = imSqh * (TMP2*wfq2[0] - TMP3*wfq2[1]);
    wfq[3] = imSqh * (TMP1*wfq2[1] - TMP4*wfq2[0]);
}

void GluonQuarktoQuark_real(const AC_FP wfg1[4], const AC_CX wfq2[6], AC_CX wfq[6])
{
    AC_FP  TMP1 = wfg1[0] + wfg1[3];
    AC_FP  TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = CMPLX(wfg1[1],  wfg1[2]);
    AC_CX TMP4 = CMPLX(wfg1[1], -wfg1[2]);
    wfq[0] = imSqh * (TMP1*wfq2[2] + TMP3*wfq2[3]);
    wfq[1] = imSqh * (TMP2*wfq2[3] + TMP4*wfq2[2]);
    wfq[2] = imSqh * (TMP2*wfq2[0] - TMP3*wfq2[1]);
    wfq[3] = imSqh * (TMP1*wfq2[1] - TMP4*wfq2[0]);
}

void GluonQuarktoQuark_coupl(const AC_CX wfg1[6], const AC_CX wfq2[6],
                             AC_CX wfq[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfg1[0] + wfg1[3];
    AC_CX TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = wfg1[1] + cImag * wfg1[2];
    AC_CX TMP4 = wfg1[1] - cImag * wfg1[2];
    // L
    wfq[0] = imSqh * (TMP1*wfq2[2] + TMP3*wfq2[3]) * coupl[0];
    wfq[1] = imSqh * (TMP2*wfq2[3] + TMP4*wfq2[2]) * coupl[0];
    // R
    wfq[2] = imSqh * (TMP2*wfq2[0] - TMP3*wfq2[1]) * coupl[1];
    wfq[3] = imSqh * (TMP1*wfq2[1] - TMP4*wfq2[0]) * coupl[1];
}

void AQuarkGluontoAQuark(const AC_CX wfq1[6], const AC_CX wfg2[6], AC_CX wfq[6])
{
    AC_CX TMP1 = wfg2[0] + wfg2[3];
    AC_CX TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = wfg2[1] + cImag * wfg2[2];
    AC_CX TMP4 = wfg2[1] - cImag * wfg2[2];
    wfq[0] = imSqh * (TMP2*wfq1[2] - TMP4*wfq1[3]);
    wfq[1] = imSqh * (TMP1*wfq1[3] - TMP3*wfq1[2]);
    wfq[2] = imSqh * (TMP1*wfq1[0] + TMP4*wfq1[1]);
    wfq[3] = imSqh * (TMP2*wfq1[1] + TMP3*wfq1[0]);
}

void AQuarkGluontoAQuark_coupl(const AC_CX wfq1[6], const AC_CX wfg2[6],
                               AC_CX wfq[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfg2[0] + wfg2[3];
    AC_CX TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = wfg2[1] + cImag * wfg2[2];
    AC_CX TMP4 = wfg2[1] - cImag * wfg2[2];
    AC_CX TMP5 = imSqh * coupl[1];
    wfq[0] = TMP5 * (TMP2*wfq1[2] - TMP4*wfq1[3]);
    wfq[1] = TMP5 * (TMP1*wfq1[3] - TMP3*wfq1[2]);
    TMP5   = imSqh * coupl[0];
    wfq[2] = TMP5 * (TMP1*wfq1[0] + TMP4*wfq1[1]);
    wfq[3] = TMP5 * (TMP2*wfq1[1] + TMP3*wfq1[0]);
}

void GluonAQuarktoAQuark(const AC_CX wfg1[6], const AC_CX wfq2[6], AC_CX wfq[6])
{
    AC_CX TMP1 = wfg1[0] + wfg1[3];
    AC_CX TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = wfg1[1] + cImag * wfg1[2];
    AC_CX TMP4 = wfg1[1] - cImag * wfg1[2];
    wfq[0] = imSqh * (TMP2*wfq2[2] - TMP4*wfq2[3]);
    wfq[1] = imSqh * (TMP1*wfq2[3] - TMP3*wfq2[2]);
    wfq[2] = imSqh * (TMP1*wfq2[0] + TMP4*wfq2[1]);
    wfq[3] = imSqh * (TMP2*wfq2[1] + TMP3*wfq2[0]);
}

void GluonAQuarktoAQuark_coupl(const AC_CX wfg1[6], const AC_CX wfq2[6],
                               AC_CX wfq[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfg1[0] + wfg1[3];
    AC_CX TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = wfg1[1] + cImag * wfg1[2];
    AC_CX TMP4 = wfg1[1] - cImag * wfg1[2];
    // L
    wfq[0] = imSqh * (TMP2*wfq2[2] - TMP4*wfq2[3]) * coupl[1];
    wfq[1] = imSqh * (TMP1*wfq2[3] - TMP3*wfq2[2]) * coupl[1];
    // R
    wfq[2] = imSqh * (TMP1*wfq2[0] + TMP4*wfq2[1]) * coupl[0];
    wfq[3] = imSqh * (TMP2*wfq2[1] + TMP3*wfq2[0]) * coupl[0];
}

void QuarkAQuarktoGluon(const AC_CX wfq1[6], const AC_CX wfq2[6],
                        AC_CX wfg[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfq1[2]*wfq2[0] + wfq1[1]*wfq2[3];
    AC_CX TMP2 = wfq1[3]*wfq2[1] + wfq1[0]*wfq2[2];
    AC_CX TMP3 = wfq1[1]*wfq2[2] - wfq1[3]*wfq2[0];
    AC_CX TMP4 = wfq1[0]*wfq2[3] - wfq1[2]*wfq2[1];
    wfg[0] = ( TMP1 + TMP2) * imSqh * coupl[0];
    wfg[1] = ( TMP4 + TMP3) * imSqh * coupl[0];
    wfg[2] = ( TMP4 - TMP3) * sqh   * coupl[0];
    wfg[3] = (-TMP1 + TMP2) * imSqh * coupl[0];
}

void AQuarkQuarktoGluon(const AC_CX wfq1[6], const AC_CX wfq2[6], AC_CX wfg[6])
{
    AC_CX TMP1 = wfq2[2]*wfq1[0] + wfq2[1]*wfq1[3];
    AC_CX TMP2 = wfq2[3]*wfq1[1] + wfq2[0]*wfq1[2];
    AC_CX TMP3 = wfq2[1]*wfq1[2] - wfq2[3]*wfq1[0];
    AC_CX TMP4 = wfq2[0]*wfq1[3] - wfq2[2]*wfq1[1];
    wfg[0] = ( TMP1 + TMP2) * imSqh;
    wfg[1] = ( TMP4 + TMP3) * imSqh;
    wfg[2] = ( TMP4 - TMP3) * sqh;
    wfg[3] = (-TMP1 + TMP2) * imSqh;
}

void LeptonALeptontoGluon(const AC_CX wfq1[6], const AC_CX wfq2[6],
                          AC_CX wfg[6], const AC_FP coupl[2])
{
    AC_CX wfg_temp[6];
    // L
    wfg_temp[0] = ( wfq1[2]*wfq2[0] + wfq1[3]*wfq2[1]) * imSqh * coupl[0];
    wfg_temp[1] = (-wfq1[3]*wfq2[0] - wfq1[2]*wfq2[1]) * imSqh * coupl[0];
    wfg_temp[2] = ( wfq1[3]*wfq2[0] - wfq1[2]*wfq2[1]) * sqh   * coupl[0];
    wfg_temp[3] = (-wfq1[2]*wfq2[0] + wfq1[3]*wfq2[1]) * imSqh * coupl[0];
    // R
    wfg[0] = ( wfq1[0]*wfq2[2] + wfq1[1]*wfq2[3]) * imSqh * coupl[1];
    wfg[1] = ( wfq1[0]*wfq2[3] + wfq1[1]*wfq2[2]) * imSqh * coupl[1];
    wfg[2] = ( wfq1[0]*wfq2[3] - wfq1[1]*wfq2[2]) * sqh   * coupl[1];
    wfg[3] = ( wfq1[0]*wfq2[2] - wfq1[1]*wfq2[3]) * imSqh * coupl[1];
    // add
    for (int i = 0; i < 4; ++i) wfg[i] += wfg_temp[i];
}

void ALeptonLeptontoGluon(const AC_CX wfq1[6], const AC_CX wfq2[6],
                          AC_CX wfg[6], const AC_FP coupl[2])
{
    AC_CX wfg_temp[6];
    // L
    wfg_temp[0] = ( wfq2[2]*wfq1[0] + wfq2[3]*wfq1[1]) * imSqh * coupl[0];
    wfg_temp[1] = (-wfq2[3]*wfq1[0] - wfq2[2]*wfq1[1]) * imSqh * coupl[0];
    wfg_temp[2] = ( wfq2[3]*wfq1[0] - wfq2[2]*wfq1[1]) * sqh   * coupl[0];
    wfg_temp[3] = (-wfq2[2]*wfq1[0] + wfq2[3]*wfq1[1]) * imSqh * coupl[0];
    // R
    wfg[0] = ( wfq2[0]*wfq1[2] + wfq2[1]*wfq1[3]) * imSqh * coupl[1];
    wfg[1] = ( wfq2[0]*wfq1[3] + wfq2[1]*wfq1[2]) * imSqh * coupl[1];
    wfg[2] = ( wfq2[0]*wfq1[3] - wfq2[1]*wfq1[2]) * sqh   * coupl[1];
    wfg[3] = ( wfq2[0]*wfq1[2] - wfq2[1]*wfq1[3]) * imSqh * coupl[1];
    // add
    for (int i = 0; i < 4; ++i) wfg[i] += wfg_temp[i];
}

void QuarkScalartoQuark(const AC_CX wfq1[6], const AC_CX wfs2[1],
                        AC_CX wfq[6], const AC_FP coupl[2])
{
    for (int i = 0; i < 4; ++i)
        wfq[i] = -imSqh * coupl[0] * wfs2[0] * wfq1[i];
}

void GluonGluontoScalar(const AC_CX wfg1[6], const AC_CX wfg2[6],
                        AC_CX wfs[1], const AC_FP coupl[2])
{
    AC_CX TMP = dot4cx(wfg1, wfg2);
    wfs[0] = imSqh * coupl[0] * TMP;
}

void ScalarGluontoGluon(const AC_CX wfs1[1], const AC_CX wfg2[6],
                        AC_CX wfg[6], const AC_FP coupl[2])
{
    for (int i = 0; i < 4; ++i)
        wfg[i] = imSqh * coupl[0] * wfs1[0] * wfg2[i];
}

void GluonScalartoGluon(const AC_CX wfg1[6], const AC_CX wfs2[1],
                        AC_CX wfg[6], const AC_FP coupl[2])
{
    for (int i = 0; i < 4; ++i)
        wfg[i] = imSqh * coupl[0] * wfs2[0] * wfg1[i];
}

void ScalarScalartoScalar(const AC_CX wfs1[1], const AC_CX wfs2[1],
                          AC_CX wfs[1], const AC_FP coupl[2])
{
    AC_CX TMP = CMPLX(1.0, 0.0);
    if (creal(coupl[1]) == -10.0 && cimag(coupl[1]) == 0.0) TMP = cImag;
    wfs[0] = imSqh * TMP * coupl[0] * wfs1[0] * wfs2[0];
}

// ─────────────────────────────────────────────────────────────────────────────
// Propagators
// ─────────────────────────────────────────────────────────────────────────────

void GluonPropagator(AC_CX wfg[6], const AC_D_FP p[4])
{
    // FP64: denominator
    AC_D_FP p2     = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX prop_d = -cImag / p2;
    // downcast 
    AC_CX prop = (AC_CX)prop_d;
    for (int i = 0; i < 4; ++i) wfg[i] *= prop;
}

void GluonPropagator_real(AC_FP wfg[4], const AC_D_FP p[4])
{
    // FP64: denominator
    AC_D_FP p2     = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_FP prop_d = 1.0 / p2;
    AC_FP prop = (AC_FP)prop_d;
    for (int i = 0; i < 4; ++i) wfg[i] *= prop;
}

void GluonPropagator_mass(AC_CX wfg[6], const AC_D_FP p[4], AC_D_FP vm, AC_D_FP vw)
{
    // FP64: denominator
    AC_D_FP p2     = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX denom  = CMPLX(p2 - vm*vm, 0.0) + cImag * vm * vw;
    AC_D_CX prop_d = -cImag / denom;
    AC_CX prop = (AC_CX)prop_d;

    // downcast 
    AC_FP pf[4] = { (AC_FP)p[0], (AC_FP)p[1], (AC_FP)p[2], (AC_FP)p[3] };
    AC_FP vmf   = (AC_FP)vm;
    AC_CX TMP   = (pf[0]*wfg[0] - pf[1]*wfg[1] - pf[2]*wfg[2] - pf[3]*wfg[3]) / (vmf * vmf);
    wfg[0] = (wfg[0] - pf[0]*TMP) * prop;
    wfg[1] = (wfg[1] - pf[1]*TMP) * prop;
    wfg[2] = (wfg[2] - pf[2]*TMP) * prop;
    wfg[3] = (wfg[3] - pf[3]*TMP) * prop;
}

void QuarkPropagator(AC_CX wfq[6], const AC_D_FP p[4], AC_D_FP fm, AC_D_FP fw)
{
    // FP64: denominator
    AC_D_FP p2      = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX denom   = CMPLX(p2 - fm*fm, 0.0) + cImag * fm * fw;
    AC_D_CX prefact_d = cImag / denom;
    AC_CX prefact = (AC_CX)prefact_d;

    // downcast 
    AC_FP pf[4] = { (AC_FP)p[0], (AC_FP)p[1], (AC_FP)p[2], (AC_FP)p[3] };
    AC_FP fmf   = (AC_FP)fm;
    AC_CX tmp[4] = {wfq[0], wfq[1], wfq[2], wfq[3]};
    AC_CX tp1 = CMPLX(pf[0] + pf[3], 0.0);
    AC_CX tp2 = CMPLX(pf[0] - pf[3], 0.0);
    AC_CX tp3 = CMPLX(pf[1],  pf[2]);
    AC_CX tp4 = CMPLX(pf[1], -pf[2]);
    wfq[0] = (tp1*tmp[2] + tp3*tmp[3] + fmf*tmp[0]) * prefact;
    wfq[1] = (tp2*tmp[3] + tp4*tmp[2] + fmf*tmp[1]) * prefact;
    wfq[2] = (tp2*tmp[0] - tp3*tmp[1] + fmf*tmp[2]) * prefact;
    wfq[3] = (tp1*tmp[1] - tp4*tmp[0] + fmf*tmp[3]) * prefact;
}

void AQuarkPropagator(AC_CX wfq[6], const AC_D_FP p[4], AC_D_FP fm, AC_D_FP fw)
{
    // FP64: denominator
    AC_D_FP p2      = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX denom   = CMPLX(p2 - fm*fm, 0.0) + cImag * fm * fw;
    AC_D_CX prefact_d = cImag / denom;
    AC_CX prefact = (AC_CX)prefact_d;

    // downcast 
    AC_FP pf[4] = { (AC_FP)p[0], (AC_FP)p[1], (AC_FP)p[2], (AC_FP)p[3] };
    AC_FP fmf   = (AC_FP)fm;
    AC_CX tmp[4] = {wfq[0], wfq[1], wfq[2], wfq[3]};
    AC_CX tp1 = CMPLX(-(pf[0] + pf[3]), 0.0);
    AC_CX tp2 = CMPLX(-(pf[0] - pf[3]), 0.0);
    AC_CX tp3 = CMPLX(-pf[1], -pf[2]);
    AC_CX tp4 = CMPLX(-pf[1],  pf[2]);
    wfq[0] = (tp2*tmp[2] - tp4*tmp[3] + fmf*tmp[0]) * prefact;
    wfq[1] = (tp1*tmp[3] - tp3*tmp[2] + fmf*tmp[1]) * prefact;
    wfq[2] = (tp1*tmp[0] + tp4*tmp[1] + fmf*tmp[2]) * prefact;
    wfq[3] = (tp2*tmp[1] + tp3*tmp[0] + fmf*tmp[3]) * prefact;
}

void ScalarPropagator(AC_CX wfs[1], const AC_D_FP p[4], AC_D_FP sm, AC_D_FP sw)
{
    // FP64: denominator 
    AC_D_FP p2     = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX denom  = CMPLX(p2 - sm*sm, 0.0) + cImag * sm * sw;
    AC_CX prop = (AC_CX)(cImag / denom);
    wfs[0] *= prop;
}
