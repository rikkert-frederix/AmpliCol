// FeynmanRules_device.cu
// Device-side implementations of all Feynman rule functions declared in
// FeynmanRules_device.cuh.
//
// Compile with relocatable device code so other .cu translation units can
// call these symbols:
//   nvcc -rdc=true -c FeynmanRules_device.cu -o FeynmanRules_device.o
//   nvcc -dlink FeynmanRules_device.o <other .o's> -o ...
//
// CUDA C++ only (nvcc). AC_CX = std::complex<double> throughout.
// Replaces FeynmanRules.h + FeynmanRules.c for device compilation.
//
// Differences from the C host version:
//   • CMPLX(r,i)   → AC_CX(r,i)       (C++ std::complex constructor)
//   • creal(z)     → (z).real()
//   • fprintf/exit → __trap()          (non-recoverable bad input)

#include "FeynmanRules_device.cuh"

// ─── constants ─────────────────────────────────────────────────────────────

static __device__ const AC_FP d_rZero  = 0.0;
static __device__ const AC_FP d_rHalf  = 0.5;
static __device__ const AC_FP d_rOne   = 1.0;
static __device__ const AC_FP d_rTwo   = 2.0;
static __device__ const AC_FP d_sqh    = 0.707106781186547524401;
static __device__ const AC_FP d_isqTwo = 1.0 / 1.414213562373095048801688724209;
static __device__ const AC_FP d_tiny   = 1e-8;

// ─── helpers ───────────────────────────────────────────────────────────────

static __device__ __forceinline__ AC_FP d_fsign(AC_FP a, AC_FP b) {
    return (b >= 0.0) ? fabs(a) : -fabs(a);
}

static __device__ __forceinline__ AC_CX d_dot4cx(const AC_CX a[6], const AC_CX b[6]) {
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
}
static __device__ __forceinline__ AC_CX d_dot4cx_fp(const AC_CX a[6], const AC_D_FP b[4]) {
    return a[0]*(AC_FP)b[0] - a[1]*(AC_FP)b[1] - a[2]*(AC_FP)b[2] - a[3]*(AC_FP)b[3];
}
static __device__ __forceinline__ AC_FP d_dot4fp(const AC_FP a[4], const AC_FP b[4]) {
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
}

// ─────────────────────────────────────────────────────────────────────────────
// External wavefunctions
// ─────────────────────────────────────────────────────────────────────────────

__device__
void ext_gluon_cmplx(const AC_D_FP p_D[4], int ihel, AC_CX wf[6])
{
    AC_FP p[4] = { (AC_FP)p_D[0], (AC_FP)p_D[1], (AC_FP)p_D[2], (AC_FP)p_D[3] };
    if (p[0] == 0.0) {
        __trap();   // zero-energy gluon: invalid input
    } else if (p[0] > 0.0) {
        AC_FP hel  = (AC_FP)ihel;
        AC_FP pp   = p[0];
        AC_FP pt   = sqrt(p[1]*p[1] + p[2]*p[2]);
        wf[0] = AC_CX(0.0, 0.0);
        wf[3] = AC_CX(hel * pt / pp * d_sqh, 0.0);
        if (pt != d_rZero) {
            AC_FP pzpt = p[3] / (pp * pt) * d_sqh * hel;
            wf[1] = AC_CX(-p[1] * pzpt, -p[2] / pt * d_sqh);
            wf[2] = AC_CX(-p[2] * pzpt,  p[1] / pt * d_sqh);
        } else {
            wf[1] = AC_CX(-hel * d_sqh, 0.0);
            wf[2] = AC_CX(d_rZero, d_fsign(d_sqh, p[3]));
        }
    } else {
        AC_FP hel  = (AC_FP)(-ihel);
        AC_FP pp   = -p[0];
        AC_FP pt   = sqrt(p[1]*p[1] + p[2]*p[2]);
        wf[0] = AC_CX(0.0, 0.0);
        wf[3] = AC_CX(hel * pt / pp * d_sqh, 0.0);
        if (pt != d_rZero) {
            AC_FP pzpt = -p[3] / (pp * pt) * d_sqh * hel;
            wf[1] = AC_CX( p[1] * pzpt,  p[2] / pt * d_sqh);
            wf[2] = AC_CX( p[2] * pzpt, -p[1] / pt * d_sqh);
        } else {
            wf[1] = AC_CX(-hel * d_sqh, 0.0);
            wf[2] = AC_CX(d_rZero, -d_fsign(d_sqh, p[3]));
        }
    }
}

__device__
void ext_gluon_real(const AC_D_FP p_D[4], int ihel, AC_FP wf[4])
{
    AC_CX wf1[6], wf0[6];
    ext_gluon_cmplx(p_D,  1, wf1);
    ext_gluon_cmplx(p_D, -1, wf0);
    if (ihel == 1) {
        for (int i = 0; i < 4; ++i)
            wf[i] = (AC_CX(0.0,1.0) * (wf1[i] + wf0[i])).real() * d_sqh;
    } else if (ihel == -1) {
        for (int i = 0; i < 4; ++i)
            wf[i] = -(wf1[i] - wf0[i]).real() * d_sqh;
    }
}

__device__
void ext_vector_mass(const AC_D_FP p_D[4], int nhel, int nsv,
                     AC_CX wf[6], AC_FP vmass)
{
    AC_FP p[4] = { (AC_FP)p_D[0], (AC_FP)p_D[1], (AC_FP)p_D[2], (AC_FP)p_D[3] };
    AC_FP hel    = (AC_FP)nhel;
    int   nsvahl = nsv * abs(nhel);
    AC_FP pt2    = p[1]*p[1] + p[2]*p[2];
    AC_FP pp     = fmin(p[0], sqrt(pt2 + p[3]*p[3]));
    AC_FP pt     = fmin(pp,   sqrt(pt2));

    if (vmass != d_rZero) {
        AC_FP hel0 = d_rOne - fabs(hel);
        if (pp == d_rZero) {
            wf[0] = AC_CX(d_rZero, 0.0);
            wf[1] = AC_CX(-hel * d_sqh, 0.0);
            wf[2] = AC_CX(d_rZero, nsvahl * d_sqh);
            wf[3] = AC_CX(hel0, 0.0);
        } else {
            AC_FP emp = p[0] / (vmass * pp);
            wf[0] = AC_CX(hel0 * pp / vmass, 0.0);
            wf[3] = AC_CX(hel0 * p[3] * emp + hel * pt / pp * d_sqh, 0.0);
            if (pt != d_rZero) {
                AC_FP pzpt = p[3] / (pp * pt) * d_sqh * hel;
                wf[1] = AC_CX(hel0*p[1]*emp - p[1]*pzpt, -nsvahl*p[2]/pt*d_sqh);
                wf[2] = AC_CX(hel0*p[2]*emp - p[2]*pzpt,  nsvahl*p[1]/pt*d_sqh);
            } else {
                wf[1] = AC_CX(-hel * d_sqh, 0.0);
                wf[2] = AC_CX(d_rZero, nsvahl * d_fsign(d_sqh, p[3]));
            }
        }
    } else {
        pp = p[0];
        pt = sqrt(p[1]*p[1] + p[2]*p[2]);
        wf[0] = AC_CX(d_rZero, 0.0);
        wf[3] = AC_CX(hel * pt / pp * d_sqh, 0.0);
        if (pt != d_rZero) {
            AC_FP pzpt = p[3] / (pp * pt) * d_sqh * hel;
            wf[1] = AC_CX(-p[1]*pzpt, -nsv*p[2]/pt*d_sqh);
            wf[2] = AC_CX(-p[2]*pzpt,  nsv*p[1]/pt*d_sqh);
        } else {
            wf[1] = AC_CX(-hel * d_sqh, 0.0);
            wf[2] = AC_CX(d_rZero, nsv * d_fsign(d_sqh, p[3]));
        }
    }
}

__device__
void ext_quark(const AC_D_FP p_D[4], int nhel, AC_CX wf[6], AC_FP fmass)
{
    AC_FP p[4] = { (AC_FP)p_D[0], (AC_FP)p_D[1], (AC_FP)p_D[2], (AC_FP)p_D[3] };
    AC_CX chi[2];

    if (p[0] > 0.0) {
        if (fabs(fmass) < d_tiny) {
            AC_FP sqp0p3;
            if (p[1] == 0.0 && p[2] == 0.0 && p[3] < 0.0)
                sqp0p3 = 0.0;
            else if (p[3] < 0.0)
                sqp0p3 = sqrt(fmax((p[1]*p[1] + p[2]*p[2]) / (p[0] - p[3]), d_rZero));
            else
                sqp0p3 = sqrt(fmax(p[0] + p[3], d_rZero));
            chi[0] = AC_CX(sqp0p3, 0.0);
            if (sqp0p3 == d_rZero)
                chi[1] = AC_CX(-nhel, 0.0) * sqrt(d_rTwo * p[0]);
            else
                chi[1] = AC_CX((double)(nhel) * p[1], -p[2]) / sqp0p3;
            if (nhel == 1) {
                wf[0] = chi[0]; wf[1] = chi[1];
                wf[2] = AC_CX(0,0); wf[3] = AC_CX(0,0);
            } else {
                wf[0] = AC_CX(0,0); wf[1] = AC_CX(0,0);
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
                pp3 = fmax((p[1]*p[1] + p[2]*p[2]) / (pp - p[3]), d_rZero);
            else
                pp3 = fmax(pp + p[3], d_rZero);
            chi[0] = AC_CX(sqrt(pp3 * 0.5 / pp), 0.0);
            if (pp3 == d_rZero)
                chi[1] = AC_CX(-nh, 0.0);
            else
                chi[1] = AC_CX((double)nh * p[1], -p[2]) / sqrt(d_rTwo * pp * pp3);
            wf[0] = sfomeg[1] * chi[im];
            wf[1] = sfomeg[1] * chi[ip];
            wf[2] = sfomeg[0] * chi[im];
            wf[3] = sfomeg[0] * chi[ip];
        }
    } else {
        if (fabs(fmass) < d_tiny) {
            AC_FP sqp0p3;
            if (p[1] == 0.0 && p[2] == 0.0 && p[3] > 0.0)
                sqp0p3 = 0.0;
            else if (p[3] > 0.0)
                sqp0p3 = -sqrt(fmax((p[1]*p[1] + p[2]*p[2]) / (p[3] - p[0]), d_rZero));
            else
                sqp0p3 = -sqrt(fmax(-(p[0] + p[3]), d_rZero));
            chi[0] = AC_CX(sqp0p3, 0.0);
            if (sqp0p3 == d_rZero)
                chi[1] = AC_CX(-nhel, 0.0) * sqrt(d_rTwo * fabs(p[0]));
            else
                chi[1] = AC_CX((double)(-nhel) * (-p[1]), -(-p[2])) / sqp0p3;
            if (-nhel == 1) {
                wf[0] = chi[0]; wf[1] = chi[1];
                wf[2] = AC_CX(0,0); wf[3] = AC_CX(0,0);
            } else {
                wf[0] = AC_CX(0,0); wf[1] = AC_CX(0,0);
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
                pp3 = fmax((p[1]*p[1] + p[2]*p[2]) / (pp + p[3]), d_rZero);
            else
                pp3 = fmax(pp - p[3], d_rZero);
            chi[0] = AC_CX(sqrt(pp3 * 0.5 / pp), 0.0);
            if (pp3 == d_rZero)
                chi[1] = AC_CX(-nh, 0.0);
            else
                chi[1] = AC_CX((double)nh * (-p[1]), -(-p[2])) / sqrt(d_rTwo * pp * pp3);
            wf[0] = sfomeg[1] * chi[im];
            wf[1] = sfomeg[1] * chi[ip];
            wf[2] = sfomeg[0] * chi[im];
            wf[3] = sfomeg[0] * chi[ip];
        }
    }
}

__device__
void ext_antiquark(const AC_D_FP p_D[4], int nhel, AC_CX wf[6], AC_FP fmass)
{
    AC_FP p[4] = { (AC_FP)p_D[0], (AC_FP)p_D[1], (AC_FP)p_D[2], (AC_FP)p_D[3] };
    AC_CX chi[2];

    if (p[0] > 0.0) {
        if (fabs(fmass) < d_tiny) {
            AC_FP sqp0p3;
            if (p[1] == 0.0 && p[2] == 0.0 && p[3] < 0.0)
                sqp0p3 = 0.0;
            else if (p[3] < 0.0)
                sqp0p3 = -sqrt(fmax((p[1]*p[1] + p[2]*p[2]) / (p[0] - p[3]), d_rZero));
            else
                sqp0p3 = -sqrt(fmax(p[0] + p[3], d_rZero));
            chi[0] = AC_CX(sqp0p3, 0.0);
            if (sqp0p3 == d_rZero)
                chi[1] = AC_CX(-nhel, 0.0) * sqrt(d_rTwo * p[0]);
            else
                chi[1] = AC_CX((double)(-nhel) * p[1], p[2]) / sqp0p3;
            if (-nhel == 1) {
                wf[0] = AC_CX(0,0); wf[1] = AC_CX(0,0);
                wf[2] = chi[0]; wf[3] = chi[1];
            } else {
                wf[0] = chi[1]; wf[1] = chi[0];
                wf[2] = AC_CX(0,0); wf[3] = AC_CX(0,0);
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
                pp3 = fmax((p[1]*p[1] + p[2]*p[2]) / (pp - p[3]), d_rZero);
            else
                pp3 = fmax(pp + p[3], d_rZero);
            chi[0] = AC_CX(sqrt(pp3 * 0.5 / pp), 0.0);
            if (pp3 == d_rZero)
                chi[1] = AC_CX(-nh, 0.0);
            else
                chi[1] = AC_CX((double)nh * p[1], p[2]) / sqrt(d_rTwo * pp * pp3);
            wf[0] = sfomeg[0] * chi[im];
            wf[1] = sfomeg[0] * chi[ip];
            wf[2] = sfomeg[1] * chi[im];
            wf[3] = sfomeg[1] * chi[ip];
        }
    } else {
        if (fabs(fmass) < d_tiny) {
            AC_FP sqp0p3;
            if (p[1] == 0.0 && p[2] == 0.0 && p[3] > 0.0)
                sqp0p3 = 0.0;
            else if (p[3] > 0.0)
                sqp0p3 = sqrt(fmax((p[1]*p[1] + p[2]*p[2]) / (p[3] - p[0]), d_rZero));
            else
                sqp0p3 = sqrt(fmax(-(p[0] + p[3]), d_rZero));
            chi[0] = AC_CX(sqp0p3, 0.0);
            if (sqp0p3 == d_rZero)
                chi[1] = AC_CX(-nhel, 0.0) * sqrt(d_rTwo * fabs(p[0]));
            else
                chi[1] = AC_CX((double)nhel * (-p[1]), (-p[2])) / sqp0p3;
            if (nhel == 1) {
                wf[0] = AC_CX(0,0); wf[1] = AC_CX(0,0);
                wf[2] = chi[0]; wf[3] = chi[1];
            } else {
                wf[0] = chi[1]; wf[1] = chi[0];
                wf[2] = AC_CX(0,0); wf[3] = AC_CX(0,0);
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
                pp3 = fmax((p[1]*p[1] + p[2]*p[2]) / (pp + p[3]), d_rZero);
            else
                pp3 = fmax(pp - p[3], d_rZero);
            chi[0] = AC_CX(sqrt(pp3 * 0.5 / pp), 0.0);
            if (pp3 == d_rZero)
                chi[1] = AC_CX(-nh, 0.0);
            else
                chi[1] = AC_CX((double)nh * (-p[1]), (-p[2])) / sqrt(d_rTwo * pp * pp3);
            wf[0] = sfomeg[0] * chi[im];
            wf[1] = sfomeg[0] * chi[ip];
            wf[2] = sfomeg[1] * chi[im];
            wf[3] = sfomeg[1] * chi[ip];
        }
    }
}

__device__
void ext_scalar(const AC_D_FP p_D[4], AC_CX wf[1])
{
    (void)p_D;
    wf[0] = AC_CX(1.0, 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Vertices
// ─────────────────────────────────────────────────────────────────────────────

__device__
void ThreeGluon(const AC_CX wf1[6], const AC_D_FP pwf1[4],
                const AC_CX wf2[6], const AC_D_FP pwf2[4],
                AC_CX wf[6])
{
    AC_CX TMP1 = d_dot4cx(wf1, wf2);
    AC_CX TMP2 = d_dot4cx_fp(wf1, pwf2);
    AC_CX TMP3 = d_dot4cx_fp(wf2, pwf1);
    for (int i = 0; i < 4; ++i) {
        wf[i] = AC_CX(0.0, d_isqTwo) * (TMP1 * AC_CX(pwf1[i] - pwf2[i], 0.0) +
                        2.0 * (TMP2 * wf2[i] - TMP3 * wf1[i]));
    }
}

__device__
void ThreeGluon_real(const AC_FP wf1[4], const AC_D_FP pwf1[4],
                     const AC_FP wf2[4], const AC_D_FP pwf2[4],
                     AC_FP wf[4])
{
    AC_FP TMP1 = d_dot4fp(wf1, wf2);
    AC_FP TMP2 = wf1[0]*pwf2[0] - wf1[1]*pwf2[1] - wf1[2]*pwf2[2] - wf1[3]*pwf2[3];
    AC_FP TMP3 = wf2[0]*pwf1[0] - wf2[1]*pwf1[1] - wf2[2]*pwf1[2] - wf2[3]*pwf1[3];
    for (int i = 0; i < 4; ++i)
        wf[i] = d_isqTwo * (TMP1 * (pwf1[i] - pwf2[i]) +
                             2.0 * (TMP2 * wf2[i] - TMP3 * wf1[i]));
}

__device__
void ThreeGluon_coupl(const AC_CX wf1[6], const AC_D_FP pwf1[4],
                      const AC_CX wf2[6], const AC_D_FP pwf2[4],
                      AC_CX wf[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = d_dot4cx(wf1, wf2);
    AC_CX TMP2 = AC_CX(0,0), TMP3 = AC_CX(0,0);
    for (int i = 0; i < 4; ++i) {
        AC_FP sign = (i == 0) ? 1.0 : -1.0;
        TMP2 += sign * wf1[i] * (2.0*pwf2[i] + pwf1[i]);
        TMP3 += sign * wf2[i] * (-2.0*pwf1[i] - pwf2[i]);
    }
    AC_CX TMP4 = AC_CX(0.0, d_sqh) * coupl[0];
    for (int i = 0; i < 4; ++i)
        wf[i] = TMP4 * (TMP1 * AC_CX(pwf1[i] - pwf2[i], 0.0) +
                         TMP2 * wf2[i] + TMP3 * wf1[i]);
}

__device__
void FourGluon(const AC_CX wf1[6], const AC_CX wf2[6],
               const AC_CX wf3[6], AC_CX wf[6])
{
    AC_CX TMP1 = d_dot4cx(wf1, wf2);
    AC_CX TMP2 = d_dot4cx(wf1, wf3);
    AC_CX TMP3 = d_dot4cx(wf2, wf3);
    for (int i = 0; i < 4; ++i)
        wf[i] = AC_CX(0.0, 0.5) * (2.0*wf2[i]*TMP2 - wf1[i]*TMP3 - wf3[i]*TMP1);
}

__device__
void TwoGluontoTensor(const AC_CX wfg1[6], const AC_CX wfg2[6], AC_CX wfT[6])
{
    wfT[0] = wfg1[0]*wfg2[1] - wfg1[1]*wfg2[0];
    wfT[1] = wfg1[0]*wfg2[2] - wfg1[2]*wfg2[0];
    wfT[2] = wfg1[0]*wfg2[3] - wfg1[3]*wfg2[0];
    wfT[3] = wfg1[1]*wfg2[2] - wfg1[2]*wfg2[1];
    wfT[4] = wfg1[1]*wfg2[3] - wfg1[3]*wfg2[1];
    wfT[5] = wfg1[2]*wfg2[3] - wfg1[3]*wfg2[2];
}

__device__
void TwoGluontoTensor_real(const AC_FP wfg1[4], const AC_FP wfg2[4], AC_FP wfT[6])
{
    wfT[0] = wfg1[0]*wfg2[1] - wfg1[1]*wfg2[0];
    wfT[1] = wfg1[0]*wfg2[2] - wfg1[2]*wfg2[0];
    wfT[2] = wfg1[0]*wfg2[3] - wfg1[3]*wfg2[0];
    wfT[3] = wfg1[1]*wfg2[2] - wfg1[2]*wfg2[1];
    wfT[4] = wfg1[1]*wfg2[3] - wfg1[3]*wfg2[1];
    wfT[5] = wfg1[2]*wfg2[3] - wfg1[3]*wfg2[2];
}

__device__
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

__device__
void TensorGluontoGluon(const AC_CX wfT1[6], const AC_CX wfg2[6], AC_CX wfg[6])
{
    wfg[0] = (wfT1[0]*wfg2[1] + wfT1[1]*wfg2[2] + wfT1[2]*wfg2[3]) * AC_CX(0.0, 0.5);
    wfg[1] = (wfT1[0]*wfg2[0] + wfT1[3]*wfg2[2] + wfT1[4]*wfg2[3]) * AC_CX(0.0, 0.5);
    wfg[2] = (wfT1[1]*wfg2[0] - wfT1[3]*wfg2[1] + wfT1[5]*wfg2[3]) * AC_CX(0.0, 0.5);
    wfg[3] = (wfT1[2]*wfg2[0] - wfT1[4]*wfg2[1] - wfT1[5]*wfg2[2]) * AC_CX(0.0, 0.5);
}

__device__
void TensorGluontoGluon_real(const AC_FP wfT1[6], const AC_FP wfg2[4], AC_FP wfg[4])
{
    wfg[0] = (wfT1[0]*wfg2[1] + wfT1[1]*wfg2[2] + wfT1[2]*wfg2[3]) * d_rHalf;
    wfg[1] = (wfT1[0]*wfg2[0] + wfT1[3]*wfg2[2] + wfT1[4]*wfg2[3]) * d_rHalf;
    wfg[2] = (wfT1[1]*wfg2[0] - wfT1[3]*wfg2[1] + wfT1[5]*wfg2[3]) * d_rHalf;
    wfg[3] = (wfT1[2]*wfg2[0] - wfT1[4]*wfg2[1] - wfT1[5]*wfg2[2]) * d_rHalf;
}

__device__
void TensorGluontoGluon_coupl(const AC_CX wfT1[6], const AC_CX wfg2[6],
                              AC_CX wfg[6], const AC_FP coupl[2])
{
    wfg[0] = (wfT1[0]*wfg2[1] + wfT1[1]*wfg2[2] + wfT1[2]*wfg2[3]) * AC_CX(0.0,0.5) * coupl[0];
    wfg[1] = (wfT1[0]*wfg2[0] + wfT1[3]*wfg2[2] + wfT1[4]*wfg2[3]) * AC_CX(0.0,0.5) * coupl[0];
    wfg[2] = (wfT1[1]*wfg2[0] - wfT1[3]*wfg2[1] + wfT1[5]*wfg2[3]) * AC_CX(0.0,0.5) * coupl[0];
    wfg[3] = (wfT1[2]*wfg2[0] - wfT1[4]*wfg2[1] - wfT1[5]*wfg2[2]) * AC_CX(0.0,0.5) * coupl[0];
}

__device__
void GluonTensortoGluon(const AC_CX wfg1[6], const AC_CX wfT2[6], AC_CX wfg[6])
{
    wfg[0] = (-wfg1[1]*wfT2[0] - wfg1[2]*wfT2[1] - wfg1[3]*wfT2[2]) * AC_CX(0.0, 0.5);
    wfg[1] = (-wfg1[0]*wfT2[0] - wfg1[2]*wfT2[3] - wfg1[3]*wfT2[4]) * AC_CX(0.0, 0.5);
    wfg[2] = (-wfg1[0]*wfT2[1] + wfg1[1]*wfT2[3] - wfg1[3]*wfT2[5]) * AC_CX(0.0, 0.5);
    wfg[3] = (-wfg1[0]*wfT2[2] + wfg1[1]*wfT2[4] + wfg1[2]*wfT2[5]) * AC_CX(0.0, 0.5);
}

__device__
void GluonTensortoGluon_real(const AC_FP wfg1[4], const AC_FP wfT2[6], AC_FP wfg[4])
{
    wfg[0] = (-wfg1[1]*wfT2[0] - wfg1[2]*wfT2[1] - wfg1[3]*wfT2[2]) * d_rHalf;
    wfg[1] = (-wfg1[0]*wfT2[0] - wfg1[2]*wfT2[3] - wfg1[3]*wfT2[4]) * d_rHalf;
    wfg[2] = (-wfg1[0]*wfT2[1] + wfg1[1]*wfT2[3] - wfg1[3]*wfT2[5]) * d_rHalf;
    wfg[3] = (-wfg1[0]*wfT2[2] + wfg1[1]*wfT2[4] + wfg1[2]*wfT2[5]) * d_rHalf;
}

__device__
void GluonTensortoGluon_coupl(const AC_CX wfg1[6], const AC_CX wfT2[6],
                              AC_CX wfg[6], const AC_FP coupl[2])
{
    wfg[0] = (-wfg1[1]*wfT2[0] - wfg1[2]*wfT2[1] - wfg1[3]*wfT2[2]) * AC_CX(0.0,0.5) * coupl[0];
    wfg[1] = (-wfg1[0]*wfT2[0] - wfg1[2]*wfT2[3] - wfg1[3]*wfT2[4]) * AC_CX(0.0,0.5) * coupl[0];
    wfg[2] = (-wfg1[0]*wfT2[1] + wfg1[1]*wfT2[3] - wfg1[3]*wfT2[5]) * AC_CX(0.0,0.5) * coupl[0];
    wfg[3] = (-wfg1[0]*wfT2[2] + wfg1[1]*wfT2[4] + wfg1[2]*wfT2[5]) * AC_CX(0.0,0.5) * coupl[0];
}

// ─── Quark/gluon interactions ────────────────────────────────────────────────

__device__
void QuarkGluontoQuark(const AC_CX wfq1[6], const AC_CX wfg2[6], AC_CX wfq[6])
{
    AC_CX TMP1 = wfg2[0] + wfg2[3];
    AC_CX TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = wfg2[1] + AC_CX(0,1) * wfg2[2];
    AC_CX TMP4 = wfg2[1] - AC_CX(0,1) * wfg2[2];
    wfq[0] = AC_CX(0.0, d_sqh) * (TMP1*wfq1[2] + TMP3*wfq1[3]);
    wfq[1] = AC_CX(0.0, d_sqh) * (TMP2*wfq1[3] + TMP4*wfq1[2]);
    wfq[2] = AC_CX(0.0, d_sqh) * (TMP2*wfq1[0] - TMP3*wfq1[1]);
    wfq[3] = AC_CX(0.0, d_sqh) * (TMP1*wfq1[1] - TMP4*wfq1[0]);
}

__device__
void QuarkGluontoQuark_real(const AC_CX wfq1[6], const AC_FP wfg2[4], AC_CX wfq[6])
{
    AC_FP  TMP1 = wfg2[0] + wfg2[3];
    AC_FP  TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = AC_CX(wfg2[1],  wfg2[2]);
    AC_CX TMP4 = AC_CX(wfg2[1], -wfg2[2]);
    wfq[0] = AC_CX(0.0, d_sqh) * (TMP1*wfq1[2] + TMP3*wfq1[3]);
    wfq[1] = AC_CX(0.0, d_sqh) * (TMP2*wfq1[3] + TMP4*wfq1[2]);
    wfq[2] = AC_CX(0.0, d_sqh) * (TMP2*wfq1[0] - TMP3*wfq1[1]);
    wfq[3] = AC_CX(0.0, d_sqh) * (TMP1*wfq1[1] - TMP4*wfq1[0]);
}

__device__
void QuarkGluontoQuark_coupl(const AC_CX wfq1[6], const AC_CX wfg2[6],
                             AC_CX wfq[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfg2[0] + wfg2[3];
    AC_CX TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = wfg2[1] + AC_CX(0,1) * wfg2[2];
    AC_CX TMP4 = wfg2[1] - AC_CX(0,1) * wfg2[2];
    AC_CX TMP5 = AC_CX(0.0, d_sqh) * coupl[0];
    wfq[0] = TMP5 * (TMP1*wfq1[2] + TMP3*wfq1[3]);
    wfq[1] = TMP5 * (TMP2*wfq1[3] + TMP4*wfq1[2]);
    TMP5   = AC_CX(0.0, d_sqh) * coupl[1];
    wfq[2] = TMP5 * (TMP2*wfq1[0] - TMP3*wfq1[1]);
    wfq[3] = TMP5 * (TMP1*wfq1[1] - TMP4*wfq1[0]);
}

__device__
void GluonQuarktoQuark(const AC_CX wfg1[6], const AC_CX wfq2[6], AC_CX wfq[6])
{
    AC_CX TMP1 = wfg1[0] + wfg1[3];
    AC_CX TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = wfg1[1] + AC_CX(0,1) * wfg1[2];
    AC_CX TMP4 = wfg1[1] - AC_CX(0,1) * wfg1[2];
    wfq[0] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[2] + TMP3*wfq2[3]);
    wfq[1] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[3] + TMP4*wfq2[2]);
    wfq[2] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[0] - TMP3*wfq2[1]);
    wfq[3] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[1] - TMP4*wfq2[0]);
}

__device__
void GluonQuarktoQuark_real(const AC_FP wfg1[4], const AC_CX wfq2[6], AC_CX wfq[6])
{
    AC_FP  TMP1 = wfg1[0] + wfg1[3];
    AC_FP  TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = AC_CX(wfg1[1],  wfg1[2]);
    AC_CX TMP4 = AC_CX(wfg1[1], -wfg1[2]);
    wfq[0] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[2] + TMP3*wfq2[3]);
    wfq[1] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[3] + TMP4*wfq2[2]);
    wfq[2] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[0] - TMP3*wfq2[1]);
    wfq[3] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[1] - TMP4*wfq2[0]);
}

__device__
void GluonQuarktoQuark_coupl(const AC_CX wfg1[6], const AC_CX wfq2[6],
                             AC_CX wfq[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfg1[0] + wfg1[3];
    AC_CX TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = wfg1[1] + AC_CX(0,1) * wfg1[2];
    AC_CX TMP4 = wfg1[1] - AC_CX(0,1) * wfg1[2];
    wfq[0] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[2] + TMP3*wfq2[3]) * coupl[0];
    wfq[1] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[3] + TMP4*wfq2[2]) * coupl[0];
    wfq[2] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[0] - TMP3*wfq2[1]) * coupl[1];
    wfq[3] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[1] - TMP4*wfq2[0]) * coupl[1];
}

__device__
void AQuarkGluontoAQuark(const AC_CX wfq1[6], const AC_CX wfg2[6], AC_CX wfq[6])
{
    AC_CX TMP1 = wfg2[0] + wfg2[3];
    AC_CX TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = wfg2[1] + AC_CX(0,1) * wfg2[2];
    AC_CX TMP4 = wfg2[1] - AC_CX(0,1) * wfg2[2];
    wfq[0] = AC_CX(0.0, d_sqh) * (TMP2*wfq1[2] - TMP4*wfq1[3]);
    wfq[1] = AC_CX(0.0, d_sqh) * (TMP1*wfq1[3] - TMP3*wfq1[2]);
    wfq[2] = AC_CX(0.0, d_sqh) * (TMP1*wfq1[0] + TMP4*wfq1[1]);
    wfq[3] = AC_CX(0.0, d_sqh) * (TMP2*wfq1[1] + TMP3*wfq1[0]);
}

__device__
void AQuarkGluontoAQuark_coupl(const AC_CX wfq1[6], const AC_CX wfg2[6],
                               AC_CX wfq[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfg2[0] + wfg2[3];
    AC_CX TMP2 = wfg2[0] - wfg2[3];
    AC_CX TMP3 = wfg2[1] + AC_CX(0,1) * wfg2[2];
    AC_CX TMP4 = wfg2[1] - AC_CX(0,1) * wfg2[2];
    AC_CX TMP5 = AC_CX(0.0, d_sqh) * coupl[1];
    wfq[0] = TMP5 * (TMP2*wfq1[2] - TMP4*wfq1[3]);
    wfq[1] = TMP5 * (TMP1*wfq1[3] - TMP3*wfq1[2]);
    TMP5   = AC_CX(0.0, d_sqh) * coupl[0];
    wfq[2] = TMP5 * (TMP1*wfq1[0] + TMP4*wfq1[1]);
    wfq[3] = TMP5 * (TMP2*wfq1[1] + TMP3*wfq1[0]);
}

__device__
void GluonAQuarktoAQuark(const AC_CX wfg1[6], const AC_CX wfq2[6], AC_CX wfq[6])
{
    AC_CX TMP1 = wfg1[0] + wfg1[3];
    AC_CX TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = wfg1[1] + AC_CX(0,1) * wfg1[2];
    AC_CX TMP4 = wfg1[1] - AC_CX(0,1) * wfg1[2];
    wfq[0] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[2] - TMP4*wfq2[3]);
    wfq[1] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[3] - TMP3*wfq2[2]);
    wfq[2] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[0] + TMP4*wfq2[1]);
    wfq[3] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[1] + TMP3*wfq2[0]);
}

__device__
void GluonAQuarktoAQuark_coupl(const AC_CX wfg1[6], const AC_CX wfq2[6],
                               AC_CX wfq[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfg1[0] + wfg1[3];
    AC_CX TMP2 = wfg1[0] - wfg1[3];
    AC_CX TMP3 = wfg1[1] + AC_CX(0,1) * wfg1[2];
    AC_CX TMP4 = wfg1[1] - AC_CX(0,1) * wfg1[2];
    wfq[0] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[2] - TMP4*wfq2[3]) * coupl[1];
    wfq[1] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[3] - TMP3*wfq2[2]) * coupl[1];
    wfq[2] = AC_CX(0.0, d_sqh) * (TMP1*wfq2[0] + TMP4*wfq2[1]) * coupl[0];
    wfq[3] = AC_CX(0.0, d_sqh) * (TMP2*wfq2[1] + TMP3*wfq2[0]) * coupl[0];
}

__device__
void QuarkAQuarktoGluon(const AC_CX wfq1[6], const AC_CX wfq2[6],
                        AC_CX wfg[6], const AC_FP coupl[2])
{
    AC_CX TMP1 = wfq1[2]*wfq2[0] + wfq1[1]*wfq2[3];
    AC_CX TMP2 = wfq1[3]*wfq2[1] + wfq1[0]*wfq2[2];
    AC_CX TMP3 = wfq1[1]*wfq2[2] - wfq1[3]*wfq2[0];
    AC_CX TMP4 = wfq1[0]*wfq2[3] - wfq1[2]*wfq2[1];
    wfg[0] = ( TMP1 + TMP2) * AC_CX(0.0, d_sqh) * coupl[0];
    wfg[1] = ( TMP4 + TMP3) * AC_CX(0.0, d_sqh) * coupl[0];
    wfg[2] = ( TMP4 - TMP3) * d_sqh              * coupl[0];
    wfg[3] = (-TMP1 + TMP2) * AC_CX(0.0, d_sqh) * coupl[0];
}

__device__
void AQuarkQuarktoGluon(const AC_CX wfq1[6], const AC_CX wfq2[6], AC_CX wfg[6])
{
    AC_CX TMP1 = wfq2[2]*wfq1[0] + wfq2[1]*wfq1[3];
    AC_CX TMP2 = wfq2[3]*wfq1[1] + wfq2[0]*wfq1[2];
    AC_CX TMP3 = wfq2[1]*wfq1[2] - wfq2[3]*wfq1[0];
    AC_CX TMP4 = wfq2[0]*wfq1[3] - wfq2[2]*wfq1[1];
    wfg[0] = ( TMP1 + TMP2) * AC_CX(0.0, d_sqh);
    wfg[1] = ( TMP4 + TMP3) * AC_CX(0.0, d_sqh);
    wfg[2] = ( TMP4 - TMP3) * d_sqh;
    wfg[3] = (-TMP1 + TMP2) * AC_CX(0.0, d_sqh);
}

__device__
void LeptonALeptontoGluon(const AC_CX wfq1[6], const AC_CX wfq2[6],
                          AC_CX wfg[6], const AC_FP coupl[2])
{
    AC_CX wfg_temp[6];
    wfg_temp[0] = ( wfq1[2]*wfq2[0] + wfq1[3]*wfq2[1]) * AC_CX(0.0, d_sqh) * coupl[0];
    wfg_temp[1] = (-wfq1[3]*wfq2[0] - wfq1[2]*wfq2[1]) * AC_CX(0.0, d_sqh) * coupl[0];
    wfg_temp[2] = ( wfq1[3]*wfq2[0] - wfq1[2]*wfq2[1]) * d_sqh              * coupl[0];
    wfg_temp[3] = (-wfq1[2]*wfq2[0] + wfq1[3]*wfq2[1]) * AC_CX(0.0, d_sqh) * coupl[0];
    wfg[0] = ( wfq1[0]*wfq2[2] + wfq1[1]*wfq2[3]) * AC_CX(0.0, d_sqh) * coupl[1];
    wfg[1] = ( wfq1[0]*wfq2[3] + wfq1[1]*wfq2[2]) * AC_CX(0.0, d_sqh) * coupl[1];
    wfg[2] = ( wfq1[0]*wfq2[3] - wfq1[1]*wfq2[2]) * d_sqh              * coupl[1];
    wfg[3] = ( wfq1[0]*wfq2[2] - wfq1[1]*wfq2[3]) * AC_CX(0.0, d_sqh) * coupl[1];
    for (int i = 0; i < 4; ++i) wfg[i] += wfg_temp[i];
}

__device__
void ALeptonLeptontoGluon(const AC_CX wfq1[6], const AC_CX wfq2[6],
                          AC_CX wfg[6], const AC_FP coupl[2])
{
    AC_CX wfg_temp[6];
    wfg_temp[0] = ( wfq2[2]*wfq1[0] + wfq2[3]*wfq1[1]) * AC_CX(0.0, d_sqh) * coupl[0];
    wfg_temp[1] = (-wfq2[3]*wfq1[0] - wfq2[2]*wfq1[1]) * AC_CX(0.0, d_sqh) * coupl[0];
    wfg_temp[2] = ( wfq2[3]*wfq1[0] - wfq2[2]*wfq1[1]) * d_sqh              * coupl[0];
    wfg_temp[3] = (-wfq2[2]*wfq1[0] + wfq2[3]*wfq1[1]) * AC_CX(0.0, d_sqh) * coupl[0];
    wfg[0] = ( wfq2[0]*wfq1[2] + wfq2[1]*wfq1[3]) * AC_CX(0.0, d_sqh) * coupl[1];
    wfg[1] = ( wfq2[0]*wfq1[3] + wfq2[1]*wfq1[2]) * AC_CX(0.0, d_sqh) * coupl[1];
    wfg[2] = ( wfq2[0]*wfq1[3] - wfq2[1]*wfq1[2]) * d_sqh              * coupl[1];
    wfg[3] = ( wfq2[0]*wfq1[2] - wfq2[1]*wfq1[3]) * AC_CX(0.0, d_sqh) * coupl[1];
    for (int i = 0; i < 4; ++i) wfg[i] += wfg_temp[i];
}

__device__
void QuarkScalartoQuark(const AC_CX wfq1[6], const AC_CX wfs2[1],
                        AC_CX wfq[6], const AC_FP coupl[2])
{
    for (int i = 0; i < 4; ++i)
        wfq[i] = -AC_CX(0.0, d_sqh) * coupl[0] * wfs2[0] * wfq1[i];
}

__device__
void GluonGluontoScalar(const AC_CX wfg1[6], const AC_CX wfg2[6],
                        AC_CX wfs[1], const AC_FP coupl[2])
{
    AC_CX TMP = d_dot4cx(wfg1, wfg2);
    wfs[0] = AC_CX(0.0, d_sqh) * coupl[0] * TMP;
}

__device__
void ScalarGluontoGluon(const AC_CX wfs1[1], const AC_CX wfg2[6],
                        AC_CX wfg[6], const AC_FP coupl[2])
{
    for (int i = 0; i < 4; ++i)
        wfg[i] = AC_CX(0.0, d_sqh) * coupl[0] * wfs1[0] * wfg2[i];
}

__device__
void GluonScalartoGluon(const AC_CX wfg1[6], const AC_CX wfs2[1],
                        AC_CX wfg[6], const AC_FP coupl[2])
{
    for (int i = 0; i < 4; ++i)
        wfg[i] = AC_CX(0.0, d_sqh) * coupl[0] * wfs2[0] * wfg1[i];
}

__device__
void ScalarScalartoScalar(const AC_CX wfs1[1], const AC_CX wfs2[1],
                          AC_CX wfs[1], const AC_FP coupl[2])
{
    AC_CX TMP = AC_CX(1.0, 0.0);
    if (coupl[1] == -10.0) TMP = AC_CX(0.0, 1.0);
    wfs[0] = AC_CX(0.0, d_sqh) * TMP * coupl[0] * wfs1[0] * wfs2[0];
}

// ─────────────────────────────────────────────────────────────────────────────
// Propagators
// ─────────────────────────────────────────────────────────────────────────────


__device__
void GluonPropagator(AC_CX wfg[6], const AC_D_FP p[4])
{
    // FP64: denominator
    AC_D_FP p2     = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX prop_d = -AC_D_CX(0.0, 1.0) / p2;
    // downcast once, rest at working precision (trivial here: a single multiply)
    AC_CX prop = AC_CX(prop_d);
    for (int i = 0; i < 4; ++i) wfg[i] *= prop;
}

__device__
void GluonPropagator_real(AC_FP wfg[4], const AC_D_FP p[4])
{
    AC_D_FP p2     = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_FP prop_d = 1.0 / p2;
    AC_FP prop = (AC_FP)prop_d;
    for (int i = 0; i < 4; ++i) wfg[i] *= prop;
}

__device__
void GluonPropagator_mass(AC_CX wfg[6], const AC_D_FP p[4], AC_D_FP vm, AC_D_FP vw)
{
    // FP64: denominator
    AC_D_FP p2     = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX denom  = AC_D_CX(p2 - vm*vm, 0.0) + AC_D_CX(0.0, 1.0) * vm * vw;
    AC_D_CX prop_d = -AC_D_CX(0.0, 1.0) / denom;
    AC_CX prop = AC_CX(prop_d);

    // downcast
    AC_FP pf[4] = { (AC_FP)p[0], (AC_FP)p[1], (AC_FP)p[2], (AC_FP)p[3] };
    AC_FP vmf   = (AC_FP)vm;
    AC_CX TMP   = (pf[0]*wfg[0] - pf[1]*wfg[1] - pf[2]*wfg[2] - pf[3]*wfg[3]) / (vmf * vmf);
    wfg[0] = (wfg[0] - pf[0]*TMP) * prop;
    wfg[1] = (wfg[1] - pf[1]*TMP) * prop;
    wfg[2] = (wfg[2] - pf[2]*TMP) * prop;
    wfg[3] = (wfg[3] - pf[3]*TMP) * prop;
}

__device__
void QuarkPropagator(AC_CX wfq[6], const AC_D_FP p[4], AC_D_FP fm, AC_D_FP fw)
{
    // FP64: denominator
    AC_D_FP p2        = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX denom     = AC_D_CX(p2 - fm*fm, 0.0) + AC_D_CX(0.0, 1.0) * fm * fw;
    AC_D_CX prefact_d = AC_D_CX(0.0, 1.0) / denom;
    AC_CX prefact = AC_CX(prefact_d);

    // downcast
    AC_FP pf[4] = { (AC_FP)p[0], (AC_FP)p[1], (AC_FP)p[2], (AC_FP)p[3] };
    AC_FP fmf   = (AC_FP)fm;
    AC_CX tmp[4] = {wfq[0], wfq[1], wfq[2], wfq[3]};
    AC_CX tp1 = AC_CX(pf[0] + pf[3], 0.0);
    AC_CX tp2 = AC_CX(pf[0] - pf[3], 0.0);
    AC_CX tp3 = AC_CX(pf[1],  pf[2]);
    AC_CX tp4 = AC_CX(pf[1], -pf[2]);
    wfq[0] = (tp1*tmp[2] + tp3*tmp[3] + fmf*tmp[0]) * prefact;
    wfq[1] = (tp2*tmp[3] + tp4*tmp[2] + fmf*tmp[1]) * prefact;
    wfq[2] = (tp2*tmp[0] - tp3*tmp[1] + fmf*tmp[2]) * prefact;
    wfq[3] = (tp1*tmp[1] - tp4*tmp[0] + fmf*tmp[3]) * prefact;
}

__device__
void AQuarkPropagator(AC_CX wfq[6], const AC_D_FP p[4], AC_D_FP fm, AC_D_FP fw)
{
    // FP64: denominator
    AC_D_FP p2        = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX denom     = AC_D_CX(p2 - fm*fm, 0.0) + AC_D_CX(0.0, 1.0) * fm * fw;
    AC_D_CX prefact_d = AC_D_CX(0.0, 1.0) / denom;
    AC_CX prefact = AC_CX(prefact_d);

    // downcast
    AC_FP pf[4] = { (AC_FP)p[0], (AC_FP)p[1], (AC_FP)p[2], (AC_FP)p[3] };
    AC_FP fmf   = (AC_FP)fm;
    AC_CX tmp[4] = {wfq[0], wfq[1], wfq[2], wfq[3]};
    AC_CX tp1 = AC_CX(-(pf[0] + pf[3]), 0.0);
    AC_CX tp2 = AC_CX(-(pf[0] - pf[3]), 0.0);
    AC_CX tp3 = AC_CX(-pf[1], -pf[2]);
    AC_CX tp4 = AC_CX(-pf[1],  pf[2]);
    wfq[0] = (tp2*tmp[2] - tp4*tmp[3] + fmf*tmp[0]) * prefact;
    wfq[1] = (tp1*tmp[3] - tp3*tmp[2] + fmf*tmp[1]) * prefact;
    wfq[2] = (tp1*tmp[0] + tp4*tmp[1] + fmf*tmp[2]) * prefact;
    wfq[3] = (tp2*tmp[1] + tp3*tmp[0] + fmf*tmp[3]) * prefact;
}

__device__
void ScalarPropagator(AC_CX wfs[1], const AC_D_FP p[4], AC_D_FP sm, AC_D_FP sw)
{
    // FP64: denominator 
    AC_D_FP p2    = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    AC_D_CX denom = AC_D_CX(p2 - sm*sm, 0.0) + AC_D_CX(0.0, 1.0) * sm * sw;
    AC_CX prop = AC_CX(AC_D_CX(0.0, 1.0) / denom);
    wfs[0] *= prop;
}
