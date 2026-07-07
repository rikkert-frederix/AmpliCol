#include "umami.h"
#include "AmpliColTypes.h"

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// External amplitude library
// ---------------------------------------------------------------------------

extern "C" void evaluate_amp(
    int const* ichan,
    int const* iint,
    double const* p,
    std::complex<double>* amps
);

// ---------------------------------------------------------------------------
// Translation of the Fortran `umami_impl` module, unnamed namespace
// ---------------------------------------------------------------------------

namespace {

static constexpr double PI = 3.14159265358979323846;

struct ChannelGroup {
    int n        = 0;   // particle count
    int nvals    = 0;   // number of processes / output values
    int n_sing   = 0;   // number of singular (EW) vertices

    std::vector<int>    namps;       // number of helicity amplitudes per iint
    std::vector<int>    iproc_start; // process start indices (summed-process mode)
    std::vector<double> col_fac;     // colour factors
    std::vector<double> iden;        // identical-particle symmetry factors
    // hel_fac[hel][iint] stored as a flat vector, row-major: index = hel + max_hel*iint
    int                 max_hel  = 0;
    int                 max_iint = 0;
    std::vector<double> hel_fac; // size max_hel * (keep_separate ? max_iint : 1)
};

// Module-level state (mirrors the Fortran module variables)
static int                      nchans               = 0;
static std::vector<ChannelGroup> chans;
static bool                     keep_processes_separate = false;
static double                   alphaEW              = 0.0;
// Max amplitudes across all (ichan,iint) pairs — evaluate_amp always writes this many values.
static int                      max_amps_global      = 0;

// ------------------------------------------------------------------
// read_amplitude_lib_light  (was a Fortran subroutine)
// ------------------------------------------------------------------

static void read_amplitude_lib_light()
{
    const char* filename = "Library/amplitudes_light.bin";
    std::ifstream f(filename, std::ios::binary);
    if (!f)
        throw std::runtime_error(std::string("Cannot open ") + filename);

    auto read_scalar = [&](auto& val) {
        f.read(reinterpret_cast<char*>(&val), sizeof(val));
    };
    auto read_array = [&](auto* ptr, std::size_t count) {
        f.read(reinterpret_cast<char*>(ptr), sizeof(*ptr) * count);
    };

    int32_t tmp = 0;
    read_scalar(tmp);
    keep_processes_separate = (tmp != 0);
    read_scalar(nchans);
    read_scalar(alphaEW);

    chans.resize(nchans);
    for (int ichan = 0; ichan < nchans; ++ichan) {
        ChannelGroup& cg = chans[ichan];
        int max_iint = 0, max_hel = 0;

        read_scalar(cg.n);
        read_scalar(cg.nvals);
        read_scalar(cg.n_sing);
        read_scalar(max_iint);
        read_scalar(max_hel);

        // namps: size max_iint if nvals==1, else size 1  (mirrors Fortran allocation)
        int namps_size = (cg.nvals == 1) ? max_iint : 1;
        cg.namps.resize(namps_size);
        read_array(cg.namps.data(), namps_size);
        for (int n : cg.namps)
            max_amps_global = std::max(max_amps_global, n);

        if (!keep_processes_separate) {
            int idim = 0;
            read_scalar(idim);
            if (idim != cg.nvals + 1) {
                fprintf(stderr, "iproc_start has wrong size %d %d %d\n",
                        ichan + 1, idim, cg.nvals + 1);
                std::exit(1);
            }
            cg.iproc_start.resize(idim);
            read_array(cg.iproc_start.data(), idim);
        }

        {
            int idim = 0;
            read_scalar(idim);
            if (idim != max_iint) {
                fprintf(stderr, "col_fac has wrong size %d %d %d\n",
                        ichan + 1, idim, max_iint);
                std::exit(1);
            }
            cg.col_fac.resize(idim);
            read_array(cg.col_fac.data(), idim);
        }

        {
            int idim = 0;
            read_scalar(idim);
            if (idim != max_iint) {
                fprintf(stderr, "iden has wrong size %d %d %d\n",
                        ichan + 1, idim, max_iint);
                std::exit(1);
            }
            cg.iden.resize(idim);
            read_array(cg.iden.data(), idim);
        }

        {
            int idim = 0;
            read_scalar(idim);
            int hel_cols = keep_processes_separate ? max_iint : 1;
            cg.max_hel  = max_hel;
            cg.max_iint = max_iint;
            std::size_t total = static_cast<std::size_t>(max_hel) *
                                static_cast<std::size_t>(hel_cols);
            cg.hel_fac.resize(total);
            read_array(cg.hel_fac.data(), total);
        }
    }
}

// ------------------------------------------------------------------
// pick_one  (was a Fortran subroutine)
// Returns a 1-based index.
// ------------------------------------------------------------------

static int pick_one(int isize, const double* vals, double rnd)
{
    double sum = 0.0;
    for (int i = 0; i < isize; ++i) sum += vals[i];
    double target = rnd * sum;

    double accum   = vals[0];
    int    ipicked = 1;
    while (accum <= target) {
        // Bounds check BEFORE the next read. The original (and the Fortran it
        // was ported from) read vals[ipicked] one past the end when target
        // lands at/just below the full sum due to floating-point rounding. On
        // the stack that over-read was harmless; with a heap std::vector it
        // corrupts allocator metadata and crashes later at free().
        if (ipicked >= isize) {
            // target rounded to (or above) the running sum: clamp to last bin.
            ipicked = isize;
            break;
        }
        accum += vals[ipicked]; // 0-based read of the next element
        ++ipicked;
    }
    return ipicked; // 1-based, in [1, isize]
}

// ------------------------------------------------------------------
// unwgt_helicities  (was a Fortran subroutine)
// hels elements are 1-based on output (matching Fortran convention).
// ------------------------------------------------------------------

static void unwgt_helicities(
    int ichan_0,   // 0-based channel index
    int iint_0,    // 0-based flavor index
    const double* amp2_hel,
    const double* rnd,
    int*          hels)
{
    const ChannelGroup& cg = chans[ichan_0];
    if (keep_processes_separate) {
        hels[0] = pick_one(cg.namps[iint_0], amp2_hel, rnd[0]);
    } else {
        for (int iproc = 0; iproc < cg.nvals; ++iproc) {
            int start = cg.iproc_start[iproc];     // 0-based after translation
            int isize = cg.iproc_start[iproc + 1] - start;
            hels[iproc] = pick_one(isize, amp2_hel + start, rnd[iproc]);
        }
    }
}

// ------------------------------------------------------------------
// test_library  (was a Fortran subroutine)
// ------------------------------------------------------------------

static void test_library()
{
    for (int ichan = 0; ichan < nchans; ++ichan) {
        const ChannelGroup& cg = chans[ichan];
        int max_iint = keep_processes_separate ? cg.nvals : 1;

        for (int iint = 0; iint < max_iint; ++iint) {
            std::string path = "Library/amp" + std::to_string(ichan + 1) +
                               "_" + std::to_string(iint + 1) + "_libc.data";
            std::ifstream f(path, std::ios::binary);
            if (!f)
                throw std::runtime_error("Cannot open test file: " + path);

            int n     = cg.n;
            int namps = cg.namps[iint];

            std::vector<double> p(4 * static_cast<std::size_t>(n));
            std::vector<std::complex<double>> amps_save(static_cast<std::size_t>(namps));
            // evaluate_amp always writes max_amps_global values regardless of (ichan,iint)
            std::vector<std::complex<double>> amps(static_cast<std::size_t>(max_amps_global));

            f.read(reinterpret_cast<char*>(p.data()),
                   static_cast<std::streamsize>(sizeof(double) * p.size()));
            f.read(reinterpret_cast<char*>(amps_save.data()),
                   static_cast<std::streamsize>(sizeof(std::complex<double>) * namps));

            int ichan1 = ichan + 1;
            int iint1  = iint  + 1;
            evaluate_amp(&ichan1, &iint1, p.data(), amps.data());

#ifdef AC_FP_SINGLE
            const double tol = 1e-2;
#else
            const double tol = 1e-6;
#endif

#ifdef AC_FP_SINGLE
            const double abs_floor = 1e-10;
#else
            const double abs_floor = 0.0;
#endif

            for (int i = 0; i < namps; ++i) {
                double denom = std::abs(amps_save[i]) + std::abs(amps[i]);
                if (denom > abs_floor &&
                    std::abs(amps_save[i] - amps[i]) / denom > tol)
                {
                    fprintf(stderr,
                            "Process library not compatible with saved amplitudes"
                            " %d %d\n", ichan + 1, iint + 1);
                    fprintf(stderr, "  %d  (%g, %g)\n", i + 1,
                            amps_save[i].real(), amps_save[i].imag());
                    fprintf(stderr, "  %d  (%g, %g)\n", i + 1,
                            amps[i].real(), amps[i].imag());
                    std::exit(1);
                }
            }
        }
    }
}

// ------------------------------------------------------------------
// umami_matrix_element_cxx  (was umami_matrix_element_impl)
// ------------------------------------------------------------------

static void umami_matrix_element_cxx(
    int           ichan1,    // 1-based
    int           iint1,     // 1-based
    const double* p,         // momenta, shape (0:3, n) column-major
    double        alphas,
    const double* rnd,       // random numbers, length nvals
    double*       amp2,      // output, length nvals
    int*          helicities // output (1-based), length nvals
) {
    int ichan0 = ichan1 - 1;
    int iint0  = iint1  - 1;
    const ChannelGroup& cg = chans[ichan0];

    int namps_cur = cg.namps[iint0];

    // evaluate_amp always writes max_amps_global values; allocate accordingly to avoid overflow.
    std::vector<std::complex<double>> amps(static_cast<std::size_t>(max_amps_global));
    std::vector<double>               amp2_hel(static_cast<std::size_t>(namps_cur));

    evaluate_amp(&ichan1, &iint1, p, amps.data());

    if (keep_processes_separate) {
        // hel_fac layout: [ih + max_hel * iint0]
        for (int ih = 0; ih < namps_cur; ++ih) {
            double hf = cg.hel_fac[static_cast<std::size_t>(ih + cg.max_hel * iint0)];
            amp2_hel[ih] = std::real(amps[ih] * cg.col_fac[iint0] *
                                     std::conj(amps[ih])) * hf;
        }
        amp2[0] = 0.0;
        for (int ih = 0; ih < namps_cur; ++ih)
            amp2[0] += amp2_hel[ih];
    } else {
        for (int v = 0; v < cg.nvals; ++v) amp2[v] = 0.0;

        int iproc = 0;
        for (int ih = 0; ih < cg.namps[0]; ++ih) {
            while (cg.iproc_start[iproc + 1] == ih)
                ++iproc;
            double hf = cg.hel_fac[static_cast<std::size_t>(ih)]; // second dim is 1
            amp2_hel[ih] = std::real(amps[ih] * cg.col_fac[iproc] *
                                     std::conj(amps[ih])) * hf;
            amp2[iproc] += amp2_hel[ih];
        }
    }

    // Pick a random helicity for each process
    unwgt_helicities(ichan0, iint0, amp2_hel.data(), rnd, helicities);

    // Coupling prefactor
    double coupling = 1.0;
    int n_qcd = cg.n - 2 - cg.n_sing;
    if (cg.n_sing < cg.n - 2)
        coupling *= std::pow(4.0 * PI * alphas, static_cast<double>(n_qcd));
    if (cg.n_sing >= 1)
        coupling *= std::pow(2.0 * 4.0 * PI * alphaEW,
                             static_cast<double>(cg.n_sing));

    for (int v = 0; v < cg.nvals; ++v)
        amp2[v] *= coupling;

    // Identical-particle symmetry factors
    if (keep_processes_separate) {
        amp2[0] /= cg.iden[iint0];
    } else {
        for (int v = 0; v < cg.nvals; ++v)
            amp2[v] /= cg.iden[v];
    }
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// C interface  (was umami.c)
// ---------------------------------------------------------------------------

static pthread_mutex_t init_mutex   = PTHREAD_MUTEX_INITIALIZER;
static int             instance_count = 0;

struct InterfaceInstance {
    unsigned int magic;
    int     particle_count;
    double* momentum_buffer;
};

UmamiStatus umami_get_meta(UmamiMetaKey meta_key, void* result)
{
    switch (meta_key) {
        case UMAMI_META_DEVICE:
            *static_cast<int*>(result) = UMAMI_DEVICE_CPU;
            break;
        case UMAMI_META_DIAGRAM_COUNT:
            *static_cast<int*>(result) = 0;
            break;
        case UMAMI_META_PARTICLE_COUNT:
            if (nchans == 0) return UMAMI_ERROR;
            *static_cast<int*>(result) = chans[0].n;
            break;
        default:
            return UMAMI_ERROR_UNSUPPORTED_META;
    }
    return UMAMI_SUCCESS;
}

UmamiStatus umami_initialize(UmamiHandle* handle, char const* /*param_card_path*/)
{
    static bool data_loaded = false;

    pthread_mutex_lock(&init_mutex);
    if (!data_loaded) {
        try {
            read_amplitude_lib_light();
            if (!keep_processes_separate) {
                fprintf(stderr, "summing over processes currently not supported\n");
                pthread_mutex_unlock(&init_mutex);
                return UMAMI_ERROR;
            }
            test_library();
            data_loaded = true;
        } catch (const std::exception& e) {
            fprintf(stderr, "umami_initialize error: %s\n", e.what());
            pthread_mutex_unlock(&init_mutex);
            return UMAMI_ERROR;
        }
    }
    ++instance_count;
    pthread_mutex_unlock(&init_mutex);

    auto* inst          = new InterfaceInstance;
    inst->magic         = 0xAC1DC0DE;
    inst->particle_count = chans[0].n;
    inst->momentum_buffer =
        new double[static_cast<std::size_t>(4 * inst->particle_count)];
    *handle = inst;

    return UMAMI_SUCCESS;
}

UmamiStatus umami_set_parameter(
    UmamiHandle  /*handle*/,
    char const*  /*name*/,
    double       /*parameter_real*/,
    double       /*parameter_imag*/
) {
    return UMAMI_ERROR_NOT_IMPLEMENTED;
}

UmamiStatus umami_get_parameter(
    UmamiHandle  /*handle*/,
    char const*  /*name*/,
    double*      /*parameter_real*/,
    double*      /*parameter_imag*/
) {
    return UMAMI_ERROR_NOT_IMPLEMENTED;
}

UmamiStatus umami_matrix_element(
    UmamiHandle          handle,
    size_t               count,
    size_t               stride,
    size_t               offset,
    size_t               input_count,
    UmamiInputKey const* input_keys,
    void const* const*   inputs,
    size_t               output_count,
    UmamiOutputKey const* output_keys,
    void* const*         outputs
) {
    const double*        momenta_in        = nullptr;
    const double*        alpha_s_in        = nullptr;
    const unsigned int*  flavor_indices_in  = nullptr;
    const unsigned int*  channel_indices_in = nullptr;
    const double*        random_helicity_in = nullptr;

    for (size_t i = 0; i < input_count; ++i) {
        const void* input = inputs[i];
        switch (input_keys[i]) {
            case UMAMI_IN_MOMENTA:
                momenta_in = static_cast<const double*>(input);
                break;
            case UMAMI_IN_ALPHA_S:
                alpha_s_in = static_cast<const double*>(input);
                break;
            case UMAMI_IN_FLAVOR_INDEX:
                flavor_indices_in = static_cast<const unsigned int*>(input);
                break;
            case UMAMI_IN_RANDOM_HELICITY:
                random_helicity_in = static_cast<const double*>(input);
                break;
            case UMAMI_IN_CHANNEL_INDEX:
                channel_indices_in = static_cast<const unsigned int*>(input);
                break;
            default:
                return UMAMI_ERROR_UNSUPPORTED_INPUT;
        }
    }
    if (!momenta_in) return UMAMI_ERROR_MISSING_INPUT;

    double* m2_out       = nullptr;
    int*    helicity_out = nullptr;
    for (size_t i = 0; i < output_count; ++i) {
        void* output = outputs[i];
        switch (output_keys[i]) {
            case UMAMI_OUT_MATRIX_ELEMENT:
                m2_out = static_cast<double*>(output);
                break;
            case UMAMI_OUT_HELICITY_INDEX:
                helicity_out = static_cast<int*>(output);
                break;
            default:
                return UMAMI_ERROR_UNSUPPORTED_OUTPUT;
        }
    }

    auto* inst = static_cast<InterfaceInstance*>(handle);

    for (size_t i = 0; i < count; ++i) {
        int    ichan  = channel_indices_in ? static_cast<int>(channel_indices_in[i + offset]) + 1 : 1;
        int    iint   = flavor_indices_in  ? static_cast<int>(flavor_indices_in [i + offset]) + 1 : 1;
        double alphas = alpha_s_in         ? alpha_s_in        [i + offset] : 0.118;
        double rnd    = random_helicity_in ? random_helicity_in[i + offset] : 0.5;

        // Repack momenta from caller layout into column-major (0:3, n) Fortran layout
        for (int j = 0; j < inst->particle_count; ++j)
            for (int k = 0; k < 4; ++k)
                inst->momentum_buffer[k + 4 * j] =
                    momenta_in[i + offset + stride * (inst->particle_count * k + j)];

        double amp2    = 0.0;
        int    helicity = 0;

        umami_matrix_element_cxx(
            ichan,
            iint,
            inst->momentum_buffer,
            alphas,
            &rnd,
            &amp2,
            &helicity
        );

        if (m2_out)       m2_out      [i + offset] = amp2;
        if (helicity_out) helicity_out[i + offset] = helicity - 1; // back to 0-based
    }

    return UMAMI_SUCCESS;
}


UmamiStatus umami_free(UmamiHandle handle)
{
    if (!handle) return UMAMI_SUCCESS;   // guard against null

    auto* inst = static_cast<InterfaceInstance*>(handle);
    if(inst->magic != 0xAC1DC0DE) return UMAMI_SUCCESS; // guard against double-free
    inst->magic = 0; // invalidate magic number
    delete[] inst->momentum_buffer;
    inst->momentum_buffer = nullptr;
    delete inst;

    // NOTE: Deliberately do NOT clear the shared module state (`chans`,
    // `nchans`, ...) here. Those are namespace-scope C++ objects with
    // non-trivial destructors that the C++ runtime already finalizes once at
    // library unload.
    pthread_mutex_lock(&init_mutex);
    if (instance_count > 0) --instance_count;
    pthread_mutex_unlock(&init_mutex);
    return UMAMI_SUCCESS;
}
