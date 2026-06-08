#include "umami.h"
#include <pthread.h>
#include <stdlib.h>

pthread_mutex_t init_mutex = PTHREAD_MUTEX_INITIALIZER;
int instance_count = 0;

typedef struct {
    int particle_count;
    double* momentum_buffer;
} InterfaceInstance;

// defined in Fortran, see umami_impl.f03
void umami_initialize_impl();
void umami_free_impl();
void umami_matrix_element_impl(
    int const* ichan,
    int const* iint,
    double const* momenta,
    double const* alphas,
    double const* rnd,
    double* amp2,
    int* helicities
);
void get_particle_count(int* count);

UmamiStatus umami_get_meta(UmamiMetaKey meta_key, void* result) {
    switch (meta_key) {
        case UMAMI_META_DEVICE:
            *((int*)result) = UMAMI_DEVICE_CPU;
            break;
        case UMAMI_META_DIAGRAM_COUNT:
            *((int*)result) = 0;
            break;
        case UMAMI_META_PARTICLE_COUNT:
            get_particle_count((int*)result);
            break;
        default:
            return UMAMI_ERROR_UNSUPPORTED_META;
    }
    return UMAMI_SUCCESS;
}

UmamiStatus umami_initialize(UmamiHandle* handle, char const* param_card_path) {
    //TODO: allow for actual independent instances
    pthread_mutex_lock(&init_mutex);
    if (instance_count == 0) {
        umami_initialize_impl();
    }
    ++instance_count;
    pthread_mutex_unlock(&init_mutex);

    InterfaceInstance* inst = malloc(sizeof(InterfaceInstance));
    *handle = inst;
    get_particle_count(&inst->particle_count);
    inst->momentum_buffer = malloc(sizeof(double) * 4 * inst->particle_count);

    return UMAMI_SUCCESS;
}

UmamiStatus umami_set_parameter(
    UmamiHandle handle,
    char const* name,
    double parameter_real,
    double parameter_imag
) {
    return UMAMI_ERROR_NOT_IMPLEMENTED;
}

UmamiStatus umami_get_parameter(
    UmamiHandle handle,
    char const* name,
    double* parameter_real,
    double* parameter_imag
) {
    return UMAMI_ERROR_NOT_IMPLEMENTED;
}

UmamiStatus umami_matrix_element(
    UmamiHandle handle,
    size_t count,
    size_t stride,
    size_t offset,
    size_t input_count,
    UmamiInputKey const* input_keys,
    void const* const* inputs,
    size_t output_count,
    UmamiOutputKey const* output_keys,
    void* const* outputs
) {
    const double* momenta_in = NULL;
    const double* alpha_s_in = NULL;
    const unsigned int* flavor_indices_in = NULL;
    const unsigned int* channel_indices_in = NULL;
    const double* random_helicity_in = NULL;

    for (size_t i = 0; i < input_count; ++i) {
        const void* input = inputs[i];
        switch (input_keys[i]) {
            case UMAMI_IN_MOMENTA:
                momenta_in = (const double*)input;
                break;
            case UMAMI_IN_ALPHA_S:
                alpha_s_in = (const double*)input;
                break;
            case UMAMI_IN_FLAVOR_INDEX:
                flavor_indices_in = (const unsigned int*)input;
                break;
            case UMAMI_IN_RANDOM_HELICITY:
                random_helicity_in = (const double*)input;
                break;
            case UMAMI_IN_CHANNEL_INDEX:
                channel_indices_in = (const unsigned int*)input;
                break;
            default:
                return UMAMI_ERROR_UNSUPPORTED_INPUT;
        }
    }
    if (!momenta_in) return UMAMI_ERROR_MISSING_INPUT;

    double* m2_out = NULL;
    int* helicity_out = NULL;
    for (size_t i = 0; i < output_count; ++i) {
        void* output = outputs[i];
        switch (output_keys[i]) {
            case UMAMI_OUT_MATRIX_ELEMENT:
                m2_out = (double*)output;
                break;
            case UMAMI_OUT_HELICITY_INDEX:
                helicity_out = (int*)output;
                break;
            default:
                return UMAMI_ERROR_UNSUPPORTED_OUTPUT;
        }
    }

    InterfaceInstance* inst = (InterfaceInstance*)handle;
    for (size_t i = 0; i < count; ++i) {
        int ichan = channel_indices_in ? channel_indices_in[i + offset] + 1 : 1;
        int iint = flavor_indices_in ? flavor_indices_in[i + offset] + 1 : 1;
        double alphas = alpha_s_in ? alpha_s_in[i + offset] : 0.118;
        double rnd = random_helicity_in ? random_helicity_in[i + offset] : 0.5;
        double amp2;
        int helicity;
        for (int j = 0; j < inst->particle_count; ++j) {
            for (int k = 0; k < 4; ++k) {
                inst->momentum_buffer[k + 4 * j] =
                    momenta_in[i + offset + stride * (inst->particle_count * k + j)];
            }
        }
        umami_matrix_element_impl(
            &ichan,
            &iint,
            inst->momentum_buffer,
            &alphas,
            &rnd,
            &amp2,
            &helicity
        );
        if (m2_out) {
            m2_out[i + offset] = amp2;
        }
        if (helicity_out) {
            helicity_out[i + offset] = helicity - 1;
        }
    }

    return UMAMI_SUCCESS;
}

UmamiStatus umami_free(UmamiHandle handle)
{
    InterfaceInstance* inst = (InterfaceInstance*)handle;
    free(inst->momentum_buffer);
    free(handle);
    //TODO: allow for actual independent instances
    pthread_mutex_lock(&init_mutex);
    --instance_count;
    if (instance_count == 0) {
        umami_free_impl();
    }
    pthread_mutex_unlock(&init_mutex);
    return UMAMI_SUCCESS;
}
