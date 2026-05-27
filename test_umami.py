import madspace as ms
import numpy as np

np.random.seed(12345)

me = ms.default_context().load_matrix_element("libamplicolmadspace.so", "")

n = 100

psmap = ms.PhaseSpaceMapping([0., 0., 0., 0., 0.], 13000.)
p_ext, x1, x2, det = psmap.map_forward([np.random.rand(n, psmap.random_dim())])
alpha_s = 0.118 * np.ones(n)
flavors = np.zeros(n, dtype=np.int32)
channels = np.zeros(n, dtype=np.int32)
rnd_hel = np.random.rand(n)

func = ms.MatrixElement(
    me,
    [
        ms.MatrixElement.momenta_in,
        ms.MatrixElement.alpha_s_in,
        ms.MatrixElement.flavor_in,
        ms.MatrixElement.channel_in,
        ms.MatrixElement.random_helicity_in
    ],
    [ms.MatrixElement.matrix_element_out, ms.MatrixElement.helicity_index_out]
)
amp2, hel = func(p_ext, alpha_s, flavors, channels, rnd_hel)
print(amp2)
print(hel)
