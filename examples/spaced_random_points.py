import tem_cpp_submodule_native
import numpy as np

n_avogadro = 6.02214076e23

def generate_carbon(vol_dims_xyz, min_atom_separation, atom_density_gcm3):
    atomic_weight = 12 # approximation!! reasonable for C
    subst_vol = np.prod(vol_dims_xyz) # [A^3]
    n_atoms = int(subst_vol * atom_density_gcm3 * 1e-24 * n_avogadro / atomic_weight)

    coords_bytes = tem_cpp_submodule_native.generate_seperated_points(n_atoms, vol_dims_xyz[0], vol_dims_xyz[1], vol_dims_xyz[2], min_atom_separation)

    return np.frombuffer(coords_bytes, dtype=np.float32).reshape((n_atoms, 3))

print(generate_carbon((1000, 1000, 1000), 1.5, 2.0).shape)
