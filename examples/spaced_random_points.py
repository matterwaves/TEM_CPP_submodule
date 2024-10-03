import tem_cpp_submodule
import numpy as np

n_avogadro = 6.02214076e23

def generate_carbon(vol_dims_xyz, min_atom_separation, atom_density_gcm3):
    atomic_weight = 12 # approximation!! reasonable for C
    subst_vol = np.prod(vol_dims_xyz) # [A^3]
    n_atoms = int(subst_vol * atom_density_gcm3 * 1e-24 * n_avogadro / atomic_weight)

    return tem_cpp_submodule.generate_seperated_points(n_atoms, *vol_dims_xyz, min_atom_separation)

print(generate_carbon((500, 500, 500), 1.5, 1.0).shape)
