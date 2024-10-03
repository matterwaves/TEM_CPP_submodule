import tem_cpp_submodule_native
import numpy as np
import typing

van_der_walls_radii_default = [
    -1.0,  1.1, -1.0, -1.0, -1.0, -1.0,  1.7, 1.55, 
    1.52, -1.0, -1.0, 2.27, 1.73, -1.0, -1.0,  1.8, 
     1.8, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.66
]

def generate_solvent_potential(
        coords: np.ndarray, 
        proton_counts: np.ndarray, 
        potential_shape: typing.Tuple[int, int, int],
        pix_size_xy: float,
        pix_size_z: float,
        r_asymptote: float = 7.5, # in angstrom
        r_probe: float = 1.7, # in angstrom
        van_der_walls_radii: typing.List[float] = van_der_walls_radii_default,
        print_progress: bool = False
        ) -> np.ndarray:
    """
    Generates a solvent potential grid based on atomic coordinates and proton counts.

    Args:
        coords (np.ndarray): A 2D array of shape (N, 3), where N is the number of atoms,
            and each row represents the x, y, z coordinates of an atom.
        proton_counts (np.ndarray): A 1D array of shape (N,), where each element represents
            the proton count for the corresponding atom in `coords`.
        potential_shape (Tuple[int, int, int]): The shape of the output potential grid as 
            (shape_x, shape_y, shape_z).
        pix_size_xy (float): The pixel size in the x and y dimensions (in angstroms).
        pix_size_z (float): The pixel size in the z dimension (in angstroms).
        r_asymptote (float, optional): The maximum distance at which an atom contributes to the solvent potential, default is 7.5 angstroms.
        r_probe (float, optional): The probe radius when generating the solvent mask (as described by Shang and Sigworth), default is 1.7 angstroms.
        van_der_walls_radii (List[float], optional): A list of van der Waals radii for each atom type.
        print_progress (bool, optional): If True, prints progress information during the computation.

    Returns:
        np.ndarray: A 3D array of shape `potential_shape` representing the computed solvent potential.
    """

    assert coords.dtype == np.float32
    assert proton_counts.dtype == np.uint32

    assert len(coords.shape) == 2
    assert coords.shape[1] == 3
    assert len(proton_counts.shape) == 1

    assert coords.shape[0] == proton_counts.shape[0]

    assert len(potential_shape) == 3

    coords_bytes = coords.tobytes()
    proton_counts_bytes = proton_counts.tobytes()
    van_der_walls_radii_bytes = np.array(van_der_walls_radii, dtype=np.float32).tobytes()

    potential_bytes = tem_cpp_submodule_native.generate_solvent_potential(
        coords_bytes, 
        proton_counts_bytes, 
        potential_shape, 
        r_asymptote, 
        r_probe, 
        pix_size_xy, 
        pix_size_z, 
        van_der_walls_radii_bytes, 
        len(van_der_walls_radii), 
        1 if print_progress else 0
    )

    return np.frombuffer(potential_bytes, dtype=np.float32).reshape(potential_shape)

