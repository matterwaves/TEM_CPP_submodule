import tem_cpp_submodule_native
import numpy as np
import typing

def generate_solvent_potential(
        coords: np.ndarray, 
        proton_counts: np.ndarray, 
        potential_shape: typing.Tuple[int, int, int],
        r_asymptote: float,
        r_probe: float,
        pix_size_rc: float,
        pix_size_z: float
        ) -> np.ndarray:
    pass
