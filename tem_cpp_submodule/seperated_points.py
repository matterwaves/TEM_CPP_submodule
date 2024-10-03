import tem_cpp_submodule_native
import numpy as np

def generate_seperated_points(n_points: int, x_max: float, y_max: float, z_max: float, min_separation: float, print_progress: bool = False) -> np.ndarray:
    """
    Generates a set of separated points in 3D space.

    Args:
       n_points (int): The number of points to generate.
       x_dim (float): The maximum value for the x-coordinate of the points.
       y_dim (float): The maximum value for the y-coordinate of the points.
       z_dim (float): The maximum value for the z-coordinate of the points.
       min_separation (float): The minimum separation between points.

    Returns:
       np.ndarray: A 2D array of shape (n_points, 3) containing the generated points.
    """
    
    coords_bytes = tem_cpp_submodule_native.generate_seperated_points(n_points, x_max, y_max, z_max, min_separation, 1 if print_progress else 0)
    return np.frombuffer(coords_bytes, dtype=np.float32).reshape((n_points, 3))