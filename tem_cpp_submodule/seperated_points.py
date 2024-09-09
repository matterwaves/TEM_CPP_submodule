import tem_cpp_submodule_native
import numpy as np

def generate_seperated_points(n_points, x_dim, y_dim, z_dim, min_separation):
    """
    Generates a set of separated points in 3D space.
    Parameters:
    - n_points (int): The number of points to generate.
    - x_dim (float): The maximum value for the x-coordinate of the points.
    - y_dim (float): The maximum value for the y-coordinate of the points.
    - z_dim (float): The maximum value for the z-coordinate of the points.
    - min_separation (float): The minimum separation between points.
    Returns:
    - numpy.ndarray: An array of shape (n_points, 3) containing the generated points.
    Example:
    >>> generate_seperated_points(100, 10.0, 10.0, 10.0, 1.0)
    array([[ 0.5,  1.2,  3.4],
           [ 2.3,  4.5,  6.7],
           [ 8.9,  7.6,  5.4],
           ...
           [ 9.1,  8.2,  7.3]])
    """
    
    coords_bytes = tem_cpp_submodule_native.generate_seperated_points(n_points, x_dim, y_dim, z_dim, min_separation)
    return np.frombuffer(coords_bytes, dtype=np.float32).reshape((n_points, 3))