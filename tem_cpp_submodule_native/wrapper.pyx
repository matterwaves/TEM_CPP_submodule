# distutils: language=c++

cdef extern from "wrapper.h":
    void generate_solvent_potential_extern(float* out_potential, int shape0, int shape1, int shape2, 
                      float* coords, void* proton_counts_raw, int atom_count, float r_asymptote,  
                      float r_probe, float pixel_size_rc, float pixel_size_z, float* vdw_dict, int vdw_dict_size,
                      int print_progress)

    void generate_seperated_points_extern(float* out, int element_num, float bound0, float bound1, float bound2, float min_atom_separation, int print_progress)

def generate_solvent_potential(
            bytes coords, 
            bytes proton_counts, 
            tuple pot_shape, 
            float r_asymptote, 
            float r_probe, 
            float pixel_size_rc, 
            float pixel_size_z,
            bytes van_der_walls_radii,
            int van_der_walls_radii_length,
            int print_progress
    ):
    assert len(pot_shape) == 3, "pot_shape must be a 3-tuple"
    assert len(coords) % (3 * sizeof(float)) == 0, "coords must be a multiple of 3 floats"
    assert len(proton_counts) % sizeof(int) == 0, "proton_counts must be a multiple of 1 int"
    assert len(coords) == len(proton_counts) * 3, "coords and proton_counts must have the same length"

    cdef bytes out_potential = bytes(pot_shape[0] * pot_shape[1] * pot_shape[2] * sizeof(float))
    cdef char* out_view = out_potential
    cdef char* coords_view = coords
    cdef char* proton_counts_view = proton_counts

    cdef char* vdw_dict = <char*>van_der_walls_radii

    generate_solvent_potential_extern(<float*>out_view, pot_shape[0], pot_shape[1], pot_shape[2],
                                      <float*>coords_view, <void*>proton_counts_view, len(proton_counts),
                                      r_asymptote, r_probe, pixel_size_rc, pixel_size_z, 
                                      <float*>vdw_dict, van_der_walls_radii_length, print_progress)
    
    return out_potential

def generate_seperated_points(int element_num, float bound0, float bound1, float bound2, float min_atom_separation, int print_progress):
    cdef bytes out_bytes = bytes(element_num * 3 * sizeof(float))
    cdef char* out_view = out_bytes
    generate_seperated_points_extern(<float*>out_view, element_num, bound0, bound1, bound2, min_atom_separation, print_progress)
    return out_bytes