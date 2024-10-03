#ifndef __SRC__ALGS_H
#define __SRC__ALGS_H

#include <stdint.h>

struct arrayF3_old {
    union {
        void* vbuff;
        float* fbuff;
    };
    int shape0;
    int shape1;
    int shape2;
};

void calc_pot(
    float* out_pot, 
    int out_pot_shape0,
    int out_pot_shape1,
    int out_pot_shape2,
    float* coords, 
    void* proton_counts_raw, 
    int atom_count,
    float* atom_pots,
    int atom_pots_shape0,
    int atom_pots_shape1,
    int atom_pots_shape2,
    float pixel_size_rc,
    float pixel_size_z,
    int print_progress);

void generate_solvent_potential_extern(float* out_potential, int shape0, int shape1, int shape2, 
                      float* coords, void* proton_counts_raw, int atom_count, float r_asymptote,  
                      float r_probe, float pixel_size_rc, float pixel_size_z, float* vdw_dict, int vdw_dict_size,
                      int print_progress);
 
void generate_seperated_points_extern(float* out, int element_num, float bound0, float bound1, float bound2, float min_atom_separation, int print_progress);

#endif