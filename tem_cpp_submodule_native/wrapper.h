#ifndef __SRC__ALGS_H
#define __SRC__ALGS_H

#include <stdint.h>

struct arrayF3 {
    union {
        void* vbuff;
        float* fbuff;
    };
    int shape0;
    int shape1;
    int shape2;
};

// void calc_pot(struct arrayF3 out_pot, float* coords, void* proton_counts, int atom_count,
//               struct arrayF3 atom_pots, float pixel_size_rc, float pixel_size_z);
   
void generate_solvent_potential_extern(float* out_potential, int shape0, int shape1, int shape2, 
                      float* coords, void* proton_counts_raw, int atom_count, float r_asymptote,  
                      float r_probe, float pixel_size_rc, float pixel_size_z, float* vdw_dict, int vdw_dict_size,
                      int print_progress);
 
void generate_seperated_points_extern(float* out, int element_num, float bound0, float bound1, float bound2, float min_atom_separation, int print_progress);

#endif