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

void calc_pot(struct arrayF3 out_pot, float* coords, void* proton_counts, int atom_count,
              struct arrayF3 atom_pots, float pixel_size_rc, float pixel_size_z);
   
void make_solvent_pot(float* dist_map, int shape0, int shape1, int shape2,
                      float* coords, void* proton_counts, int atom_count, 
                      float r_asymptote, float r_probe, float pixel_size_rc, float pixel_size_z);

 
void generate_seperated_points(float* out, int element_num, float bound0, float bound1, float bound2, float min_atom_separation);

#endif