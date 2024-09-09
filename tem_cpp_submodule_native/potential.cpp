#include "wrapper.h"
#include <stdio.h>
#include <math.h>

static int min_c(int a, int b) {
    if(a < b) return a;
    return b;
}

static int max_c(int a, int b) {
    if(a > b) return a;
    return b;
}

// Function to calculate 3d potential from atom list
void calc_pot(struct arrayF3 out_pot, float* coords, void* proton_counts_raw, int atom_count,
              arrayF3 atom_pots, float pixel_size_rc, float pixel_size_z, int print_progress)
{
    // Convert pointers (python doesn't like having pointers to non standard types)
    int32_t* proton_counts = (int32_t*)proton_counts_raw;

    // Loop through atoms
    for(int i = 0; i < atom_count; i++) {
        // Get atom id and check that it is within range
        int atom_id = proton_counts[i];
        if(atom_id < 0 || atom_id > atom_pots.shape0) continue;

        // Get the potential of this atom
        float* atom_pot = &atom_pots.fbuff[atom_id * atom_pots.shape1 * atom_pots.shape2];

        // The atom positions
        float pos0 = coords[3*i + 0];
        float pos1 = coords[3*i + 1];
        float pos2 = coords[3*i + 2];

        // The indicies of the atom positions in the potential array
        int ind0 = (int)roundf(pos0/pixel_size_rc + out_pot.shape0/2.0f);
        int ind1 = (int)roundf(pos1/pixel_size_rc + out_pot.shape1/2.0f);
        int ind2 = (int)roundf(pos2/pixel_size_z  + out_pot.shape2/2.0f);
        
        if(ind2 < 0 || ind2 >= out_pot.shape2) continue;

        int start0 = ind0 - atom_pots.shape1/2;
        int start1 = ind1 - atom_pots.shape2/2;

        // These are the bounding indicies of where we will do our computation
        int s0 = max_c(start0, 0);
        int e0 = min_c(start0 + atom_pots.shape1, out_pot.shape0);

        int s1 = max_c(start1, 0);
        int e1 = min_c(start1 + atom_pots.shape2, out_pot.shape1);

        // Embed atom potential into the 3D potential
        for(int i0 = s0; i0 < e0; i0++) {
            int my_ind0 = i0 - start0;
            for(int i1 = s1; i1 < e1; i1++) {
                int my_ind1 = i1 - start1;
                out_pot.fbuff[(i0*out_pot.shape1 + i1)*out_pot.shape2 + ind2] += atom_pot[my_ind0 * atom_pots.shape2 + my_ind1];
            }
        }

        if(print_progress != 0 && i % 10000 == 0)
            printf("\r... ... Embedded %d/%d atoms", i, atom_count);
    }

    if(print_progress != 0)
        printf("\n");
}
