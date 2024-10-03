#include "wrapper.h"

#include <cmath>
#include <vector>

static int min_c(int a, int b) {
    if(a < b) return a;
    return b;
}

static int max_c(int a, int b) {
    if(a > b) return a;
    return b;
}

// Shang and Sigworth Density (Polar)
float ssdPolar(float rad_dist) {
    float num1 =  0.5f  * erff( (rad_dist - 0.5f) / sqrtf(2));
    float num2 =  0.2f  * expf(-(rad_dist - 1.7f)*(rad_dist - 1.7f) / (2 * 1.77f * 1.77f));
    float num3 = -0.15f * expf(-(rad_dist - 1.7f)*(rad_dist - 1.7f) / (2 * 1.06f * 1.06f));
    return 0.5 + num1 + num2 + num3;
}

// Shang and Sigworth Density (Non-Polar)
float ssdNonPolar(float rad_dist) {
    float num1 =  0.5f  * erff( (rad_dist - 1) / sqrtf(2));
    float num2 =  0.15f * expf(-(rad_dist - 2.2f)*(rad_dist - 2.2f) / (2 * 1.77f * 1.77f));
    float num3 = -0.12f * expf(-(rad_dist - 3.6f)*(rad_dist - 3.6f) / (2 * 0.85f * 0.85f));
    return 0.5 + num1 + num2 + num3;
}

 
void generate_solvent_potential_extern(float* out_potential, int shape0, int shape1, int shape2, 
                      float* coords, void* proton_counts_raw, int atom_count, float r_asymptote,  
                      float r_probe, float pixel_size_rc, float pixel_size_z, float* vdw_dict, int vdw_dict_size,
                      int print_progress)
{
    printf("Van der wall radius of atom 0: %f %d\n", vdw_dict[0], vdw_dict_size);

    // convert pointers (python doesn't like having pointers to non standard types)
    int32_t* proton_counts = (int32_t*)proton_counts_raw;
    
    std::vector<bool> solvent_mask;
    solvent_mask.resize(shape0 * shape1 * shape2);
    std::fill(solvent_mask.begin(), solvent_mask.end(), true);

    std::vector<int8_t> proton_map;
    proton_map.resize(shape0 * shape1 * shape2);

    float min_dist = 0;

    // Now we loop over the atoms in the list to create the initial masks
    for(int i = 0; i < atom_count; i++) {
        // Skip any atoms with proton counts that are negative or havier than gold
        int atom_id = proton_counts[i];
        if(atom_id < 0 || atom_id > vdw_dict_size - 1) continue;

        // `vdw_rad` is the actuall vad der walls radius of the atom
        float vdw_rad = vdw_dict[atom_id]; 
        if(vdw_rad < 0) continue; // Skip atoms that we don't support

        // r_asymptote in pixel units
        int rsrc = (int)r_asymptote/pixel_size_rc;
        int rsz  = (int)r_asymptote/pixel_size_z;

        // The atom positions
        float pos0 = coords[3*i + 0];
        float pos1 = coords[3*i + 1];
        float pos2 = coords[3*i + 2];

        // The indicies of the atom positions in the potential array
        int ind0 = (int)ceil(pos0/pixel_size_rc + shape0/2.0f);
        int ind1 = (int)ceil(pos1/pixel_size_rc + shape1/2.0f);
        int ind2 = (int)ceil(pos2/pixel_size_z  + shape2/2.0f);

        // These are the bounding indicies of where we will do our computation
        // atom_position.xyz +- r_asymptote.xyz
        int s0 = max_c(ind0 - rsrc, 0);
        int e0 = min_c(ind0 + rsrc, shape0);

        int s1 = max_c(ind1 - rsrc, 0);
        int e1 = min_c(ind1 + rsrc, shape1);

        int s2 = max_c(ind2 - rsz, 0);
        int e2 = min_c(ind2 + rsz, shape2);

        // Actually perform the per atom calculation
        for(int i0 = s0; i0 < e0; i0++) {
            // c0 = x position of atom relative to this pixel
            float c0 = (i0 - shape0/2)*pixel_size_rc - pos0;
            for(int i1 = s1; i1 < e1; i1++) {
                // c1 = y position of atom relative to this pixel
                float c1 = (i1 - shape1/2)*pixel_size_rc - pos1;
                for(int i2 = s2; i2 < e2; i2++) {
                    // c2 = z position of atom relative to this pixel
                    float c2 = (i2 - shape2/2)*pixel_size_z - pos2;

                    float curr_r = sqrtf(c0*c0 + c1*c1 + c2*c2) - vdw_rad;

                    float curr_dist = out_potential[(i0*shape1 + i1)*shape2 + i2];

                    if(curr_dist > curr_r) {
                        out_potential[(i0*shape1 + i1)*shape2 + i2] = curr_r;
                        proton_map[(i0*shape1 + i1)*shape2 + i2] = atom_id;

                        if(min_dist > curr_r) min_dist = curr_r;
                    }

                    if(r_probe > curr_r)
                        solvent_mask[(i0*shape1 + i1)*shape2 + i2] = false;
                }
            }
        }

        if(print_progress == 1 && i % 10000 == 0) printf("\rAnalyzed %d/%d atoms", i, atom_count);
    }

    if(print_progress == 1)
        printf("\n");

    // define a custom Index3D type just for convenience
    struct Index3D{
        int32_t i0, i1, i2;
    };

    // Now we find the boundary points
    std::vector<struct Index3D> boundary_inds;
    
    // Go through mask to find all boundary points
    for(int i0 = 1; i0 < shape0-1; i0++) {
        for(int i1 = 1; i1 < shape1-1; i1++) {
            for(int i2 = 1; i2 < shape2-1; i2++) {
                // If the mask is equal to zero, we skip it
                if(solvent_mask[(i0*shape1 + i1)*shape2 + i2] == false) continue;

                // Check if any neighboring pixels are equal to zero
                if(solvent_mask[((i0 - 1)*shape1 + i1)*shape2 + i2] == false ||
                   solvent_mask[((i0 + 1)*shape1 + i1)*shape2 + i2] == false ||
                   solvent_mask[(i0*shape1 + (i1 - 1))*shape2 + i2] == false ||
                   solvent_mask[(i0*shape1 + (i1 + 1))*shape2 + i2] == false ||
                   solvent_mask[(i0*shape1 + i1)*shape2 + (i2 - 1)] == false ||
                   solvent_mask[(i0*shape1 + i1)*shape2 + (i2 + 1)] == false)
                {
                    struct Index3D ind;
                    ind.i0 = i0;
                    ind.i1 = i1;
                    ind.i2 = i2;

                    // if so, we add these coords to the vector
                    boundary_inds.push_back(ind);
                }
            }
        }
    }

    // Finally we fill in the areas around the boundary points
    int i = 0;
    for(struct Index3D ind: boundary_inds) {
        int ind0 = ind.i0;
        int ind1 = ind.i1;
        int ind2 = ind.i2;

        // Get atom id and check that it's within bounds
        int atom_id = proton_map[(ind0*shape1 + ind1)*shape2 + ind2];
        if(atom_id < 1 || atom_id > vdw_dict_size-1) continue;

        // Get van der wals radius, and skip if it isn't defined
        float vdw_rad = vdw_dict[atom_id];
        if(vdw_rad < 0) continue;
        
        float r_shrink = vdw_rad + out_potential[(ind0*shape1 + ind1)*shape2 + ind2];

        // r_shrink in units of pixels instead of angstroms
        int rsrc = (int)(r_shrink/pixel_size_rc);
        int rsz  = (int)(r_shrink/pixel_size_z);

        // factors to convert from units of pixels^2 to angrom^2
        float prc = pixel_size_rc*pixel_size_rc;
        float pz = pixel_size_z*pixel_size_z;

        // precompute r_shrink^2 just because
        float rs2 = r_shrink*r_shrink;

        // These are the bounding indicies of where we will do our computation
        // boundary_point.xyz +- r_shrink.xyz
        int s0 = max_c(ind0 - rsrc - 1, 0);
        int e0 = min_c(ind0 + rsrc + 1, shape0);

        int s1 = max_c(ind1 - rsrc - 1, 0);
        int e1 = min_c(ind1 + rsrc + 1, shape1);

        int s2 = max_c(ind2 - rsz - 1, 0);
        int e2 = min_c(ind2 + rsz + 1, shape2);

        // Actually perform the per boundary point calculation
        for(int i0 = s0; i0 < e0; i0++) {
            // dd0: dist to boundary point squared in pixels^2 (only on axis0)
            int dd0 = (i0 - ind0)*(i0 - ind0);
            for(int i1 = s1; i1 < e1; i1++) {
                // dd1: dist to boundary point squared in pixels^2 (axis0 + axis1)
                int dd1 = dd0 + (i1 - ind1)*(i1 - ind1);
                for(int i2 = s2; i2 < e2; i2++) {
                    // check if dist is smaller than r_shrink, if so set mask to 1
                    if(dd1*prc + (i2 - ind2)*(i2 - ind2)*pz < rs2)
                        solvent_mask[(i0*shape1 + i1)*shape2 + i2] = true;
                }
            }
        }

        if(print_progress == 1 && i % 10000 == 0) printf("\rAnalyzed %d/%ld boundary points", i, boundary_inds.size());
        i++;
    }

    if(print_progress == 1)
        printf("\n");

    const int SAMPLES = 10000;
    float dist_range = r_asymptote - min_dist;

    std::vector<float> polar_potential;
    polar_potential.resize(SAMPLES+1);
    
    std::vector<float> nonpolar_potential;
    nonpolar_potential.resize(SAMPLES+1);

    for(int i = 0; i < SAMPLES+1; i++) {
        float i_dist = min_dist + dist_range*i/((float)SAMPLES);
        polar_potential[i] = 3.6f * pixel_size_z * ssdPolar(i_dist);
        nonpolar_potential[i] = 3.6f * pixel_size_z * ssdNonPolar(i_dist);
    }
    
    if(print_progress == 1)
        printf("Calculating final potential...");
    
    // Go through mask to find all boundary points
    for(int i0 = 0; i0 < shape0; i0++) {
        for(int i1 = 0; i1 < shape1; i1++) {
            for(int i2 = 0; i2 < shape2; i2++) {
                int index = (i0*shape1 + i1)*shape2 + i2;

                // If the mask is equal to zero, we skip it
                if(solvent_mask[index] == false) {
                    out_potential[index] = 0;
                    continue;
                }

                float pot_index_f = SAMPLES*(out_potential[index] - min_dist)/dist_range;

                int pot_index_l = (int)floorf(pot_index_f);
                int pot_index_h = (int) ceilf(pot_index_f);
                
                if(pot_index_l < 0) pot_index_l = 0;
                if(pot_index_h < 0) pot_index_h = 0;
                
                if(pot_index_l > SAMPLES) pot_index_l = SAMPLES;
                if(pot_index_h > SAMPLES) pot_index_h = SAMPLES;

                float lerp_val = pot_index_f - pot_index_l;
                
                if(proton_map[index] == 6)
                    out_potential[index] = (1 - lerp_val) * nonpolar_potential[pot_index_l] + lerp_val * nonpolar_potential[pot_index_h];
                else
                    out_potential[index] = (1 - lerp_val) * polar_potential[pot_index_l] + lerp_val * polar_potential[pot_index_h];
            }
        }
    }
    
    if(print_progress == 1)
        printf("done!\n");
}