import tem_cpp_submodule
import numpy as np
import os

from matplotlib import pyplot as plt

from Bio import PDB

proton_count_dict = {'H': 1, 'C': 6, 'N': 7, 'O': 8,
          'NA': 11, 'MG': 12, 'P': 15, 'S': 16,
          'FE': 26, 'ZN': 30, 'AU': 79}

def read_pdb(filename):
    # create file parser
    
    parser = PDB.PDBParser(QUIET=True)

    # read structure
    structure = parser.get_structure('struct', filename)
    
    coords = []
    proton_counts = []

    # get atom data from file
    for _, atom in enumerate(structure.get_atoms()):
        if atom.element == "C" or atom.element == "N" or atom.element == "O":
            raw_coord = atom.get_coord()
            coords.append([-raw_coord[1], raw_coord[0], -raw_coord[2]])
            proton_counts.append(proton_count_dict[atom.element])

    coords = np.array(coords)
    coords -= coords.mean(axis=0)
    
    return np.array(coords).astype(np.float32), np.array(proton_counts).astype(np.uint32)

my_file_directory = os.path.dirname(os.path.realpath(__file__))

coords, proton_counts = read_pdb(os.path.join(my_file_directory, "sample_ribosome_structure.pdb"))

potential_shape = (512, 512, 512)

protein_potential = tem_cpp_submodule.generate_protein_potential(
    coords, 
    proton_counts,
    potential_shape,
    1, 1,
    print_progress=True
)

solvent_potential = tem_cpp_submodule.generate_solvent_potential(
    coords, 
    proton_counts,
    potential_shape,
    1, 1,
    print_progress=True
)

full_potential = solvent_potential + protein_potential

plt.imshow(full_potential.sum(axis=0))
plt.show()