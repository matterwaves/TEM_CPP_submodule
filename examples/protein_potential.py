import tem_cpp_submodule
import numpy as np

from matplotlib import pyplot as plt

coords = np.array([
    [-2, -2, 6],
    [6, 9, 6],
    [10, 10, 6]
], dtype=np.float32)

proton_counts = np.array([6, 7, 8], dtype=np.uint32)


solvent_potential = tem_cpp_submodule.generate_protein_potential(
    coords, 
    proton_counts, 
    (20, 20, 20), 1, 1
)

plt.imshow(solvent_potential.sum(axis=2))
plt.show()
