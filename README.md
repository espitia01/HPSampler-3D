# Sampling the conformal space of an HP protein using Monte Carlo Methods
# Protein Folding Simulation in 3D HP Model

This repository contains a Python-based simulation for protein folding using the 3D Hydrophobic-Polar (HP) model. The HP model is a simplified representation of proteins where amino acids are classified as either hydrophobic (H) or polar (P). The simulation explores protein conformations in a 3D lattice and calculates the energy of each conformation based on the non-covalent interactions between hydrophobic residues.

## Programs

### `HPDistance.py`

This module calculates the cumulative sum of movements in a 3D lattice to determine the positions of amino acids in a protein chain. It generates a distance matrix representing the distances between all pairs of amino acids in the structure.

### `HPEnergy.py`

This module computes the energy of a given protein conformation within the 3D HP model. It uses a distance matrix to identify non-covalent contacts between hydrophobic residues that are adjacent but not sequentially connected in the protein chain.

### `HPShow.py`

This visualization tool uses `matplotlib` to plot the 3D structure of the protein. It represents hydrophobic residues in blue, polar residues in red, and marks the first amino acid in green. It also displays the current energy, temperature, and time of the simulation.

### `HPFold.py`

The main simulation engine that iteratively explores different protein conformations by performing Monte Carlo moves. It uses simulated annealing, varying the temperature to escape local minima and find a conformation with minimal energy. It records the simulation data, including energy levels and structures, and outputs this data to a CSV file.

### `HPMove.py` (Not provided, but mentioned)

Presumably, this module is responsible for generating new conformations by making random moves in the protein structure, which are then evaluated by the simulation in `HPFold.py`.

## Usage

To run the simulation, ensure you have Python and the required libraries installed, then execute the `HPFold.py` script with the desired parameters.
