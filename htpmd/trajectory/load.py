"""
Module for trajectory loaders.
"""
from abc import ABC
from abc import abstractmethod
import os
import json
import numpy as np

from htpmd.trajectory.base import ExtendedLAMMPSTrajectoryFile


REQUIRED_METADATA = {
    'mol_smiles',  # Smiles for molecule. E.g. '[Cu]CCO[Au]' for PEO
    'poly_smiles',  # Smiles for polymer.  E.g. 'CCOCCOCCOCCOCCOCCO' for PEO
    'force_field',  # Force field. E.g. 'PCFF+'
    'material_group',  # Type of materials. E.g. polymer
    'temperature',  # Temperature in K. E.g. 353
    'time_step',  # Time step between saved trajectory frames in ps. E.g. 2.
    'cation_raw_type',  # Raw atom type in LAMMPS used to identify cation diffusivity. E.g. 90
    'anion_raw_type',  # Raw atom type in LAMMPS used to identfiy anion diffusivity. E.g. 93
    'polymer_raw_type_range',  # The range of raw atom types of all polymer atoms. E.g. [0, 90]
    'polymer_solvate_types',  # A list of atom types used to define polymer diffusivity. E.g. [7, 8, 16]
}


class TrajectoryLoader(ABC):

    @abstractmethod
    def load(self, source):
        raise NotImplementedError


class LammpsTrajectoryLoader(TrajectoryLoader):

    def load(self, trajectory_path):
        file_name = os.path.join(trajectory_path, 'traj.lammpstrj')
        trajectory = ExtendedLAMMPSTrajectoryFile(file_name)
        trajectory.read()
        trajectory.determine_lattice_and_coords()
        return trajectory


def get_metadata(dir_name):
    with open(os.path.join(dir_name, 'meta.json')) as f:
        metadata = json.load(f)
    assert set(metadata.keys()).issuperset(REQUIRED_METADATA)
    return metadata


def get_population_matrix(dir_name):
    """Load the population matrix computed by lammps."""
    pop_mat = np.loadtxt(os.path.join(dir_name, 'population.txt'))
    return pop_mat
