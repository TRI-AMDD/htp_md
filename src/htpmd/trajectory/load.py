"""
Module for trajectory loaders.
"""
from abc import ABC
from abc import abstractmethod
import os
import glob
import json
import numpy as np

from htpmd.trajectory.base import ExtendedLAMMPSTrajectoryFile
from htpmd.shared.population_matrix import generate_population_matrix, compute_pop_matrix_from_pizza_dump

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
    'salt_cation',
    'salt_anion',
}


class TrajectoryLoader(ABC):

    @abstractmethod
    def load(self, source):
        raise NotImplementedError


class LammpsTrajectoryLoader(TrajectoryLoader):

    def load(self, trajectory_path):
        traj_file = glob.glob(trajectory_path+'/'+'*.lammpstrj')
        file_name = traj_file[0]
        trajectory = ExtendedLAMMPSTrajectoryFile(file_name)
        trajectory.read()
        trajectory.determine_lattice_and_coords()
        return trajectory


def get_metadata(dir_name):
    metadata_path = os.path.join(dir_name, 'meta.json')
    with open(metadata_path) as f:
        metadata = json.load(f)
    if not set(metadata.keys()).issuperset(REQUIRED_METADATA):
        fields_missing = [x for x in REQUIRED_METADATA if x not in metadata.keys()]
        fields_missing_str = ','.join(fields_missing)
        error_message = f'The following fields are missing in metadata file {metadata_path} : {fields_missing_str}'
        raise ValueError(error_message)
    return metadata


def get_population_matrix(dir_name):
    """Loads the population matrix computed by lammps if population.txt is available; otherwise, computes population
    matrix if aggregates.txt file is available. Computes population matrix using generate_population_matrix if non of
    population.txt and aggregates.txt """
    file_ref_path = os.path.join(dir_name, 'population.txt')
    path_to_aggregate_file = os.path.join(dir_name, 'aggregates.txt')
    if os.path.exists(file_ref_path):
        pop_mat = np.loadtxt(file_ref_path)
        stacked_pop_mat = np.empty((0, 1))
    elif os.path.exists(path_to_aggregate_file):
        pop_mat = compute_pop_matrix_from_pizza_dump(dir_name)
        np.savetxt(os.path.join(dir_name, 'population.txt'), pop_mat)
        stacked_pop_mat = np.empty((0, 1))
    else:
        stacked_pop_mat, pop_mat = generate_population_matrix(dir_name)
        # save the population matrices in the local directory
        np.savetxt(os.path.join(dir_name, 'population.txt'), pop_mat)
        np.save(os.path.join(dir_name, 'stacked_pop'), stacked_pop_mat)

    return stacked_pop_mat, pop_mat
