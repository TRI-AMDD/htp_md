"""
Module for trajectory loaders.
"""
from abc import ABC
from abc import abstractmethod
import os

from htpmd.trajectory.base import ExtendedLAMMPSTrajectoryFile


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
