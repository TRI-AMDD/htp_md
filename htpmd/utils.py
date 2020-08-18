from __future__ import print_function, division

import os
import json
import numpy as np
from collections import defaultdict
from mdtraj.formats import LAMMPSTrajectoryFile
from mdtraj.formats.lammpstrj import _EOF
from mdtraj.utils.six.moves import xrange
from scipy.spatial import KDTree, cKDTree
import itertools
import heapq


_TOL = 1e-2
_ATOM_MASSES = [
    0.0, 1.008, 4.002602, 6.94, 9.0121831, 10.81, 12.011, 14.007, 15.999,
    18.998403163, 20.1797, 22.98976928, 24.305, 26.9815385, 28.085,
    30.973761998, 32.06, 35.45, 39.948, 39.0983, 40.078, 44.955908,
    47.867, 50.9415, 51.9961, 54.938044, 55.845, 58.933194, 58.6934, 63.546,
    65.38, 69.723, 72.63, 74.921595, 78.971, 79.904, 83.798, 85.4678, 87.62,
    88.90584, 91.224, 92.90637, 95.95, 97.90721, 101.07, 102.9055, 106.42,
    107.8682, 112.414, 114.818, 118.71, 121.76, 127.6, 126.90447, 131.293,
    132.90545196, 137.327, 138.90547, 140.116, 140.90766, 144.242, 144.91276,
    150.36, 151.964, 157.25, 158.92535, 162.5, 164.93033, 167.259, 168.93422,
    173.045, 174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217,
    195.084, 196.966569, 200.592, 204.38, 207.2, 208.9804, 209.0, 210.0,
    222.0, 223.0, 226.0, 227.0, 232.0377, 231.03588, 238.02891, 237.0, 244.0,
    243.0, 247.0, 247.0, 251.0, 252.0]


def load_lammps(lammps_file, use_mass=False, tol=_TOL):
    """Load a lammpstraj file.

    Args:
        lammps_file: file path
        use_mass: whether to use mass to determine atom types

    Returns:
        coords: shape (F, N, 3)
        lattices: shape (F, 3)
        types: shape (N,)
        unwrapped_coords: shape (F, N, 3)
    """
    # Loads lammpstraj file.
    with ExtendedLAMMPSTrajectoryFile(lammps_file) as f:
        coords, lattices, angles, types, ixyz, masses = f.read()
    assert ixyz is not None, 'needs to have ixyz'
    assert np.allclose(angles, 90.), 'trajs are not in orthogonal boxes'
    assert np.isclose(types, types[0]).all(), 'atom type changes'
    assert np.isclose(masses, masses[0]).all(), 'atom mass changes'
    if not np.isclose(lattices, lattices[0]).all():
        warnings.warn('box size changes.')
    ixyz = ixyz.astype('int32')
    types = types[0]
    masses = masses[0]
    assert types.shape[0] == coords.shape[1] == ixyz.shape[1] == masses.shape[0]
    assert coords.shape[0] == lattices.shape[0] == ixyz.shape[0]

    unwrapped_coords = coords + lattices[:, np.newaxis] * ixyz

    # Wrap coords
    imag_loc = np.floor_divide(coords, lattices[:, np.newaxis])
    wrapped_coords = coords - imag_loc * lattices[:, np.newaxis]
    assert np.alltrue((np.max(wrapped_coords, axis=1) -
                       np.min(wrapped_coords, axis=1)) <= lattices)

    if use_mass:
        new_types = []
        for mass in masses:
            diffs = np.abs(mass - _ATOM_MASSES)
            atomic_number = np.argmin(diffs)
            assert diffs[atomic_number] <= tol, 'diff is {} for {}'.format(
                diffs[atomic_number], mass)
            new_types.append(atomic_number)
        new_types = np.array(new_types, dtype=np.int32)
        return wrapped_coords, lattices, new_types, unwrapped_coords
    else:
        return wrapped_coords, lattices, types, unwrapped_coords


def distance_pbc(x0, x1, dimensions):
    """
    Distance between atoms in periodic boundary conditions.

    x0: (3, ) numpy array
    x1: (m, 3) numpy array
    dimensions: (3, ) numpy array
    """
    assert x0.shape[0] == 3 and x1.shape[1] == 3
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


def batch_pdist_pbc(x0, x1, dimensions):
    """
    Calculate pairwise distance between atoms in periodic boundary conditions
    in batch.

    x0: shape (B, N0, 3)
    x1: shape (B, N1, 3)
    dimensions: shape (B, 3)

    returns
    pdist: shape (B, N1, N0)
    """
    assert x0.shape[0] == x1.shape[0]
    assert x0.shape[-1] == x1.shape[-1] == 3
    distances = np.abs(x0[:, np.newaxis] - x1[:, :, np.newaxis])
    dimensions = dimensions[:, np.newaxis, np.newaxis]
    distances = np.where(distances > 0.5 * dimensions, distances - dimensions,
                         distances)
    return np.sqrt((distances ** 2).sum(axis=-1))


def relative_coord_pbc(x0, x1, dimensions):
    """
    The coordinate of atoms in x1 relative to x0 in periodic boundary
    conditions.

    x0: (3, ) numpy array
    x1: (m, 3) numpy array
    dimensions: (3, ) numpy array
    """
    assert x0.shape[0] == 3 and x1.shape[1] == 3
    delta = x1 - x0
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    delta = np.where(delta < -0.5 * dimensions, delta + dimensions, delta)
    return delta


class ExtendedLAMMPSTrajectoryFile(LAMMPSTrajectoryFile):
    def __init__(self, filename, mode='r', force_overwrite=True):
        super(ExtendedLAMMPSTrajectoryFile, self).__init__(
            filename, mode, force_overwrite)

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read data from a lammpstrj file.

        Parameters
        ----------
        n_frames : int, None
            The number of frames you would like to read from the file.
            If None, all of the remaining frames will be loaded.
        stride : np.ndarray, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates
            from the file.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
        cell_lengths : np.ndarray, None
            The lengths (a,b,c) of the unit cell for each frame, or None if
            the information is not present in the file.
        cell_angles : np.ndarray, None
            The angles (\alpha, \beta, \gamma) defining the unit cell for
            each frame, or None if  the information is not present in the file.
        """
        if not self._mode == 'r':
            raise ValueError('read() is only available when file is opened '
                             'in mode="r"')

        if n_frames is None:
            frame_counter = itertools.count()
        else:
            frame_counter = xrange(n_frames)

        if stride is None:
            stride = 1

        all_coords, all_lengths, all_angles, all_types, all_ixyz, all_mass =\
            [], [], [], [], [], []
        for _ in frame_counter:
            try:
                frame_coords, frame_lengths, frame_angles, frame_types,\
                    frame_ixyz, frame_mass = self._read()
                if atom_indices is not None:
                    raise NotImplementedError
            except _EOF:
                break

            all_coords.append(frame_coords)
            all_lengths.append(frame_lengths)
            all_angles.append(frame_angles)
            all_types.append(frame_types)
            all_ixyz.append(frame_ixyz)
            all_mass.append(frame_mass)

            for j in range(stride - 1):
                # throw away these frames
                try:
                    self._read()
                except _EOF:
                    break

        all_coords = np.array(all_coords)
        all_lengths = np.array(all_lengths, dtype=np.float32)
        all_angles = np.array(all_angles, dtype=np.float32)
        all_types = np.array(all_types, dtype=np.int32)
        if all_ixyz[0] is None:
            all_ixyz = None
        else:
            all_ixyz = np.array(all_ixyz, dtype=np.int32)
        if all_mass[0] is None:
            all_mass = None
        else:
            all_mass = np.array(all_mass, dtype=np.float32)
        return (all_coords, all_lengths, all_angles, all_types, all_ixyz,
                all_mass)

    def _read(self):
        """Read a single frame. """

        # --- begin header ---
        first = self._fh.readline()  # ITEM: TIMESTEP
        if first == '':
            raise _EOF()
        self._fh.readline()  # timestep
        self._fh.readline()  # ITEM: NUMBER OF ATOMS
        self._n_atoms = int(self._fh.readline())  # num atoms

        box_header = self._fh.readline().split()  # ITEM: BOX BOUNDS
        self._line_counter += 5
        if len(box_header) == 9:
            lengths, angles = self.parse_box('triclinic')
        elif len(box_header) == 6:
            lengths, angles = self.parse_box('orthogonal')
        else:
            raise IOError('lammpstrj parse error on line {0:d} of "{1:s}". '
                          'This file does not appear to be a valid '
                          'lammpstrj file.'.format(self._line_counter,
                                                   self._filename))

        column_headers = self._fh.readline().split()[2:]  # ITEM: ATOMS ...
        if self._frame_index == 0:
            # Detect which columns the atom index, type and coordinates are.
            columns = {header: idx for idx, header
                       in enumerate(column_headers)}

            # Make sure the file contains an x, y, and z-coordinate of the same
            # style.
            coord_keywords = [('x', 'y', 'z'),  # unscaled
                              ('xs', 'ys', 'zs'),  # scaled
                              ('xu', 'yu', 'zu'),  # unwrapped
                              ('xsu', 'ysu', 'zsu')]  # scaled and unwrapped
            for keywords in coord_keywords:
                if set(keywords).issubset(column_headers):
                    break
            else:
                raise IOError('Invalid .lammpstrj file. Must contain x, y, '
                              'and z coordinates that all adhere to the same '
                              'style.')

            try:
                self._atom_index_column = columns['id']
                if 'element' in columns:
                    self._atom_type_column = columns['element']
                else:
                    self._atom_type_column = columns['type']
                self._xyz_columns = [columns[keywords[0]],
                                     columns[keywords[1]],
                                     columns[keywords[2]]]
                if 'ix' in columns:
                    self._ixyz_columns = [columns['ix'],
                                          columns['iy'],
                                          columns['iz']]
                else:
                    self._ixyz_columns = None
                if 'mass' in columns:
                    self._mass_columns = columns['mass']
                else:
                    self._mass_columns = None
            except KeyError:
                raise IOError("Invalid .lammpstrj file. Must contain 'id', "
                              "'type', 'x*', 'y*' and 'z*' entries.")
        self._line_counter += 4
        # --- end header ---

        xyz = np.empty(shape=(self._n_atoms, 3))
        types = np.empty(shape=self._n_atoms, dtype='int')
        if self._ixyz_columns is not None:
            ixyz = np.empty(shape=(self._n_atoms, 3), dtype='int')
        else:
            ixyz = None
        if self._mass_columns is not None:
            mass = np.empty(shape=self._n_atoms, dtype='float')
        else:
            mass = None

        # --- begin body ---
        prev_atom_index = 0
        for idx in xrange(self._n_atoms):
            line = self._fh.readline()
            if line == '':
                raise _EOF()
            split_line = line.split()
            try:
                atom_index = int(split_line[self._atom_index_column])
                assert atom_index > prev_atom_index, 'atom_index is not sorted'
                prev_atom_index = atom_index
                elem = split_line[self._atom_type_column]
                types[idx] = int(elem)
                xyz[idx] = [float(split_line[column])
                            for column in self._xyz_columns]
                if self._ixyz_columns is not None:
                    ixyz[idx] = [int(split_line[column])
                                 for column in self._ixyz_columns]
                if self._mass_columns is not None:
                    mass[idx] = float(split_line[self._mass_columns])
            except Exception:
                raise IOError('lammpstrj parse error on line {0:d} of "{1:s}".'
                              ' This file does not appear to be a valid '
                              'lammpstrj file.'.format(self._line_counter,
                                                       self._filename))
            self._line_counter += 1
        # --- end body ---

        self._frame_index += 1
        return xyz, lengths, angles, types, ixyz, mass
