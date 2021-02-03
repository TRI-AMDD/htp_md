"""
Module for definition of Trajectory and specific subtypes.
"""
import numpy as np
from mdtraj.formats import LAMMPSTrajectoryFile
from mdtraj.formats.lammpstrj import _EOF
from mdtraj.utils.six.moves import xrange
import itertools

from htpmd.constants import ATOM_MASSES


class Trajectory:

    @property
    def all_coords(self):
        return self._all_coords

    @all_coords.setter
    def all_coords(self, val):
        self._all_coords = val

    @property
    def all_lengths(self):
        return self._all_lengths

    @all_lengths.setter
    def all_lengths(self, val):
        self._all_lengths = val

    @property
    def all_angles(self):
        return self._all_angles

    @all_angles.setter
    def all_angles(self, val):
        self._all_angles = val

    @property
    def all_types(self):
        return self._all_types

    @all_types.setter
    def all_types(self, val):
        self._all_types = val

    @property
    def all_ixyz(self):
        return self._all_ixyz

    @all_ixyz.setter
    def all_ixyz(self, val):
        self._all_ixyz = val

    @property
    def all_masses(self):
        return self._all_masses

    @all_masses.setter
    def all_masses(self, val):
        self._all_masses = val

    @property
    def unwrapped_coords(self):
        return self._unwrapped_coords

    @unwrapped_coords.setter
    def unwrapped_coords(self, val):
        self._unwrapped_coords = val

    @property
    def wrapped_coords(self):
        return self._wrapped_coords

    @wrapped_coords.setter
    def wrapped_coords(self, val):
        self._wrapped_coords = val

    @property
    def raw_types(self):
        return self._raw_types

    @raw_types.setter
    def raw_types(self, val):
        self._raw_types = val

    @property
    def atom_types(self):
        return self._atom_types

    @atom_types.setter
    def atom_types(self, val):
        self._atom_types = val

    @property
    def lattices(self):
        return self.all_lengths


class ExtendedLAMMPSTrajectoryFile(LAMMPSTrajectoryFile, Trajectory):

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

        self.all_coords = np.array(all_coords)
        self.all_lengths = np.array(all_lengths, dtype=np.float32)
        self.all_angles = np.array(all_angles, dtype=np.float32)
        self.all_types = np.array(all_types, dtype=np.int32)
        if all_ixyz[0] is None:
            self.all_ixyz = None
        else:
            self.all_ixyz = np.array(all_ixyz, dtype=np.int32)
        if all_mass[0] is None:
            self.all_masses = None
        else:
            self.all_masses = np.array(all_mass, dtype=np.float32)
        return self.all_coords, self.all_lengths, self.all_angles, self.all_types, self.all_ixyz, self.all_masses

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

    def determine_lattice_and_coords(self):

        self.unwrapped_coords = self.all_coords + self.all_lengths[:, np.newaxis] * self.all_ixyz

        imag_loc = np.floor_divide(self.all_coords, self.all_lengths[:, np.newaxis])
        self.wrapped_coords = self.all_coords - imag_loc * self.all_lengths[:, np.newaxis]

        self.raw_types = self.all_types[0]

        new_types = []
        for mass in self.all_masses[0]:
            diffs = np.abs(mass - ATOM_MASSES)
            atomic_number = np.argmin(diffs)
            new_types.append(atomic_number)
        self.atom_types = np.array(new_types, dtype=np.int32)

    def remove_drift(self):
        self.unwrapped_coords = _get_coords_without_drift(
            self.unwrapped_coords, self.atom_types)


def _compute_center_of_mass(coords, atom_types):
    element_masses = np.array(ATOM_MASSES)
    atom_masses = element_masses[atom_types]
    return (np.sum(coords * atom_masses[np.newaxis, :, np.newaxis], axis=1) /
            np.sum(atom_masses))


def _get_coords_without_drift(coords, atom_types):
    cos_coord = _compute_center_of_mass(coords, atom_types)
    return coords - cos_coord[:, np.newaxis]