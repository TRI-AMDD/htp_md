from __future__ import print_function, division
import numpy as np

from htpmd.trajectory.base import ExtendedLAMMPSTrajectoryFile
from htpmd.constants import ATOM_MASSES

_TOL = 1e-2


def load_lammps(lammps_file, tol=_TOL):
    """Load a lammpstraj file.

    Args:
        lammps_file: file path

    Returns:
        coords: shape (F, N, 3)
        lattices: shape (F, 3)
        types: shape (N,), the raw type of atoms in lammps
        atom_types: shape (N,), the atom types by atomic number
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

    new_types = []
    for mass in masses:
        diffs = np.abs(mass - ATOM_MASSES)
        atomic_number = np.argmin(diffs)
        assert diffs[atomic_number] <= tol, 'diff is {} for {}'.format(
            diffs[atomic_number], mass)
        new_types.append(atomic_number)
    new_types = np.array(new_types, dtype=np.int32)
    return wrapped_coords, lattices, types, new_types, unwrapped_coords


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
