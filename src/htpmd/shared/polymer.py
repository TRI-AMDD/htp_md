"""
Module for methods to compute properties for the polymer electrolytes.

This module implements several analysis functions for the polymer electrolyte
MD data described in paper [1]. Currently, they cannot be directly used in
other datasets because some data-specific constants are written in the
functions. However, users can still reuse these functions after minor
modifications. We are working on making the functions more general and users
are encouraged to contribute.

[1] Xie, et al. arXiv preprint arXiv:2101.05339 (2021).
"""
import numpy as np

from htpmd.shared.utils import check_params
from htpmd.constants import ATOM_MASSES, \
    FARADAY_CONSTANT, \
    BOLTZMANN_CONSTANT, \
    ANGSTROM, \
    CENTIMETER, \
    NANOSECOND, \
    PICOSECOND, \
    KILOGRAM, \
    RawType, \
    AtomType
from pymatgen.core.structure import Structure


def compute_diffusivity(trajectory, **params):
    """
    Description:
        Diffusivity of a specified atom type (unit: cm^2/s).

        Example:
        `li_diffusivity = compute_diffusivity(trajectory, **{'target_type': 90})`
        or
        `li_diffusivity = compute_diffusivity(trajectory, target_type=90)`

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            target_type (int)

    Returns:
        float:                                          diffusivity (for given target_type)

    """
    required_parameters = ('target_type', 'time_step')
    check_params(required_parameters, params)
    delta_t = params['time_step'] * PICOSECOND

    target_idx = np.nonzero(trajectory.raw_types == params['target_type'])[0]
    target_coords = trajectory.unwrapped_coords[:, target_idx]
    msd = np.mean(np.sum((target_coords[-1] - target_coords[0])**2, axis=-1))
    diffusivity = msd / (len(target_coords) - 1) / 6 / delta_t  # A^2/s
    diffusivity = diffusivity * (ANGSTROM / CENTIMETER)**2  # cm^2/s
    return diffusivity


def compute_diffusivity_array(trajectory, **params):
    """
    Description:
        Diffusivity of a specified atom type (unit: cm^2/s).

        Example:
        `li_diffusivity = compute_diffusivity(trajectory, **{'target_type': 90})`
        or
        `li_diffusivity = compute_diffusivity(trajectory, target_type=90)`

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            target_type (int)

    Returns:
        float:                                          diffusivity (for given target_type)

    """
    required_parameters = ('target_type', 'time_step')
    check_params(required_parameters, params)
    delta_t = params['time_step'] * PICOSECOND

    target_idx = np.nonzero(trajectory.raw_types == params['target_type'])[0]
    target_coords = trajectory.unwrapped_coords[:, target_idx]

    msd = [np.mean(np.sum((target_coords[t] - target_coords[0]) ** 2, axis=-1)) for t in range(1, len(target_coords))]
    idx_len = np.arange(1, len(msd) + 1, 1)

    diffusivity = msd / idx_len / 6 / delta_t  # A^2/s
    diffusivity = diffusivity * (ANGSTROM / CENTIMETER) ** 2  # cm^2/s
    return diffusivity


def compute_polymer_diffusivity(trajectory, **params):
    """
    Description:
        Diffusivity of the polymer, defined by the average diffusivity of
        N, O, S atoms in the polymer chain.

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:

    Returns:
        float:                                          diffusivity (for polymer)

    """
    required_parameters = ('time_step', 'polymer_raw_type_range', 'polymer_solvate_types')
    check_params(required_parameters, params)
    delta_t = params['time_step'] * PICOSECOND

    solvate_types_list = [trajectory.atom_types == atom_type for atom_type in params['polymer_solvate_types']]
    solvate_types = np.logical_or.reduce(solvate_types_list)
    poly_types = np.logical_and(trajectory.raw_types >= params['polymer_raw_type_range'][0],
                                trajectory.raw_types <= params['polymer_raw_type_range'][1])
    poly_solvate_types = poly_types & solvate_types
    poly_solvate_idx = np.nonzero(poly_solvate_types)[0]
    target_coords = trajectory.unwrapped_coords[:, poly_solvate_idx]
    msd = np.mean(np.sum((target_coords[-1] - target_coords[0])**2, axis=-1))
    diffusivity = msd / (len(target_coords) - 1) / 6 / delta_t  # A^2/s
    diffusivity = diffusivity * (ANGSTROM / CENTIMETER)**2  # cm^2/s
    return diffusivity


def compute_polymer_diffusivity_array(trajectory, **params):
    """
    Description:
        Diffusivity of the polymer, defined by the average diffusivity of
        N, O, S atoms in the polymer chain.

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:

    Returns:
        float:                                          diffusivity (for polymer)

    """
    required_parameters = ('time_step', 'polymer_raw_type_range', 'polymer_solvate_types')
    check_params(required_parameters, params)
    delta_t = params['time_step'] * PICOSECOND

    solvate_types_list = [trajectory.atom_types == atom_type for atom_type in params['polymer_solvate_types']]
    solvate_types = np.logical_or.reduce(solvate_types_list)
    poly_types = np.logical_and(trajectory.raw_types >= params['polymer_raw_type_range'][0],
                                trajectory.raw_types <= params['polymer_raw_type_range'][1])
    poly_solvate_types = poly_types & solvate_types
    poly_solvate_idx = np.nonzero(poly_solvate_types)[0]
    target_coords = trajectory.unwrapped_coords[:, poly_solvate_idx]
    msd = [np.mean(np.sum((target_coords[t] - target_coords[0]) ** 2, axis=-1)) for t in range(1, len(target_coords))]
    idx_len = np.arange(1, len(msd) + 1, 1)

    diffusivity = msd / idx_len / 6 / delta_t  # A^2/s
    diffusivity = diffusivity * (ANGSTROM / CENTIMETER) ** 2  # cm^2/s
    return diffusivity


def compute_molality(trajectory, **params):
    """
    Description:
        Molality of the polymer/salt mixture (unit: mol Li / kg polymer).

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:

    Returns:
        float:                                          molality

    """
    required_parameters = ('polymer_raw_type_range', 'cation_raw_type')
    check_params(required_parameters, params)

    poly_types = np.logical_and(trajectory.raw_types >= params['polymer_raw_type_range'][0],
                                trajectory.raw_types <= params['polymer_raw_type_range'][1])
    poly_idx = np.nonzero(poly_types)[0]
    li_idx = np.nonzero(trajectory.raw_types == params['cation_raw_type'])[0]

    atom_masses = np.array(ATOM_MASSES)

    poly_mass = np.sum(atom_masses[trajectory.atom_types[poly_idx]])

    return float(len(li_idx)) / poly_mass * KILOGRAM


def compute_conductivity(trajectory, **params):
    """
    Description:
        Compute the conductivity and transference number using the
        cluster-Nernst-Einstein described in the following paper.

        France-Lanord and Grossman. "Correlations from ion pairing and the
        Nernst-Einstein equation." Physical review letters 122.13 (2019): 136001.

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            pop_mat: np.array

    Returns:
        float:                                          conductivity (unit: S/cm)
        float:                                          transference_number

    """
    required_parameters = ('pop_mat', 'time_step', 'temperature', 'cation_raw_type', 'anion_raw_type')
    check_params(required_parameters, params)

    max_cluster = 10
    T = params['temperature']
    pop_mat = params['pop_mat']
    z_i, z_j = 1, 1  # charges carried by cation and anions

    li_diff = compute_diffusivity(trajectory, target_type=params['cation_raw_type'], **params)  # cm^2/s
    tfsi_diff = compute_diffusivity(trajectory, target_type=params['anion_raw_type'], **params)  # cm^2/s

    assert np.isclose(trajectory.lattices[0:1], trajectory.lattices).all()

    V = np.prod(trajectory.lattices[0]) * (ANGSTROM / CENTIMETER)**3  # cm^3

    cond = 0.
    total_ion = 0.
    tn_numerator, tn_denominator = 0., 0.

    for i in range(max_cluster):
        for j in range(max_cluster):
            if i > j:
                cond += FARADAY_CONSTANT**2 / V / BOLTZMANN_CONSTANT / T * \
                    (i * z_i - j * z_j)**2 * pop_mat[i, j] * li_diff
                tn_numerator += i * z_i * (i * z_i - j * z_j) * pop_mat[i, j] * li_diff
                tn_denominator += (i * z_i - j * z_j)**2 * pop_mat[i, j] * li_diff
            elif i < j:
                cond += FARADAY_CONSTANT**2 / V / BOLTZMANN_CONSTANT / T * \
                    (i * z_i - j * z_j)**2 * pop_mat[i, j] * tfsi_diff
                tn_numerator += i * z_i * (i * z_i - j * z_j) * pop_mat[i, j] * tfsi_diff
                tn_denominator += (i * z_i - j * z_j)**2 * pop_mat[i, j] * tfsi_diff
            else:
                pass
            total_ion += (i + j) * pop_mat[i, j]
    tn = tn_numerator / tn_denominator

    return cond, tn


def compute_conductivity_array(trajectory, **params):
    """
    Description:
        Compute the conductivity and transference number using the
        cluster-Nernst-Einstein described in the following paper.

        France-Lanord and Grossman. "Correlations from ion pairing and the
        Nernst-Einstein equation." Physical review letters 122.13 (2019): 136001.

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            pop_mat: np.array

    Returns:
        float:                                          conductivity (unit: S/cm)
        float:                                          transference_number

    """
    required_parameters = ('pop_mat', 'time_step', 'temperature', 'cation_raw_type', 'anion_raw_type')
    check_params(required_parameters, params)

    max_cluster = 10
    T = params['temperature']
    pop_mat = params['pop_mat']
    z_i, z_j = 1, 1  # charges carried by cation and anions

    li_diff = compute_diffusivity_array(trajectory, target_type=params['cation_raw_type'], **params)  # cm^2/s
    tfsi_diff = compute_diffusivity_array(trajectory, target_type=params['anion_raw_type'], **params)  # cm^2/s

    assert np.isclose(trajectory.lattices[0:1], trajectory.lattices).all()

    V = np.prod(trajectory.lattices[0]) * (ANGSTROM / CENTIMETER)**3  # cm^3

    cond = 0.
    total_ion = 0.
    tn_numerator, tn_denominator = 0., 0.

    for i in range(max_cluster):
        for j in range(max_cluster):
            if i > j:
                cond += FARADAY_CONSTANT**2 / V / BOLTZMANN_CONSTANT / T * \
                    (i * z_i - j * z_j)**2 * pop_mat[i, j] * li_diff
                tn_numerator += i * z_i * (i * z_i - j * z_j) * pop_mat[i, j] * li_diff
                tn_denominator += (i * z_i - j * z_j)**2 * pop_mat[i, j] * li_diff
            elif i < j:
                cond += FARADAY_CONSTANT**2 / V / BOLTZMANN_CONSTANT / T * \
                    (i * z_i - j * z_j)**2 * pop_mat[i, j] * tfsi_diff
                tn_numerator += i * z_i * (i * z_i - j * z_j) * pop_mat[i, j] * tfsi_diff
                tn_denominator += (i * z_i - j * z_j)**2 * pop_mat[i, j] * tfsi_diff
            else:
                pass
            total_ion += (i + j) * pop_mat[i, j]
    tn = tn_numerator / tn_denominator

    return cond, tn


def compute_msd_curve(trajectory, **params):
    """
    Description:
        Computes Mean Squared Dispacement curve (unit: ns, A^2)

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            target_type (int)

    Returns:
        ts:                                              time series (np.array)
        msds:                                            mean squared displacement (np.array)
    """

    required_parameters = ('target_type', 'time_step')
    check_params(required_parameters, params)
    delta_t = params['time_step'] * PICOSECOND

    target_idx = np.nonzero(trajectory.raw_types == params['target_type'])[0]
    target_coords = trajectory.unwrapped_coords[:, target_idx]

    ts = np.linspace(1, target_coords.shape[0] - 1, 100, dtype=int)
    msds = np.array([
        np.mean(np.sum((target_coords[t:] - target_coords[:-t])**2, axis=-1))
        for t in ts])
    # Convert to ns
    ts = ts * delta_t * NANOSECOND / PICOSECOND
    return ts, msds


def compute_non_avg_msd_curve(trajectory, **params):
    """
    Description:
        Computes Mean Squared Displacement curve (unit: ns, A^2) without time averaging

    Version: 1.0.0

    Author:
        Name:                                           Arash Khajeh
        Affiliation:                                    TRI
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            target_type (int)

    Returns:
        ts:                                              time series (np.array)
        msds:                                            mean squared displacement (np.array)
    """

    required_parameters = ('target_type', 'time_step')
    check_params(required_parameters, params)
    delta_t = params['time_step'] * PICOSECOND

    target_idx = np.nonzero(trajectory.raw_types == params['target_type'])[0]
    target_coords = trajectory.unwrapped_coords[:, target_idx]

    ts = np.linspace(1, target_coords.shape[0] - 1, target_coords.shape[0] - 1, dtype=int)
    msds = np.array([
        np.mean(np.sum((target_coords[t] - target_coords[0])**2, axis=-1))
        for t in ts])
    # Convert to ns
    ts = ts * delta_t * NANOSECOND / PICOSECOND
    return ts, msds


def compute_ngp_curve(trajectory, **params):
    """
    Description:
        Computes non-Gaussian parameter curve (unit: ns, dimensionless)
    Version: 1.0.0
    Author:
        Name:                                           Arthur France-Lanord
        Affiliation:                                    MIT
        Email:                                          <optional>
    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            target_type (int)
    Returns:
        ts:                                              time series (np.array)
        ngps:                                            non-Gaussian parameter (np.array)
    """

    required_parameters = ('target_type', 'time_step')
    check_params(required_parameters, params)
    delta_t = params['time_step'] * PICOSECOND

    target_idx = np.nonzero(trajectory.raw_types == params['target_type'])[0]
    target_coords = trajectory.unwrapped_coords[:, target_idx]

    ts = np.linspace(1, target_coords.shape[0] - 1, 100, dtype=int)
    msds = np.array([
        np.mean(np.sum((target_coords[t:] - target_coords[:-t])**2, axis=-1))
        for t in ts])
    mfds = np.array([
        np.mean(np.sum((target_coords[t:] - target_coords[:-t])**2, axis=-1)**2)
        for t in ts])
    ngps = np.array([
        (3*mfds)/(5*msds*msds)-1.])

    # Convert to ns
    ts = ts * delta_t * NANOSECOND / PICOSECOND
    return ts, ngps


def get_cif_at_frame(trajectory, **params):
    """
    Description:
        Computes the cif text representation of the structure at frame k.

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            k (int)

    Returns:

    """
    required_parameters = ('k', )
    check_params(required_parameters, params)

    structure = Structure(
        lattice=np.diag(trajectory.lattices[params['k']]),
        species=trajectory.atom_types,
        coords=trajectory.wrapped_coords[params['k']],
        coords_are_cartesian=True,
        to_unit_cell=True
    )
    return structure.to('cif')


def compute_displacement(trajectory, **params):
    """
    Description:
        Compute the mean or max displacement a specified atom type (unit: A).

        Example:
        `li_mean_disp = compute_displacement(trajectory, **{'target_type': 90, 'type': 'mean'})`

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            target_type (int)
                                                            type: 'mean' or 'max'

    Returns:
        float:                                          diffusivity (for given target_type)

    """
    required_parameters = ('target_type', 'type')
    check_params(required_parameters, params)

    target_idx = np.nonzero(trajectory.raw_types == params['target_type'])[0]
    target_coords = trajectory.unwrapped_coords[:, target_idx]
    disp = np.sqrt(np.sum((target_coords[-1] - target_coords[0])**2, axis=-1))
    if params['type'] == 'mean':
        disp = np.mean(disp)
    elif params['type'] == 'max':
        disp = np.max(disp)
    else:
        raise NotImplementedError
    return disp


def compute_simulation_length(trajectory, **params):
    """
    Description:
        Compute the total length of simulation in ns

        Example:
        `simulation_length = compute_displacement(trajectory, **{'time_step': 2.0})`

    Version: 1.0.0

    Author:
        Name:                                           Tian Xie
        Affiliation:                                    MIT
        Email:                                          <optional>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:

    Returns:
        float:                                          total simulation length in ns

    """
    required_parameters = ('time_step',)
    check_params(required_parameters, params)
    delta_t = params['time_step'] * PICOSECOND

    total_t = (trajectory.unwrapped_coords.shape[0] - 1) * delta_t

    return total_t / NANOSECOND
