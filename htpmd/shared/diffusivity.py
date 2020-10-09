"""
Module for methods to compute diffusivity.
"""
import numpy as np

from htpmd.shared.utils import check_params


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
    required_parameters = ('target_type', )
    check_params(required_parameters, params)

    target_idx = np.nonzero(trajectory.raw_types == params['target_type'])[0]
    target_coords = trajectory.unwrapped_coords[:, target_idx]
    msd = np.mean(np.sum((target_coords[-1] - target_coords[0])**2, axis=-1))
    diffusivity = msd / (len(target_coords) - 1) / 6 * 5e-5  # cm^2/s
    return diffusivity
