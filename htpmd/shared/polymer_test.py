import unittest
import json
import os
import pytest
import numpy as np
from htpmd.shared.polymer import (
    compute_diffusivity, compute_polymer_diffusivity, compute_conductivity,
    compute_molarity)
from htpmd.trajectory.load import (LammpsTrajectoryLoader, get_population_matrix)


def approx_equal(val1, val2):
    return pytest.approx(val1, rel=1e-3) == val2


@pytest.mark.parametrize(
    'dir_name,li_diff,tfsi_diff',
    [
        ('test_data/9-0-246295613-0', 1.2874560352233332e-06, 2.5612469517951516e-06),
        ('test_data/9-0-413610210-0', 1.035696582236192e-06, 2.3270058381337847e-06),
    ])
def test_compute_diffusivity(dir_name, li_diff, tfsi_diff):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(li_diff, compute_diffusivity(trajectory, target_type=90, time_step=2.))
    assert approx_equal(tfsi_diff, compute_diffusivity(trajectory, target_type=93, time_step=2.))


@pytest.mark.parametrize(
    'dir_name,diff',
    [
        ('test_data/9-0-246295613-0', 1.4832525239039016e-06),
        ('test_data/9-0-413610210-0', 9.208230905677645e-07),
    ])
def test_compute_polymer_diffusivity(dir_name, diff):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(diff, compute_polymer_diffusivity(
        trajectory, time_step=2., polymer_raw_type_range=[0, 89], polymer_solvate_types=[7, 8, 16]))


@pytest.mark.parametrize(
    'dir_name,cond,tn',
    [
        ('test_data/9-0-246295613-0', 0.004507580946462734, 0.18290430863167545),
        ('test_data/9-0-413610210-0', 0.005354250978187845, -0.061131167509913646),
    ])
def test_compute_conductivity(dir_name, cond, tn):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    pop_mat = get_population_matrix(dir_name)
    cond_result, tn_result = compute_conductivity(
        trajectory, pop_mat=pop_mat, time_step=2., temperature=353.0, cation_raw_type=90, anion_raw_type=93)
    assert approx_equal(cond, cond_result)
    assert approx_equal(tn, tn_result)


@pytest.mark.parametrize(
    'dir_name,mod',
    [
        ('test_data/9-0-246295613-0', 1.4631208109388065),
        ('test_data/9-0-413610210-0', 1.3711409170333055),
    ])
def test_compute_conductivity_2(dir_name, mod):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(mod, compute_molarity(
        trajectory, polymer_raw_type_range=[0, 89], cation_raw_type=90))
