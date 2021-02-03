import unittest
import json
import os
import pytest
import numpy as np
from htpmd.shared.polymer import (
    compute_diffusivity, compute_polymer_diffusivity, compute_conductivity,
    compute_molarity)
from htpmd.trajectory.load import LammpsTrajectoryLoader


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
    assert approx_equal(li_diff, compute_diffusivity(trajectory, target_type=90))
    assert approx_equal(tfsi_diff, compute_diffusivity(trajectory, target_type=93))


@pytest.mark.parametrize(
    'dir_name,diff',
    [
        ('test_data/9-0-246295613-0', 1.4832525239039016e-06),
        ('test_data/9-0-413610210-0', 9.208230905677645e-07),
    ])
def test_compute_polymer_diffusivity(dir_name, diff):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(diff, compute_polymer_diffusivity(trajectory))


@pytest.mark.parametrize(
    'dir_name,cond',
    [
        ('test_data/9-0-246295613-0', 0.004507580946462734),
        ('test_data/9-0-413610210-0', 0.005354250978187845),
    ])
def test_compute_conductivity(dir_name, cond):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    pop_mat = get_population_matrix(dir_name)
    assert approx_equal(cond, compute_conductivity(trajectory, pop_mat=pop_mat))


@pytest.mark.parametrize(
    'dir_name,mod',
    [
        ('test_data/9-0-246295613-0', 1.4631208109388065),
        ('test_data/9-0-413610210-0', 1.3711409170333055),
    ])
def test_compute_conductivity(dir_name, mod):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(mod, compute_molarity(trajectory))
