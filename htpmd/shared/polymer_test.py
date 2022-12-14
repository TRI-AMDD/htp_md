from curses import meta
import pytest
from htpmd.shared.polymer import (
    compute_diffusivity, compute_polymer_diffusivity, compute_conductivity,
    compute_molality, compute_displacement, compute_simulation_length, compute_diffusivity_array,
    compute_conductivity_array, compute_polymer_diffusivity_array, compute_density,
    compute_degree_polymerization)
from htpmd.trajectory.load import (LammpsTrajectoryLoader, get_population_matrix, get_metadata)


def approx_equal(val1, val2):
    return pytest.approx(val1, rel=1.0e-4) == val2


@pytest.mark.parametrize(
    'dir_name,li_diff,tfsi_diff',
    [
        ('test_data/9-0-246295613-0', 1.2874560352233334e-06, 2.561246951795152e-06),
        ('test_data/9-0-413610210-0', 1.0356965822361922e-06, 2.327005838133785e-06),
    ])
def test_compute_diffusivity(dir_name, li_diff, tfsi_diff):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(li_diff, compute_diffusivity(trajectory, target_type=90, time_step=2.))
    assert approx_equal(tfsi_diff, compute_diffusivity(trajectory, target_type=93, time_step=2.))


@pytest.mark.parametrize(
    'dir_name,li_diff,tfsi_diff',
    [
        ('test_data/9-0-246295613-0', 1.2874560352233334e-06, 2.561246951795152e-06),
        ('test_data/9-0-413610210-0', 1.0356965822361922e-06, 2.327005838133785e-06),
    ])
def test_compute_diffusivity_array(dir_name, li_diff, tfsi_diff):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(li_diff, compute_diffusivity_array(trajectory, target_type=90, time_step=2.)[-1])
    assert approx_equal(tfsi_diff, compute_diffusivity_array(trajectory, target_type=93, time_step=2.)[-1])


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
    'dir_name,diff',
    [
        ('test_data/9-0-246295613-0', 1.4832525239039016e-06),
        ('test_data/9-0-413610210-0', 9.208230905677645e-07),
    ])
def test_compute_polymer_diffusivity_array(dir_name, diff):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(diff, compute_polymer_diffusivity_array(
        trajectory, time_step=2., polymer_raw_type_range=[0, 89], polymer_solvate_types=[7, 8, 16])[-1])


@pytest.mark.parametrize(
    'dir_name,cond,tn',
    [
        ('test_data/9-0-246295613-0', 0.004798255867131224, 0.12111878931985638),
        ('test_data/9-0-413610210-0', 0.0055822299628502615, -0.13037903387492417),
    ])
def test_compute_conductivity(dir_name, cond, tn):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    stacked_population, pop_mat = get_population_matrix(dir_name)
    cond_result, tn_result = compute_conductivity(
        trajectory, pop_mat=pop_mat, time_step=2., temperature=353.0, cation_raw_type=90, anion_raw_type=93)
    assert approx_equal(cond, cond_result)
    assert approx_equal(tn, tn_result)


@pytest.mark.parametrize(
    'dir_name,cond,tn',
    [
        ('test_data/9-0-246295613-0', 0.004798255867131224, 0.12111878931985638),
        ('test_data/9-0-413610210-0', 0.0055822299628502615, -0.13037903387492417),
    ])
def test_compute_conductivity_array(dir_name, cond, tn):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    stacked_population, pop_mat = get_population_matrix(dir_name)
    cond_result, tn_result = compute_conductivity_array(
        trajectory, pop_mat=pop_mat, time_step=2., temperature=353.0, cation_raw_type=90, anion_raw_type=93)
    assert approx_equal(cond, cond_result[-1])
    assert approx_equal(tn, tn_result[-1])


@pytest.mark.parametrize(
    'dir_name,mod',
    [
        ('test_data/9-0-246295613-0', 1.4631208109388065),
        ('test_data/9-0-413610210-0', 1.3711409170333055),
    ])
def test_compute_molality(dir_name, mod):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(mod, compute_molality(
        trajectory, polymer_raw_type_range=[0, 89], cation_raw_type=90))


@pytest.mark.parametrize(
    'dir_name,type',
    [
        ('test_data/9-0-246295613-0', 'mean'),
        ('test_data/9-0-413610210-0', 'mean'),
        ('test_data/9-0-246295613-0', 'max'),
        ('test_data/9-0-413610210-0', 'max'),
    ])
def test_compute_displacement(dir_name, type):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert compute_displacement(trajectory, target_type=90, type=type) > 0.
    assert compute_displacement(trajectory, target_type=93, type=type) > 0.


@pytest.mark.parametrize(
    'dir_name,total_length',
    [
        ('test_data/9-0-246295613-0', 0.012),
        ('test_data/9-0-413610210-0', 0.014),
    ])
def test_compute_simulation_length(dir_name, total_length):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(compute_simulation_length(trajectory, time_step=2.0), total_length)
    assert approx_equal(compute_simulation_length(trajectory, time_step=2.0), total_length)


@pytest.mark.parametrize(
    'dir_name,true_density',
    [
        ('test_data/9-0-246295613-0', 1.1921),
        ('test_data/9-0-413610210-0', 1.3461),
    ])
def test_compute_simulation_length(dir_name, true_density):
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert approx_equal(compute_density(trajectory), true_density)


@pytest.mark.parametrize(
    'dir_name,true_degree_polymerization',
    [
        ('test_data/9-0-246295613-0', 14),
        ('test_data/9-0-413610210-0', 26),
    ])
def test_compute_simulation_length(dir_name, true_degree_polymerization):
    metadata = get_metadata(dir_name)
    trajectory = LammpsTrajectoryLoader().load(dir_name)
    assert compute_degree_polymerization(trajectory, **metadata) == true_degree_polymerization
