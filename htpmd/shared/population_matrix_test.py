import pytest
import os
import numpy as np
from htpmd.shared.population_matrix import compute_pop_matrix_from_pizza_dump


def approx_equal(val1, val2):
    return pytest.approx(val1, rel=1.0e-6) == val2


@pytest.mark.parametrize(
    'dir_name',
    [
        ('test_data/9-0-246295613-0'),
    ])
def test_compute_pop_matrix_from_pizza_dump(dir_name):
    pop_matrix = compute_pop_matrix_from_pizza_dump(dir_name, 90, 94)
    path_to_pop_matrix = os.path.join(dir_name, 'population.txt')
    assert approx_equal(pop_matrix, np.loadtxt(path_to_pop_matrix))
