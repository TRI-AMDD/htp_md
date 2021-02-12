import json
import os
import pytest
from htpmd.analysis import get_all_properties


@pytest.mark.parametrize(
    'dir_name',
    [
        ('test_data/9-0-246295613-0'),
        ('test_data/9-0-413610210-0'),
    ])
def test_get_all_properties(dir_name):
    results = get_all_properties(dir_name)
    property_list = [
        'li_diffusivity', 'tfsi_diffusivity', 'poly_diffusivity',
        'conductivity', 'molarity', 'li_msd_curve', 'tfsi_msd_curve',
        'structure', 'mol_smiles', 'poly_smiles', 'force_field',
        'material_group', 'temperature', 'time_step']
    assert set(property_list) == set(list(results.keys()))