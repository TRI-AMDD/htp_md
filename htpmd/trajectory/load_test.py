import pytest

from htpmd.trajectory.load import get_metadata


@pytest.mark.parametrize(
    'dir_name',
    [
        ('test_data/9-0-246295613-0'),
        ('test_data/9-0-413610210-0'),
    ])
def test_get_metadata(dir_name):
    metadata_list = [
        'mol_smiles', 'poly_smiles', 'force_field', 'material_group',
        'temperature', 'time_step']
    meta_data = get_metadata(dir_name)
    assert set(metadata_list) == set(list(meta_data.keys()))