import pytest
from os.path import join, dirname
from htpmd.analysis import get_all_properties

TEST_DATA_FOLDER = join(dirname(__file__), 'test_data')


@pytest.mark.parametrize(
    'dir_name',
    [
        (join(TEST_DATA_FOLDER, '9-0-246295613-0')),
        (join(TEST_DATA_FOLDER, '9-0-413610210-0')),
    ])
def test_get_all_properties(dir_name):
    results = get_all_properties(dir_name)
    assert results is not None
