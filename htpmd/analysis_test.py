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
