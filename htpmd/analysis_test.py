import unittest
import json
import os
from htpmd.analysis import get_all_properties


class TestGetAllProperties(unittest.TestCase):
    def test_success_run(self):
        dir_lists = [
            '9-0-246295613-0',
            '9-0-413610210-0']
        for dir_name in dir_lists:
            with self.subTest():
                results = get_all_properties(
                    os.path.join('test_data', dir_name))
                with open(os.path.join(
                    'test_data', 'test_results', dir_name + '.json')) as f:
                    correct_results = json.load(f)
                # Only compare the subset that is recorded
                subset_results = {k: results[k] for k in correct_results.keys()}
                self.assertDictEqual(subset_results, correct_results)


if __name__ == '__main__':
    unittest.main()
