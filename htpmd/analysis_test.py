import unittest
from htpmd.analysis import get_all_properties


class TestGetAllProperties(unittest.TestCase):
    def test_success_run(self):
        dir_lists = [
            'test_data/9-0-246295613-0',
            'test_data/9-0-413610210-0']
        for dir_name in dir_lists:
            with self.subTest():
                get_all_properties(dir_name)


if __name__ == '__main__':
    unittest.main()
