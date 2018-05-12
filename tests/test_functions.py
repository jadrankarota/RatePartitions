import os
from unittest import TestCase

from rate_partitions import read_input_file


class TestFunctions(TestCase):
    def setUp(self):
        self.test_data_path = os.path.dirname(__file__) + "/TestData/"
        self.input_file = self.test_data_path +  "AS1-1.phy_r8s.txt"

    def test_read_input_file(self):
        result = read_input_file(self.input_file)
        self.assertIn(0.222972, result)

