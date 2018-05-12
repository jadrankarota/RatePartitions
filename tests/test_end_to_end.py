import os
from unittest import TestCase

from rate_partitions import run


class TestRatePartitions(TestCase):
    def setUp(self):
        self.test_data_path = os.path.dirname(__file__) + "/TestData/"

    def test_run(self):
        """End-to-end test

        Test that the output data of running the rate_partitions script is the
        expected file contents
        """
        infile = self.test_data_path + "AS1-1.phy_r8s.txt"
        outfile = self.test_data_path + "AS1-1.phy_r8s.txt_2.5.txt"
        divnum = 2.5
        result = run(infile, divnum)

        with open(outfile) as handle:
            expected_result = handle.read()
            self.assertEqual(expected_result, result)

