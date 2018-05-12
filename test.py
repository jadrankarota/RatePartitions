from unittest import TestCase

from rate_partitions import run


class TestRatePartitions(TestCase):
    def test_run(self):
        """End-to-end test

        Test that the output data of running the rate_partitions script is the
        expected file contents
        """
        infile = "TestData/AS1-1.phy_r8s.txt"
        divnum = 2.5
        result = run(infile, divnum)
        with open("TestData/AS1-1.phy_r8s.txt_2.5.txt", "r") as handle:
            expected_result = handle.read()
            self.assertEqual(expected_result, result)

