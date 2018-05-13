import os
from unittest import TestCase
from mock import patch

from rate_partitions import read_input_file, verify_divnum, generate_sites, \
    generate_sites_last_partition


class TestFunctions(TestCase):
    def setUp(self):
        self.test_data_path = os.path.dirname(__file__) + "/TestData/"
        self.input_file = self.test_data_path +  "AS1-1.phy_r8s.txt"
        self.input_data = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

    def test_read_input_file(self):
        result = read_input_file(self.input_file)
        self.assertIn(0.222972, result)
        self.assertEqual(4000, len(result))

    @patch("rate_partitions.sys.exit")
    def test_verify_divnum__bad_input(self, mock_sys):
        """Test when the given divfactor is lower than 1.1

        We will exit with error and print a message for the user
        """
        verify_divnum(0.001)
        self.assertTrue(mock_sys.called)

    @patch("rate_partitions.sys.exit")
    def test_verify_divnum__good_input(self, mock_sys):
        """Test when the given divfactor is greater than 1.1

        We will not exit with error
        """
        verify_divnum(2.5)
        self.assertFalse(mock_sys.called)

    def test_generate_sites__first_partition(self):
        # upper_value is the max_rate of input data
        result = generate_sites(
            self.input_data,
            partition=1,
            upper_value=max(self.input_data),
            min_rate=min(self.input_data),
            divnum=2.5,
        )
        expected_sites, expected_lower_value = [1, 2, 3, 4], 0.64
        self.assertEqual(expected_sites, result[0])
        self.assertEqual(expected_lower_value, result[1])

    def test_generate_sites_last_partition(self):
        result = generate_sites_last_partition(
            self.input_data, lower_value=min(self.input_data), upper_value=0.64)
        expected = [5, 6, 7, 8, 9, 10]
        self.assertEqual(expected, result)

    def test_generate_sites__middle_partition(self):
        # assume this is a second partition
        result = generate_sites(
            self.input_data,
            partition=2,
            upper_value=0.64,
            min_rate=min(self.input_data),
            divnum=2.5,
        )
        expected_sites, expected_lower_value = [5, 6], 0.47
        self.assertEqual(expected_sites, result[0])
        self.assertAlmostEqual(expected_lower_value, round(result[1], 2))
