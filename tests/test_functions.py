import os
from unittest import TestCase
from mock import patch

from rate_partitions import read_input_file, verify_divnum, clean_string


class TestFunctions(TestCase):
    def setUp(self):
        self.test_data_path = os.path.dirname(__file__) + "/TestData/"
        self.input_file = self.test_data_path +  "AS1-1.phy_r8s.txt"

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

    def test_clean_string(self):
        result = clean_string("this text,,[]'")
        expected = "this text,,"
        self.assertEqual(expected, result)

    def test_clean_string__extra_char(self):
        result = clean_string("this text,,[]'", additional_char=",")
        expected = "this text"
        self.assertEqual(expected, result)
