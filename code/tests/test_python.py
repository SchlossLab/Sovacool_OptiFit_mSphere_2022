import unittest
from ..py.combine_open_lists import mothurList
from ..py.uc_to_list import uc_to_list


class testMothurList(unittest.TestCase):
    def test_combine(self):
        ml1 = mothurList("userLabel", ["a,b", "c,d", "e"])
        ml2 = mothurList("userLabel", ["f,g", "h"])
        ml1.combine(ml2)
        self.assertEqual(ml1, mothurList("userLabel", ["a,b", "c,d", "e", "f,g", "h"]))

    def test_write(
        self,
        list_filename_test="code/tests/data/test_open.list",
        list_filename_oracle="code/tests/data/oracle_open.list",
    ):
        ml1 = mothurList("userLabel", ["a,b", "c,d", "e", "f,g", "h"])
        ml1.write(list_filename_test)

        with open(list_filename_test, "r") as infile:
            test_list = infile.readlines()
        with open(list_filename_oracle, "r") as infile:
            oracle_list = infile.readlines()

        self.assertEqual(test_list, oracle_list)


class testUcList(unittest.TestCase):
    def test_convert(
        self,
        uc_filename="code/tests/data/closed.uc",
        list_filename_test="code/tests/data/test_closed.list",
        list_filename_oracle="code/tests/data/oracle_closed.list",
    ):
        uc_to_list(uc_filename, list_filename_test)

        with open(list_filename_test, "r") as infile:
            test_list = infile.readlines()
        with open(list_filename_oracle, "r") as infile:
            oracle_list = infile.readlines()

        self.assertEqual(test_list, oracle_list)


if __name__ == "__main__":
    unittest.main()
