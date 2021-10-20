import unittest
from ..py.combine_open_lists import mothurList
#import ..code.py.uc_to_list


class testMothurList(unittest.TestCase):

    def test_combine(self):
        ml1 = mothurList('userLabel', ['a,b', 'c,d', 'e'])
        ml2 = mothurList('userLabel', ['f,g', 'h'])
        ml1.combine(ml2)
        self.assertEqual(ml1, mothurList('userLabel', ['a,b', 'c,d', 'e', 'f,g', 'h']))


class testUcList(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
