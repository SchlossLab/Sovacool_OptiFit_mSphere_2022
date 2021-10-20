import unittest


class testMothurList(unittest.TestCase):
    def test_combine(self):
        self.assertEqual(['a', 'b'] + ['c'], ['a', 'b', 'c'])


if __name__ == "__main__":
    unittest.main()
