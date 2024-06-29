from sage.all import *
from index_calc import solveDLP
from helpers import modpow
import unittest


class UnitTests(unittest.TestCase):
    def test_basic(self):
        """
        Most basic test: 1^0 = 1
        """
        self.assertEqual(0, solveDLP(1, 1, 2, 5))


    def test_small(self):
        """
        Small numbers test: easily verifiable by eye
        """
        self.assertEqual(3, solveDLP(3, 6, 7, 5))


    def test_textbook(self):
        """
        The test given from Example 3.58 in HPS:
        """
        self.assertEqual(8500, solveDLP(37, 211, 18443, 5))


    def test_sq_free_raises_exception(self):
        """
        Assert that breaking the assumption that $p-1$ is square-free raises an exception
        """
        p = 2^2 * 3 + 1 
        with self.assertRaises(Exception):
            solveDLP(3, 6, p, 5)

    def test_large(self):
        """
        Large Example
        """
        p = 14087 
        g = 5 
        x = 11608 
        h = modpow(g, x, p)
        self.assertEqual(x, solveDLP(g, h, p, 15))


