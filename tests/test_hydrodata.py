# -*- coding: utf-8 -*-
"""Test HydroData public methods
"""

import unittest
import numpy as np
import freesif as fs


class TestHydroData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # establish HydroData instances and associated verification data
        cls._data1 = fs.open_sif('../test_files/slowdrift_G1.SIF')
        cls._data1_verified = np.load(
            '../test_files/slowdrift_G1_verified_arrs.npz')

    @classmethod
    def tearDownClass(cls):
        cls._data1.close()
        cls._data1_verified.close()

    def test_periods(self):
        pers = self._data1.get_periods()[::-1]
        pers_verified = self._data1_verified['periods']
        self.assertTrue(np.allclose(pers, pers_verified))

    def test_directions_deg(self):
        dirs = self._data1.get_directions()  # 'degrees' is default
        dirs_verified = self._data1_verified['directions_deg']
        self.assertTrue(np.allclose(dirs, dirs_verified))

    def test_directions_rad(self):
        dirs = self._data1.get_directions(unit='radians')
        dirs_verified = self._data1_verified['directions_rad']
        self.assertTrue(np.allclose(dirs, dirs_verified))

if __name__ == '__main__':
    unittest.main()
