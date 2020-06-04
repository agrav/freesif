# -*- coding: utf-8 -*-
"""Test HydroData public methods
"""

import os
import unittest
import numpy as np
import freesif as fs


FILES = os.path.join(os.path.dirname(__file__), 'files')


class TestHydroData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # establish HydroData instances and associated verification data
        cls._data1 = fs.open_sif(os.path.join(FILES, 'hydro', 'slowdrift_G1.SIF'))
        cls._data1_verified = np.load(os.path.join(FILES, 'hydro', 'slowdrift_G1_verified_arrs.npz'))

    @classmethod
    def tearDownClass(cls):
        cls._data1.close()
        cls._data1_verified.close()

    def test_periods(self):
        pers = self._data1.get_periods()
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

    def test_motion_raos(self):
        res = self._data1.get_motion_raos()
        res_verified = self._data1_verified['motion_raos']
        self.assertTrue(np.allclose(res, res_verified))

    def test_addedmass(self):
        res = self._data1.get_addedmass()
        res_verified = self._data1_verified['addedmass']
        self.assertTrue(np.allclose(res, res_verified))

    def test_angular_freqs(self):
        res = self._data1.get_angular_freqs()
        res_verified = self._data1_verified['angular_freqs']
        self.assertTrue(np.allclose(res, res_verified))

    def test_bodymass(self):
        res = self._data1.get_bodymass()
        res_verified = self._data1_verified['bodymass']
        self.assertTrue(np.allclose(res, res_verified))

    def test_excitationforce_raos(self):
        res = self._data1.get_excitationforce_raos()
        res_verified = self._data1_verified['excitationforce_raos']
        self.assertTrue(np.allclose(res, res_verified))

    def test_horiz_meandrift(self):
        res = self._data1.get_horiz_meandrift()
        res_verified = self._data1_verified['horiz_meandrift']
        self.assertTrue(np.allclose(res, res_verified))

    def test_hydrostatic_restoring(self):
        res = self._data1.get_hydrostatic_restoring()
        res_verified = self._data1_verified['hydrostatic_restoring']
        self.assertTrue(np.allclose(res, res_verified))

    def test_meandrift(self):
        res = self._data1.get_meandrift()
        res_verified = self._data1_verified['meandrift']
        self.assertTrue(np.allclose(res, res_verified))

    def test_mooring_restoring(self):
        res = self._data1.get_mooring_restoring()
        res_verified = self._data1_verified['mooring_restoring']
        self.assertTrue(np.allclose(res, res_verified))

    def test_potentialdamping(self):
        res = self._data1.get_potentialdamping()
        res_verified = self._data1_verified['potentialdamping']
        self.assertTrue(np.allclose(res, res_verified))

    def test_viscousdamping(self):
        res = self._data1.get_viscousdamping()
        res_verified = self._data1_verified['viscousdamping']
        self.assertTrue(np.allclose(res, res_verified))


if __name__ == '__main__':
    unittest.main()
