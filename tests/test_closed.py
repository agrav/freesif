# -*- coding: utf-8 -*-
"""Test that closed File and SifData instances behave as expected.
"""

import os
import unittest
import freesif as fs


FILES = os.path.join(os.path.dirname(__file__), 'files')


hydrodata_methods_and_args = [('get_addedmass', ()),
                              ('get_angular_freqs', ()),
                              ('get_bodymass', ()),
                              ('get_bodyproperties', ()),
                              ('get_directions', ()),
                              ('get_excitationforce_raos', ()),
                              ('get_fluidkinematics_raos', ()),
                              ('get_hydrostatic_restoring', ()),
                              ('get_mooring_restoring', ()),
                              ('get_motion_raos', ()),
                              ('get_periods', ()),
                              ('get_points', ()),
                              ('get_potentialdamping', ()),
                              ('get_sectionforce_raos', ()),
                              ('get_sections', ()),
                              ('get_timesteps', ()),
                              ('get_viscousdamping', ()),
                              ('get_meandrift', ()),
                              ('get_horiz_meandrift', ())]

strucdata_methods_and_args = [('get_setnames', ()),
                              ('get_nodes', ()),
                              ('get_nodenumbers', ()),
                              ('get_noderesults', ('displacement',)),
                              ('get_elements', ()),
                              ('get_elementnumbers', ()),
                              ('get_elementresults', ('generalstress',))]

file_methods_and_args = [('__contains__', ('key',)),
                         ('__getitem__', ('key',)),
                         ('__iter__', ()),
                         ('get', ('key',)),
                         ('items', ()),
                         ('keys', ()),
                         ('values', ())]


class TestHydroData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # establish a closed HydroData instance
        cls._data = fs.open_sif(os.path.join(FILES, 'hydro', 'slowdrift_G1.SIF'))
        cls._data.close()

    def test_methods(self):
        # all methods should raise ClosedFileError
        for method, args in hydrodata_methods_and_args:
            self.assertRaises(fs.exceptions.ClosedFileError,
                              getattr(self._data, method),
                              *args)

    def test_close(self):
        # calling close() on an already closed SifData should return None
        self.assertIsNone(self._data.close())


class TestStrucData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # establish a closed StrucData instance (single sup. elem.)
        cls._data = fs.open_sif(os.path.join(FILES, 'struc', 'single_super_elem', 'test01_1stord_linstat_R1.SIU'))
        cls._data.close()

    def test_methods(self):
        # all methods should raise ClosedFileError
        for method, args in strucdata_methods_and_args:
            self.assertRaises(fs.exceptions.ClosedFileError,
                              getattr(self._data, method),
                              *args)

    def test_close(self):
        # calling close() on an already closed SifData should return None
        self.assertIsNone(self._data.close())


class TestFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # establish a closed File instance
        cls._fname = fs.sif2hdf5(os.path.join(FILES, 'hydro', 'slowdrift_G1.SIF'))
        cls._file = fs.open_hdf5(cls._fname)
        cls._file.close()

    @classmethod
    def tearDownClass(cls):
        os.remove(cls._fname)

    def test_methods(self):
        # all methods should raise ClosedFileError
        for method, args in file_methods_and_args:
            self.assertRaises(fs.exceptions.ClosedFileError,
                              getattr(self._file, method),
                              *args)

    def test_close(self):
        # calling close() on an already closed File should return None
        self.assertIsNone(self._file.close())


if __name__ == '__main__':
    unittest.main()
