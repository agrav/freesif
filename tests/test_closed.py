# -*- coding: utf-8 -*-
"""Test that closed File and SifData instances behave as expected.
"""

import unittest
import freesif as fs

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

file_methods_and_args = [('__contains__', ('key',)),
                         ('__getitem__', ('key',)),
                         ('__iter__', ()),
                         ('get', ('key',)),
                         ('items', ()),
                         ('iteritems', ()),
                         ('iterkeys', ()),
                         ('itervalues', ()),
                         ('keys', ()),
                         ('values', ())]


class TestHydroData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # establish a closed HydroData instance
        cls._data = fs.open_sif('../test_files/slowdrift_G1.SIF')
        cls._data.close()

#    def setUp(self):
#
#        # establish a closed HydroData instance
#        self._data = fs.open_sif('../test_files/slowdrift_G1.SIF')
#        self._data.close()

    def test_methods(self):
        # all methods should raise ClosedFileError
        for method, args in hydrodata_methods_and_args:
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
        cls._file = fs.open_hdf5('../test_files/slowdrift_G1.h5')
        cls._file.close()

#    def setUp(self):
#        # establish a closed File instance
#        self._file = fs.open_hdf5('../test_files/slowdrift_G1.h5')
#        self._file.close()

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
