# -*- coding: utf-8 -*-
"""Test sif2hdf5 function for different filetypes.
"""

import unittest
import freesif as fs
import os


class TestSIF2HDF5(unittest.TestCase):
    """Test sif2hdf5 function for different filetypes. It is only checked that
    the correct file is created, content on the file is not verified. Content
    is verified in the HydroData and StrucData tests.
    """

    @classmethod
    def setUpClass(cls):
        # file names
        cls._in_files = {'siu_single': './files/struc/single_super_elem/test01_2ndord_linstat_R1.SIU',
                         'siu_assembly': './files/struc/assembly/R100.SIU',
                         'fem_single': './files/struc/single_super_elem/test01_2ndord_linstat_T1.FEM',
                         'fem_assembly': './files/struc/assembly/T100.FEM',
                         'sif_hydro': './files/hydro/slowdrift_G1.SIF'}
        cls._out_files = {'siu_single': './files/tmp/siu_single_R1.h5',
                          'siu_assembly': './files/tmp/siu_assembly_R100.h5',
                          'fem_single': './files/tmp/fem_single_T1.h5',
                          'fem_assembly': './files/tmp/fem_assembly_T100.h5',
                          'sif_hydro': './files/tmp/sif_hydro_G1.h5'}

        # create /tmp folder if it doesn't exist
        if not os.path.isdir('./files/tmp'):
            os.mkdir('./files/tmp')

        # remove h5-files if exist
        for fname in cls._out_files.values():
            if os.path.isfile(fname):
                os.remove(fname)

    @classmethod
    def tearDownClass(cls):
        # remove h5-files
        for fname in cls._out_files.values():
            if os.path.isfile(fname):
                os.remove(fname)

    def test_SIU_single(self):
        fs.sif2hdf5(self._in_files['siu_single'], hdf5name=self._out_files['siu_single'])
        self.assertTrue(os.path.isfile(self._out_files['siu_single']))

    def test_SIU_assembly(self):
        fs.sif2hdf5(self._in_files['siu_assembly'], hdf5name=self._out_files['siu_assembly'])
        self.assertTrue(os.path.isfile(self._out_files['siu_assembly']))

    def test_FEM_single(self):
        fs.sif2hdf5(self._in_files['fem_single'], hdf5name=self._out_files['fem_single'])
        self.assertTrue(os.path.isfile(self._out_files['fem_single']))

    def test_FEM_assembly(self):
        fs.sif2hdf5(self._in_files['fem_assembly'], hdf5name=self._out_files['fem_assembly'])
        self.assertTrue(os.path.isfile(self._out_files['fem_assembly']))

    def test_SIF_hydro(self):
        fs.sif2hdf5(self._in_files['sif_hydro'], hdf5name=self._out_files['sif_hydro'])
        self.assertTrue(os.path.isfile(self._out_files['sif_hydro']))


if __name__ == '__main__':
    unittest.main()
