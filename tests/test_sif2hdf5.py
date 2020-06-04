# -*- coding: utf-8 -*-
"""
Test sif2hdf5 function for different filetypes.
"""
import unittest
import freesif as fs
import os
import shutil


FILES = os.path.join(os.path.dirname(__file__), 'files')


class TestSIF2HDF5(unittest.TestCase):
    """Test *sif2hdf5* function for different filetypes. It is only checked that
    the correct file is created, that it opens with the *open_hdf5* function and
    that the type and number of records on the file is correct. The actual data
    on the files is not verified, this is done in the HydroData and StrucData
    tests.
    """

    @classmethod
    def setUpClass(cls):
        # file names
        cls._in_files = dict(
            siu_single=os.path.join(FILES, 'struc', 'single_super_elem', 'test01_2ndord_linstat_R1.SIU'),
            siu_assembly=os.path.join(FILES, 'struc', 'assembly', 'R100.SIU'),
            fem_single=os.path.join(FILES, 'struc', 'single_super_elem', 'test01_2ndord_linstat_T1.FEM'),
            fem_assembly=os.path.join(FILES, 'struc', 'assembly', 'T100.FEM'),
            sif_hydro=os.path.join(FILES, 'hydro', 'slowdrift_G1.SIF')
        )
        cls._out_files = dict(
            siu_single=os.path.join(FILES, 'tmp', 'siu_single_R1.h5'),
            siu_assembly=os.path.join(FILES, 'tmp', 'siu_assembly_R100.h5'),
            fem_single=os.path.join(FILES, 'tmp', 'fem_single_T1.h5'),
            fem_assembly=os.path.join(FILES, 'tmp', 'fem_assembly_T100.h5'),
            sif_hydro=os.path.join(FILES, 'tmp', 'sif_hydro_G1.h5')
        )

        # create a clean /tmp directory
        try:
            shutil.rmtree(os.path.join(FILES, 'tmp'))
        except FileNotFoundError:
            # ok, so it was already removed
            pass
        finally:
            # create empty directory
            os.mkdir(os.path.join(FILES, 'tmp'))

    @classmethod
    def tearDownClass(cls):
        # remove '\tmp' directory with H5 files
        try:
            shutil.rmtree(os.path.join(FILES, 'tmp'))
        except FileNotFoundError:
            # ok, so it was already removed
            pass

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

        # TODO: check that type/number of records on h5 file is correct inside
        # the following tests

    def test_open_SIU_single(self):
        f = fs.open_hdf5(self._out_files['siu_single'])
        f.close()

    def test_open_SIU_assembly(self):
        f = fs.open_hdf5(self._out_files['siu_assembly'])
        f.close()

    def test_open_FEM_single(self):
        f = fs.open_hdf5(self._out_files['fem_single'])
        f.close()

    def test_open_FEM_assembly(self):
        f = fs.open_hdf5(self._out_files['fem_assembly'])
        f.close()

    def test_open_SIF_hydro(self):
        f = fs.open_hdf5(self._out_files['sif_hydro'])
        f.close()


if __name__ == '__main__':
    unittest.main()
