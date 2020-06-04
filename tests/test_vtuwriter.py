# -*- coding: utf-8 -*-
"""
Created on Sun May 13 22:22:07 2018

@author: audun
"""

import unittest
import freesif as fs
import os


FILES = os.path.join(os.path.dirname(__file__), 'files')


class TestVtuWriter(unittest.TestCase):
    """Test the VtuWriter class. It is only checked that the correct file is
    created without errors.
    """

    @classmethod
    def setUpClass(cls):
        # file names
        cls._in_files = dict(siu_single=os.path.join(FILES, 'struc', 'single_super_elem', 'test01_2ndord_linstat_R1.SIU'))
        cls._out_files = dict(siu_single=os.path.join(FILES, 'tmp', 'siu_single.vtu'))

        # create /tmp folder if it doesn't exist
        if not os.path.isdir(os.path.join(FILES, 'tmp')):
            os.mkdir(os.path.join(FILES, 'tmp'))

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

    def test_SIU_single_noderesults(self):
        """Write some node results for whole model with disconnected=False
        """
        data = fs.open_sif(self._in_files['siu_single'])
        nodes = data.get_nodes()
        elements = data.get_elements()
        disp = data.get_noderesults('displacement', rescases=1)
        vtu = fs.utils.VtuWriter(self._out_files['siu_single'])
        vtu.new_piece(nodes, elements)
        vtu.start_pointdata()
        vtu.add_data(disp, 'displacement')
        vtu.close()
        self.assertTrue(os.path.isfile(self._out_files['siu_single']))


if __name__ == '__main__':
    unittest.main()

