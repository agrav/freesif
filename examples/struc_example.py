# -*- coding: utf-8 -*-
"""
Created on Thu May 24 07:36:25 2018

@author: audun
"""

import freesif as fs

# convert to HDF5
h5name = fs.sif2hdf5('../tests/files//struc/single_super_elem/test01_2ndord_linstat_R1.SIU')

# open HDF5 file
f = fs.open_hdf5(h5name)  # f is a 'File' object

# get 'StrucData' object
data = f['test01_2ndord_linstat_R1']

# get node and element data for shells
nodes = data.get_nodes()
elems = data.get_elements()

# get element results
res = data.get_noderesults('displacement', rescases=1)


# write data to Paraview
vtu = fs.utils.VtuWriter('vonmises.vtu')
vtu.new_piece(nodes, elems)
vtu.start_pointdata()
vtu.add_data(res[:,:3], 'displacement')  # write x,y,z
vtu.close()

f.close()
