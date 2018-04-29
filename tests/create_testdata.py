# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 07:51:19 2018

@author: audun
"""

# create datasets (arrays and .vtu files) for verification
# against Xtract and as subsequent basis for unit testing

# f_data['test01_1stord_linstat_R1/LD_plates/elemres/nodes']
# f_data['test01_1stord_linstat_R1/LD_plates/elemres/elems']
# f_data['test01_1stord_linstat_R1/LD_plates/elemres/generalstress']
# f_data['test01_1stord_linstat_R1/LD_plates/noderes/nodes']
# f_data['test01_1stord_linstat_R1/LD_plates/noderes/elems']
# f_data['test01_1stord_linstat_R1/LD_plates/noderes/displacement']

import freesif as fs
import h5py

# file to store datasets
f_data = h5py.File('files/verified_testdata.h5', 'w')


### single superelem, linear static, 1st order, simple LC ###

# open sif-data
data = fs.open_sif('files/struc/single_super_elem/test01_1stord_linstat_R1.SIU')

# LD plates, elemres
gr = f_data.create_group('test01_1stord_linstat_R1/LD_plates/elemres')
nodes = data.get_nodes(sets='LD_plates', kind='shell', disconnected=True)
elems = data.get_elements(sets='LD_plates', kind='shell', disconnected=True)
generalstress = data.get_elementresults('generalstress', rescases=1, sets='LD_plates')

gr.create_dataset('nodes', data=nodes)
gr.create_dataset('connectivity', data=elems[0])
gr.create_dataset('offset', data=elems[1])
gr.create_dataset('eltyp', data=elems[2])
gr.create_dataset('generalstress', data=generalstress)

vtu = fs.utils.VtuWriter('files/vtu/test01_1stord_linstat_R1_LD_plates_elemres.vtu')
vtu.new_piece(nodes, elems)
vtu.start_pointdata()
vtu.add_data_cellpoints(generalstress[:,0,:], 'generalstress_lower',compnames=['SIGXX', 'SIGYY', 'TAUXY'])
vtu.add_data_cellpoints(generalstress[:,1,:], 'generalstress_upper',compnames=['SIGXX', 'SIGYY', 'TAUXY'])

# LD plates, noderes
gr = f_data.create_group('test01_1stord_linstat_R1/LD_plates/noderes')
nodes = data.get_nodes(sets='LD_plates', kind='shell', disconnected=False)
elems = data.get_elements(sets='LD_plates', kind='shell', disconnected=False)
displacement = data.get_noderesults('displacement', rescases=1, sets='LD_plates',
                                    disconnected=False)
gr.create_dataset('nodes', data=nodes)
gr.create_dataset('connectivity', data=elems[0])
gr.create_dataset('offset', data=elems[1])
gr.create_dataset('eltyp', data=elems[2])
gr.create_dataset('displacement', data=displacement)

vtu = fs.utils.VtuWriter('files/vtu/test01_1stord_linstat_R1_LD_plates_noderes.vtu')
vtu.new_piece(nodes, elems)
vtu.start_pointdata()
vtu.add_data(displacement, 'displacement',compnames=['X', 'Y', 'Z', 'RX', 'RY', 'RZ'])
vtu.add_data(displacement[:,:3], 'displacement_vec',compnames=['X', 'Y', 'Z'])

# LD beams, elemres
sets = ['LD_chords', 'LD_beams']
gr = f_data.create_group('test01_1stord_linstat_R1/LD_beams/elemres')
nodes = data.get_nodes(sets=sets, kind='beam', disconnected=True)
elems = data.get_elements(sets=sets, kind='beam', disconnected=True)
beamforce = data.get_elementresults('beamforce', rescases=1, sets=sets)

gr.create_dataset('nodes', data=nodes)
gr.create_dataset('connectivity', data=elems[0])
gr.create_dataset('offset', data=elems[1])
gr.create_dataset('eltyp', data=elems[2])
gr.create_dataset('beamforce', data=beamforce)

vtu = fs.utils.VtuWriter('files/vtu/test01_1stord_linstat_R1_LD_beams_elemres.vtu')
vtu.new_piece(nodes, elems)
vtu.start_pointdata()
vtu.add_data_cellpoints(beamforce, 'beamforce',
                        compnames=['NXX', 'NXY', 'NXZ','MXX', 'MXY', 'MXZ'])


# plate element data quads
#nodes_quad4 = data.get_nodes(sets='MD_plates', kind='shell', disconnected=True)
#elems_quad4 = data.get_elements(sets='MD_plates', kind='shell', disconnected=True)
#generalstress_quad4 = data.get_elementresults('generalstress', rescases=1, sets='MD_plates')
#




#vtu.new_piece(nodes_quad4, elems_quad4)
#vtu.start_pointdata()
#vtu.add_data_cellpoints(generalstress_quad4[:,0,:], 'generalstress_lower',
#                        compnames=['sxx', 'syy', 'sxy'])

#vtu.new_piece(nodes_bm2, elems_bm2)
#vtu.start_pointdata()
#vtu.add_data_cellpoints(beamforce_bm2, 'beamforce',
#                        compnames=['NXX', 'NXY', 'NXZ','MXX', 'MXY', 'MXZ'])

vtu.close()




data.close()


### single superelem, linear static, 2nd order, simple LC ###







f_data.close()
