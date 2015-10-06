# -*- coding: utf-8 -*-
"""Read and write vtf files
"""

import numpy as np
from itertools import count
from collections import namedtuple, OrderedDict
#from elementresults import ElementResults

# ---Global variables------------------------------------------------------- #

block_types = set(['*NODES', '*ELEMENTS', '*RESULTS', '*GLVIEWSCALAR',
                   '*GLVIEWGEOMETRY'])

nodes_directives = set(['%NO_ID', '%WITH_ID'])

elem_directives = set(['%NAME', '%DESCRIPTION', '%NODES', '%NO_ID',
                       '%WITH_ID'])

elem_types = set(['%BEAMS', '%BEAMS_3', '%TRIANGLES','%TRIANGLES_6',
                  '%QUADS', '%QUADS_8'])

res_directives = set(['%DIMENSION', '%WITH_ID', '%NO_ID'])

res_types = set(['%PER_NODE', '%PER_ELEMENT', '%PER_ELEMENT_NODE',
                 '%PER_ELEMENT_FACE', '%PER_ELEMENT_FACE_NODE',
                 '%PER_FACE #ID'])

geom_directives = set(['%NAME', '%DESCRIPTION'])

scalar_directives = set(['%NAME', '%DESCRIPTION', '%STEP'])

# ---Global dicts----------------------------------------------------------- #

# element name map, elem_name_map[sesam_name] = vtf_name
elem_name_map = {'BEAS': '%BEAMS',
                 'BTSS': '%BEAMS_3',
                 'FTRS': '%TRIANGLES',
                 'SCTS': '%TRIANGLES_6',
                 'FQUS': '%QUADS',
                 'SCQS': '%QUADS_8'}

# and the other way: elem_name_map[vtf_name] = sesam_name
elem_name_map2 = {y:x for x,y in elem_name_map.items()}

number_of_nodes = {'%BEAMS':2,
                   '%BEAMS_3':3,
                   '%TRIANGLES':3,
                   '%TRIANGLES_6':6,
                   '%QUADS':4,
                   '%QUADS_8':8}

# ---Helper functions------------------------------------------------------- #

def _reformatdata():
    """Convert data from freesif format to format used here.
                              freesif:                 here:
    nodes:             2d array (nnodes, 3)      2d array (nnodes, 3/4)

    elements:
    element results:

    """
    pass

# ---Main class------------------------------------------------------------- #


class vtfdata(object):
    """Read data from and to vtf-files.

    Beam and shell results should be written to separate files.
    """

    # why named tuple??

    Nodes = namedtuple('Nodes', ['nodes_arr', 'nodeid_arr', 'with_id'])
    Scalar = namedtuple('Scalar', ['name', 'description', 'results_ids'])
    Elements = namedtuple('Elements', ['name', 'description', 'nodes_id',
                                       'elems_dict', 'results_ids',
                                       'elem_ids', 'with_id, id'])
    Results = namedtuple('Results', ['res_type', 'res_dict', 'elems_id',
                                     'dimension'])

    def __init__(self, name='Model', description=''):

        self.geom_id = 1
        self.name = name
        self.description = description

        self.nodes_id_counter = count(1)
        self.nodes = OrderedDict()

        self.elems_id_counter = count(1)
        self.elems = OrderedDict()

        self.scalar_id_counter = count(1)
        self.scalars = OrderedDict()

        # should probably have a 'node_results' and a 'element_results'..
        self.results_id_counter = count(1)
        self.results = OrderedDict()

    def read(self, filename):

        # will currently fail if blank lines are present in file...

        lines = open(filename).readlines()

        lines_iter = iter(lines)
        line = lines_iter.next()

        try:
            while True:

                # skip to first block
                while line.split()[0] not in block_types:
                    line = lines_iter.next()
#                line = _skip_lines(line)


                # encountered a block definition
                block_name, block_id = line.split()
                block_id = int(block_id)

                if block_name == '*NODES':
                    line = lines_iter.next()

                    # skip NO_ID/WITH_ID directive if any
                    if line.split()[0] in nodes_directives:
                        line = lines_iter.next()

                    # read node data until next block def
                    node_list = []
                    while line.split()[0] not in block_types:
                        node_list.append([float(x) for x in line.split()])
                        line = lines_iter.next()

                    self.nodes[block_id] = np.array(node_list, np.float32)

                elif block_name == '*ELEMENTS':
                    name = ''
                    description = ''
                    nodes_id = 0
                    elems_dict = OrderedDict()
                    elem_ids = OrderedDict()
                    with_id = False
                    results_ids = []

                    line = lines_iter.next()

                    # read directives
                    while line.split()[0] in elem_directives:
                        toks = line.split()
                        if toks[0] == '%NAME':
                            if len(toks) == 2:
                                name = toks[1][1:-1]
                            else:
                                name = ' '.join(toks[1:])[1:-1]
                        elif toks[0] == '%DESCRIPTION':
                            if len(toks) == 2:
                                description = toks[1][1:-1]
                            else:
                                description = ' '.join(toks[1:])[1:-1]
                        elif toks[0] == '%NODES':
                            nodes_id = int(toks[1][1:])
                        elif toks[0] == '%WITH_ID':
                            with_id = True
                        line = lines_iter.next()

                    # read element defs for each element type
                    while line.strip() in elem_types:
                        elem_type = line.strip()
                        line = lines_iter.next()
                        elem_list = []
                        elem_id_list = []
                        stop_if_in = elem_types.union(block_types)
                        while line.split()[0] not in stop_if_in:
                            toks = line.split()
                            if with_id:
                                elem_id_list.append(toks[0])
                                elem_list.append(toks[1:])
                            else:
                                elem_list.append(toks)
                            line = lines_iter.next()

                        elems_dict[elem_type] = np.array(elem_list, np.int32)
                        elem_ids[elem_type] = np.array(elem_id_list, np.int32)

                    self.elems[block_id] = vtfdata.Elements(name,
                                                             description,
                                                             nodes_id,
                                                             elems_dict,
                                                             results_ids,
                                                             elem_ids)

                elif block_name == '*RESULTS':
                    res_type = ''
                    res_dict = OrderedDict()
                    elems_id = 0
                    dimension = 1
                    with_id = False

                    line = lines_iter.next()

                    # read directives
                    while line.split()[0] in res_directives.union(res_types):
                        toks = line.split()
                        if toks[0] == '%DIMENSION':
                            dimension = int(toks[1])
                        elif toks[0] == '%WITH_ID':  # only for nodes ?
                            with_id = True
                        elif toks[0] in res_types:
                            res_type = toks[0]
                            elems_id = int(toks[1][1:])

                        line = lines_iter.next()

                    # read results for element types in the referenced
                    # elements block.
                    # assume only one res_type per block
                    # assume scalar results

                    if res_type == '%PER_ELEMENT_NODE':
                        # get elements
                        elems = self.elems[elems_id]
                        elems.results_ids.append(block_id)

                        for elem_type, arr in elems.elems_dict.iteritems():
                            n_res = number_of_nodes[elem_type]
                            resinds = range(n_res)
                            res = []
                            for i in arr:
                                elem_res = []
                                for j in resinds:
                                    elem_res.append(float(line))
                                    line = lines_iter.next()
                                res.append(elem_res)

                            # create array
                            res_arr = np.array(res, dtype=np.float32)
#                            res_arr.shape = (1, len(arr), n_res, 1, 1)

                            # add to res_dict
                            res_dict[elem_type] = res_arr

                            # check if end of block
                            if line.split()[0] in block_types:
                                break

                        self.results[block_id] = vtfdata.Results(res_type[1:],
                            res_dict, elems_id, dimension)

                    else:
                        raise NotImplementedError('result type {} not '
                                                  'implemented'.format(
                                                  res_type))


                elif block_name == '*GLVIEWSCALAR':
                    line = lines_iter.next()
                    # read directives
                    while line.split()[0] in scalar_directives:
                        toks = line.split()
                        if toks[0] == '%NAME':
                            if len(toks) == 2:
                                name = toks[1][1:-1]
                            else:
                                name = ' '.join(toks[1:])[1:-1]
                        elif toks[0] == '%DESCRIPTION':
                            if len(toks) == 2:
                                description = toks[1][1:-1]
                            else:
                                description = ' '.join(toks[1:])[1:-1]
                        elif toks[0] == '%STEP':
                            # read result block ids (comma separated)
                            line = lines_iter.next()
                            res_ids = [int(x.strip()) for x in line.split(',')]
                        try:
                            line = lines_iter.next()
                        except StopIteration:
                            self.scalars[block_id] = vtfdata.Scalar(
                                name, description, res_ids)
                            raise

                    self.scalars[block_id] = vtfdata.Scalar(
                        name, description, res_ids)


                elif block_name == '*GLVIEWGEOMETRY':
                    self.geom_id = block_id
                    line = lines_iter.next()
                    # read directives
                    while line.split()[0] in geom_directives:
                        toks = line.split()
                        if toks[0] == '%NAME':
                            if len(toks) == 2:
                                self.name = toks[1][1:-1]
                            else:
                                self.name = ' '.join(toks[1:])[1:-1]
                        elif toks[0] == '%DESCRIPTION':
                            if len(toks) == 2:
                                self.description = toks[1][1:-1]
                            else:
                                self.description = ' '.join(toks[1:])[1:-1]
                        line = lines_iter.next()


        except StopIteration:
            pass

#    def get_nodes(self, node_id=None):
#        if node_id:
#            return self.nodes[node_id]
#        else:
#            return self.nodes.values()[0]
#
#    def get_elements(self, elem_id=None):
#        if elem_id:
#            return self.elems[elem_id]
#        else:
#            return self.elems.values()[0]
#
#    def get_results(self, res_id=None):
#        if res_id:
#            return self.results[res_id]
#        else:
#            return self.results.values()[0]
#
#    def get_scalars(self):
#        pass

    def add_nodes(self, nodes_arr, nodes_id=None):

        if not nodes_id:
            nodes_id = self.nodes_id_counter.next()
            # need to make sure the id number is not taken:
            while nodes_id in self.nodes.keys():
                nodes_id = self.nodes_id_counter.next()

        self.current_nodes_id = nodes_id
        self.nodes[nodes_id] = nodes_arr

    def add_elements(self, name, elems_dict, description='', geom_name=None,
                     nodes_id=None, elems_id=None):

        if not nodes_id:
            nodes_id = self.current_nodes_id
        if not elems_id:
            elems_id = self.elems_id_counter.next()
            # need to make sure the id number is not taken:
            while elems_id in self.elems.keys():
                elems_id = self.elems_id_counter.next()

        self.current_elems_id = elems_id
        self.elems[elems_id] = vtfdata.Elements(name, description, nodes_id,
                                                 elems_dict, [])

    def add_scalar(self, name, description='', scalar_id=None):
        """name is either a str or a sequence of str's.
        """

        if not scalar_id:
            scalar_id = self.scalar_id_counter.next()

        self.scalars[scalar_id] = vtfdata.Scalar(name, description, [])
        self.current_scalar_id = scalar_id

        return scalar_id

    def add_element_scalar_results(self, res_type, res_dict, results_id=None,
                                   elems_id=None, scalar_id=None):
        """The required shape of arrays in res_dict will depend on res_type.
        """

        if not results_id:
            results_id = self.results_id_counter.next()
            # need to make sure the id number is not taken:
            while results_id in self.results.keys():
                results_id = self.results_id_counter.next()

        # associate results_id with scalar
        if not scalar_id:
            scalar_id = self.current_scalar_id
        self.scalars[scalar_id].results_ids.append(results_id)

         # assiciate results_id with elements
        if not elems_id:
            elems_id = self.current_elems_id
        self.elems[elems_id].results_ids.append(results_id)

        # add results data
        self.results[results_id] = vtfdata.Results(res_type, res_dict,
                                                    elems_id, 1)

    def write(self, filename):

        f_out = open(filename, 'w')

        # header
        f_out.write('*VTF-1.00\n')

        # nodes
        for nodes_id, nodes_arr in self.nodes.iteritems():
            f_out.write('*NODES {}\n'.format(nodes_id))
            f_out.write('%WITH_ID\n')
            string ='{:<8.0f}' +  3*'{:>12.5f}' + '\n'
            for coords in nodes_arr:
                f_out.write(string.format(*coords))

        # elements
        for elems_id, elems_data in self.elems.iteritems():
            f_out.write('*ELEMENTS {}\n'.format(elems_id))
            f_out.write('%NAME "{}"\n'.format(elems_data.name))
            f_out.write('%DESCRIPTION "{}"\n'.format(elems_data.description))
            f_out.write('%NODES #{}\n'.format(elems_data.nodes_id))
            f_out.write('%NO_ID\n')

            # get order for which element types shall be written
            if elems_data.results_ids:

                # TODO:
                # here, we should verify that all res_dicts associated with
                # this element block has identical key lists (elem types).
                # skip for now, and just use first ref..
                result_keys = self.results[elems_data.results_ids[0]] \
                    .res_dict.keys()

                # create an OrderedDict and add element types for which
                # results exist:
                ordered_elems_dict = OrderedDict()
                for key in result_keys:
                    ordered_elems_dict[key] = elems_data.elems_dict[key]

                # then add element types without results
                existing_keys = ordered_elems_dict.keys()
                for key, value in elems_data.elems_dict.iteritems():
                    if key not in existing_keys:
                        ordered_elems_dict[key] = value

            else:  # no results associated with this elements block
                ordered_elems_dict = elems_data.elems_dict

            # write element defs in correct order
            for elem_type, elems_arr in ordered_elems_dict.iteritems():
                # element type:
                f_out.write('{}\n'.format(elem_type))

                # re-arrange node order
#                if sesam_name == 'SCQS':
#                    elems_arr = np.concatenate(
#                        (elems_arr[:,:7:2],
#                         elems_arr[:,1::2]),
#                         axis=1)

                # write elements
                nnodes = elems_arr.shape[1]
                string = nnodes*'{:<10}' + '\n'
                for node_refs in elems_arr:
                    f_out.write(string.format(*node_refs))

        # geometry
        f_out.write('*GLVIEWGEOMETRY {}\n'.format(self.geom_id))
        f_out.write('%NAME "{}"\n'.format(self.name))
        f_out.write('%DESCRIPTION "{}"\n'.format(self.description))
        f_out.write('%ELEMENTS\n')
        string = '{}' + (len(self.elems)-1) * ',{}'
        f_out.write(string.format(*self.elems.keys()))
        f_out.write('\n')

        # results
        for results_id, results_data in self.results.iteritems():

            f_out.write('*RESULTS {}\n'.format(results_id))
            f_out.write('%DIMENSION {}\n'.format(results_data.dimension))
            f_out.write('%{} #{}\n'.format(results_data.res_type,
                                          results_data.elems_id))

            # TODO: handle different res_types, currently assumes
            # PER_ELEMENT_NODE?

            for elem_type, res_arr in results_data.res_dict.iteritems():
                for node_values in res_arr:
                    for node_value in node_values:
                        f_out.write('{}\n'.format(node_value))

        # scalars
        for scalar_id, scalar_data in self.scalars.iteritems():
            f_out.write('*GLVIEWSCALAR {}\n'.format(scalar_id))
            f_out.write('%NAME "{}"\n'.format(scalar_data.name))
            f_out.write('%DESCRIPTION "{}"\n'.format(
                        scalar_data.description))
            f_out.write('%STEP 1\n')  # implement several steps per scalar?

            resid_str = ''
            for results_id in scalar_data.results_ids:
                resid_str += '{} '.format(results_id)
            f_out.write(resid_str[:-1] + '\n')

        f_out.close()
