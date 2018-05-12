# -*- coding: utf-8 -*-
"""Read and write vtf files
"""

import numpy as np
from itertools import count
from collections import OrderedDict
from ..exceptions import VtfError

# ---Global variables------------------------------------------------------- #

block_types = {'*NODES', '*ELEMENTS', '*RESULTS', '*GLVIEWSCALAR',
               '*GLVIEWGEOMETRY'}

nodes_directives = {'%NO_ID', '%WITH_ID'}

elem_directives = {'%NAME', '%DESCRIPTION', '%NODES', '%NO_ID',
                   '%WITH_ID'}

elem_types = {'%BEAMS', '%BEAMS_3', '%TRIANGLES','%TRIANGLES_6',
              '%QUADS', '%QUADS_8'}

res_directives = {'%DIMENSION', '%WITH_ID', '%NO_ID'}

res_types = {'%PER_NODE', '%PER_ELEMENT', '%PER_ELEMENT_NODE',
             '%PER_ELEMENT_FACE', '%PER_ELEMENT_FACE_NODE',
             '%PER_FACE'}

geom_directives = {'%NAME', '%DESCRIPTION'}

scalar_directives = {'%NAME', '%DESCRIPTION', '%STEP'}

# ---Global dicts----------------------------------------------------------- #

elem_type_map = {15: '%BEAMS',
                 23: '%BEAMS_3',
                 25: '%TRIANGLES',
                 26: '%TRIANGLES_6',
                 24: '%QUADS',
                 28: '%QUADS_8'}

number_of_nodes = {'%BEAMS':2,
                   '%BEAMS_3':3,
                   '%TRIANGLES':3,
                   '%TRIANGLES_6':6,
                   '%QUADS':4,
                   '%QUADS_8':8}

# ---Helper functions------------------------------------------------------- #

def _reformatdata():
    """Convert data from freesif format to format used here.
    """

#                              freesif:                 here:
#    nodes:             2d array (nnodes, 3)      2d array (nnodes, 3/4)
#
#    elements:
#
#    element results:

    pass

# ---Data classes----------------------------------------------------------- #

class Data(object):
    pass


class Nodes(Data):
    """
    """

    def __init__(self, nodes_arr, nodeid_arr=None):
        self.nodes_arr = nodes_arr
        self.nodeid_arr = nodeid_arr
        self.with_id = False if nodeid_arr is None else True


class Scalar(Data):
    """
    """

    def __init__(self, name, description, results_ids=[]):
        self.name = name
        self.description = description
        self.results_ids = results_ids
        # self.step ??


class Elements(Data):
    """
    """

    def __init__(self, name, description, nodes_id, elems_dict, elem_ids=None):
        self.name = name
        self.description = description
        self.nodes_id = nodes_id
        self.elems_dict = elems_dict
        self.with_id = True if elem_ids else False
        self.elem_ids = elem_ids
        self.results_ids = []


class Results(Data):
    """
    """

    # should this be ElementResults?

    def __init__(self, res_type, res_dict, elems_id, dimension):
        self.res_type = res_type
        self.res_dict = res_dict
        self.elems_id = elems_id
        self.dimension = dimension


# ---Main class------------------------------------------------------------- #


class VtfData(object):
    """Read/Write data from/to vtf-files.

    Beam and shell results should generally be written to separate files.
    """

    def __init__(self, name='Model', description=''):

        self.geom_id = 1
        self.name = name
        self.description = description

        self.nodes_id_counter = count(1)
        self.nodes = OrderedDict()
        self.current_nodes_id = None

        self.elems_id_counter = count(1)
        self.elems = OrderedDict()

        self.scalar_id_counter = count(1)
        self.scalars = OrderedDict()

        # should probably have a 'node_results' and a 'element_results'..
        self.results_id_counter = count(1)
        self.results = OrderedDict()

        self.counters = {'nodes': self.nodes_id_counter,
                         'elems': self.elems_id_counter,
                         'scalars': self.scalar_id_counter,
                         'results': self.results_id_counter}

        self.containers = {'nodes': self.nodes,
                           'elems': self.elems,
                           'scalars': self.scalar,
                           'results': self.results}

        self.current_ids = {'nodes': self.current_nodes_id}
        # TODO: finish this with elems and scalars

    def _get_new_id(self, objtype):
        obj_id = self.counters[objtype].next()
        while obj_id in self.containers[objtype].keys():
            obj_id = self.counters[objtype].next()

    def _get_current_id(self, objtype):
        pass

    def read(self, filename):

        lines = open(filename).readlines()

        lines_iter = iter(lines)

        # handle blank lines
        def get_next_line():
            line = lines_iter.next()
            while line.strip() == '':  # empty line
                line = lines_iter.next()
            return line

        line = get_next_line()

        try:
            while True:

                # skip to first block
                while line.split()[0] not in block_types:
                    line = get_next_line()
#                line = _skip_lines(line)


                # encountered a block definition
                block_name, block_id = line.split()
                block_id = int(block_id)

                if block_name == '*NODES':
                    line = get_next_line()
                    with_id = False  # Default

                    # NO_ID/WITH_ID directive if any
                    toks = line.split()
                    if toks[0] == '%NO_ID':
                        line = get_next_line()
                    elif toks[0] == '%WITH_ID':
                        with_id = True
                        line = get_next_line()

                    # read node data until next block def
                    node_list = []
                    nodeid_list = []
                    while line.split()[0] not in block_types:
                        if with_id:
                            toks = line.split()
                            node_list.append(map(float, toks[1:]))
                            nodeid_list.append(int(toks[0]))
                        else:
                            node_list.append(map(float, line.split()))
                        line = get_next_line()

                    nodes_arr = np.array(node_list)
                    if not with_id:
                        nodeid_arr = None
                    else:
                        nodeid_arr = np.array(nodeid_list)

                    self.nodes[block_id] = Nodes(nodes_arr, nodeid_arr)

                elif block_name == '*ELEMENTS':
                    name = ''
                    description = ''
                    nodes_id = 0
                    elems_dict = OrderedDict()
                    elem_ids = OrderedDict()
                    with_id = False

                    line = get_next_line()

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
                        line = get_next_line()

                    # read element defs for each element type
                    while line.strip() in elem_types:
                        elem_type = line.strip()
                        line = get_next_line()
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
                            line = get_next_line()

                        elems_dict[elem_type] = np.array(elem_list, np.int32)
                        elem_ids[elem_type] = np.array(elem_id_list, np.int32)

                    if not with_id:
                        elem_ids = None

                    self.elems[block_id] = Elements(name,
                                                    description,
                                                    nodes_id,
                                                    elems_dict,
                                                    elem_ids)

                elif block_name == '*RESULTS':
                    res_type = ''
                    res_dict = OrderedDict()
                    elems_id = 0
                    dimension = 1
                    with_id = False

                    line = get_next_line()

                    # read directives
                    while line.split()[0] in res_directives.union(res_types):
                        toks = line.split()
                        if toks[0] == '%DIMENSION':
                            dimension = int(toks[1])
                        elif toks[0] == '%WITH_ID':  # only for node results
                            with_id = True
                        elif toks[0] in res_types:
                            res_type = toks[0]
                            elems_id = int(toks[1][1:])

                        line = get_next_line()

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
                                    line = get_next_line()
                                res.append(elem_res)

                            # create array
                            res_arr = np.array(res, dtype=np.float32)
#                            res_arr.shape = (1, len(arr), n_res, 1, 1)

                            # add to res_dict
                            res_dict[elem_type] = res_arr

                            # check if end of block
                            if line.split()[0] in block_types:
                                break

                        self.results[block_id] = Results(res_type[1:],
                                                         res_dict,
                                                         elems_id,
                                                         dimension)

                    else:
                        raise NotImplementedError(
                            'result type {} not implemented'.format(res_type))


                elif block_name == '*GLVIEWSCALAR':
                    line = get_next_line()
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
                            line = get_next_line()
                            res_ids = [int(x.strip()) for x in line.split(',')]
                        try:
                            line = get_next_line()
                        except StopIteration:
                            self.scalars[block_id] = Scalar(
                                name, description, res_ids)
                            raise

                    self.scalars[block_id] = Scalar(
                        name, description, res_ids)


                elif block_name == '*GLVIEWGEOMETRY':
                    self.geom_id = block_id
                    line = get_next_line()
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
                        line = get_next_line()


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

    def add_nodes(self, nodes, node_numbers=None, nodes_id=None):
        """
        """
        if not nodes_id:
            nodes_id = self._get_new_id('nodes')
        self.current_nodes_id = nodes_id
        self.nodes[nodes_id] = Nodes(nodes, node_numbers)

    def add_elements(self, name, elements, elem_numbers=None,
                     description='', nodes_id=None, elems_id=None):
        """
        """
        has_elem_numbers = True if elem_numbers is not None else False

        if not nodes_id:
#            nodes_id = self.current_nodes_id
            nodes_id = self._get_current_id()
        if not elems_id:
            elems_id = self.elems_id_counter.next()
            while elems_id in self.elems.keys():
                elems_id = self.elems_id_counter.next()

        # need to convert data format:

        con, offset, eltyp = elements

        # change con if nodenumbers
        if self.nodes[nodes_id].with_id:
            con = self.nodes[nodes_id].nodeid_arr.take(con)

        unique_eltyp = np.unique(eltyp)
        elems_dict = OrderedDict()
        elem_ids = OrderedDict()
        for et in unique_eltyp:
            elems_dict[et] = []
            elem_ids[et] = []

        o0 = 0
        for o1, et in zip(offset, eltyp):
            elems_dict[et].append(con[o0:o1])
            o0 = o1

        for et in elems_dict.keys():
            elems_dict[elem_type_map[et]] = np.array(elems_dict[et])
            if not self.nodes[nodes_id].with_id:
                elems_dict[elem_type_map[et]] += 1
            del elems_dict[et]

        if has_elem_numbers:
            o0 = 0
            for elno, o1, et in zip(elem_numbers, offset, eltyp):
                elem_ids[et].append(elno)
                o0 = o1

            for et in elem_ids.keys():
                elem_ids[elem_type_map[et]] = np.array(elem_ids[et])
                del elem_ids[et]

        # need to reorder node sequence for certain element types
        for et, arr in elems_dict.items():
            if et == '%QUADS_8':
                arr[...] = np.concatenate((arr[:,:7:2], arr[:,1::2]), axis=1)


        # store offset and eltyp in Elements object
        # (needed when adding element results)

        self.current_elems_id = elems_id
        if not has_elem_numbers:
            elem_ids = None
        self.elems[elems_id] = Elements(name, description, nodes_id,
                                        elems_dict, elem_ids)

    def add_scalar(self, name, description='', scalar_id=None):
        """
        """

        # currently (time)steps are not implemented...

        if not scalar_id:
            scalar_id = self.scalar_id_counter.next()

        self.scalars[scalar_id] = Scalar(name, description)
        self.current_scalar_id = scalar_id

        return scalar_id

    def add_element_scalar_results(self, res_type, results, elems_id=None,
                                   results_id=None, scalar_id=None,
                                   scalar_name=None, scalar_description=None):
        """results_id, scalar_id, scalar_name and scalar_description can either
        be int, int, str and str respectively or sequences of those types. If
        multiple result components are given, i.e. the provided results array
        have shape (..., ncomp), then the sequences must have lenght ncomp. If
        scalar_id(s) are provided, scalar_name and scalar_description shall not
        be provided, and vice versa. If neither of these are provided, default
        behaviour is to use the current scalar_id, if one exists, and the
        provided results are for a single component, otherwise new scalars with
        default names are created. scalar_description is optional.
        """

        # Relevant result types (res_type):
        # %PER_ELEMENT: results.shape = (nres, nelems, ncomps)
        # %PER_ELEMENT_NODE: results.shape = (nres, nrespts, ncomps)
        # %PER_ELEMENT_FACE: results.shape = (nres, nelems, 2, ncomps)
        # %PER_ELEMENT_FACE_NODE: results.shape = (nres, nrespts, 2, ncomps)

        # add result cases (nres) as steps ?

        # alternatives:

        # - scalar_id provided as single int:
        #   - results array must contain result data for a single component
        #
        # - Nor scalar_id or scalar_name is provided:
        #   - if single result component:
        #     - try scalar_id = self.current_scalar_id
        #     - else scalar_id = self.add_scalar('rescomp1')
        #     - register results_id with scalar and element group...
        #     - create Results object
        #   - if multiple result components:
        #     - loop comp axis and:
        #       - get new scalar_id (self.add_scalar)
        #       - get new results_id (self._get_new_id('results'))
        #       - register results_id with scalar and element group...
        #       - create Results object
        #
        # - scalar_name (and opt scalar_description) provided as single str:
        #   - results array must contain result data for a single component


        if not results_id:
            results_id = self.results_id_counter.next()
            # need to make sure the id number is not taken:
            while results_id in self.results.keys():
                results_id = self.results_id_counter.next()

         # assiciate results_id with elements
        if not elems_id:
            elems_id = self.current_elems_id
        self.elems[elems_id].results_ids.append(results_id)

        # associate results_id with scalar
        if not scalar_id:
            scalar_id = self.current_scalar_id
        self.scalars[scalar_id].results_ids.append(results_id)


        # add results data
        self.results[results_id] = Results(res_type, results, elems_id, 1)

    def write(self, filename):

        f_out = open(filename, 'w')

        # header
        f_out.write('*VTF-1.00\n')

        # nodes
        for nodes_id, nodes_data in self.nodes.iteritems():
            f_out.write('*NODES {}\n'.format(nodes_id))
            if nodes_data.with_id:
                f_out.write('%WITH_ID\n')
                string ='{:<8.0f}' +  3*'{:>12.5f}' + '\n'
                for nodeid, coords in zip(nodes_data.nodeid_arr,
                                          nodes_data.nodes_arr):
                    f_out.write(string.format(nodeid, *coords))
            else:
                f_out.write('%NO_ID\n')
                string = 3*'{:>12.5f}' + '\n'
                for coords in nodes_data.nodes_arr:
                    f_out.write(string.format(*coords))


        # elements
        for elems_id, elems_data in self.elems.iteritems():
            f_out.write('*ELEMENTS {}\n'.format(elems_id))
            f_out.write('%NAME "{}"\n'.format(elems_data.name))
            f_out.write('%DESCRIPTION "{}"\n'.format(elems_data.description))
            f_out.write('%NODES #{}\n'.format(elems_data.nodes_id))
            f_out.write('%WITH_ID\n' if elems_data.with_id else '%NO_ID\n')

            # in order to write element results for only certain element types,
            # e.g. only beams or only shells, the order in which data for the
            # different element types are written matter. The element types
            # which have results should be written first, and in the same order
            # as on the associated results dicts.

            # get order for which element types shall be written
            if elems_data.results_ids:

                result_keys = self.results[elems_data.results_ids[0]] \
                    .res_dict.keys()

                # verify that all res_dicts associated with
                # this element block has identical key lists (elem types):
                for results_id in elems_data.results_ids:
                    if self.results[results_id].res_dict.keys() != result_keys:
                        raise VtfError(
                            'Element block "{}": All results associated with a'
                            ' single element block should apply to the same '
                            'element types. Basically don\'t mix beam and '
                            'shell results for the same element block in the '
                            'same vtf file.'.format(elems_data.name))

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
            if elems_data.with_id:
                for elem_ids_arr, (elem_type, elems_arr) in zip(
                    elems_data.elem_ids.itervalues(),
                    ordered_elems_dict.iteritems()):
                    # element type:
                    f_out.write('{}\n'.format(elem_type))

                    # write elements
                    nnodes = elems_arr.shape[1]
                    string = (1 + nnodes)*'{:<10}' + '\n'
                    for elem_id, node_refs in zip(elem_ids_arr, elems_arr):
                        f_out.write(string.format(elem_id, *node_refs))
            else:
                for elem_type, elems_arr in ordered_elems_dict.iteritems():
                    # element type:
                    f_out.write('{}\n'.format(elem_type))

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
