# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""
"""

import numpy as np
from struct import pack
from base64 import b64encode

from ..exceptions import VtkError


# index arrays to reorder element connectivity - keys are SESAM element id's
sesam2vtk_connectivity = {23 : np.array([0,2,1]),
                          28 : np.array([0,2,4,6,1,3,5,7])}

# index arrays to reorder element node results - keys are SESAM element id's
sesam2vtk_elemresults = {26 : np.array([0,2,4,1,3,5])}

# get vtk element type id
sesam2vtk_elementtypeid = {15 : 3,
                           23 : 21,
                           24 : 9,
                           25 : 5,
                           26 : 22,
                           28 : 23}

# get vtk ids from sesam ids
typeidmapper = np.vectorize(lambda x: sesam2vtk_elementtypeid[x])

# get vtk datatype name from numpy datatype name
numpy2vtk_datatype = {'int8'    : 'Int8',
                      'uint8'   : 'UInt8',
                      'int16'   : 'Int16',
                      'uint16'  : 'UInt16',
                      'int32'   : 'Int32',
                      'uint32'  : 'UInt32',
                      'int64'   : 'Int64',
                      'uint64'  : 'UInt64',
                      'float32' : 'Float32',
                      'float64' : 'Float64'}

def _reorder_data(data, offsets, types, indexarr_map):
    start = 0
    for stop, typ in zip(offsets, types):
        indarr = indexarr_map.get(typ)
        if indarr is not None:
            data[start:stop] = data[start:stop][indarr]
        start = stop

def _reorder_connectivity_vtk(cells):
    """Perform inplace reordering of connectivity array
    """
    con, offsets, types = cells
    _reorder_data(con, offsets, types, sesam2vtk_connectivity)

def _reorder_cellpointdata_vtk(resultarr, offsets, types):
    """Perform inplace reordering of element node results
    """
    _reorder_data(resultarr, offsets, types, sesam2vtk_elemresults)


class VtuWriter(object):
    """write data to vtk xml format (UnstructuredGrid)

    Examples
    --------
    ::

        import freesif as fs
        f = fs.open_hdf5('R1.h5')
        t1 = f['T1']
        nodes = t1.get_nodes()
        elems = t1.get_elements()
        disp = t1.get_noderesults('displacement')

        vtu = fs.utils.vtuwriter('out.vtu')
        vtu.new_piece(nodes, elems)
        vtu.start_pointdata()
        vtu.add_data(disp, 'displacement')
        vtu.close()
        f.close()

    """

    def __init__(self, filename):
        """
        """
        self._file = open(filename, 'wb')
        self._piece_open = False
        self._pointdata_open = False
        self._celldata_open = False
        self._elements = []
        self._current_offsets = None
        self._current_types = None

        self._open_element('VTKFile', type='UnstructuredGrid', version='1.0',
                           byte_order='LittleEndian', header_type='UInt64')
        self._open_element('UnstructuredGrid')


    def new_piece(self, points, cells):
        """
        """
        self._close_data()
        self._close_piece()
        self._piece_open = True

        _reorder_connectivity_vtk(cells)
        connectivity, offsets, types = cells

        # keep reference to offsets and types
        self._current_offsets = offsets
        self._current_types = types

        types = typeidmapper(types)

        self._open_element('Piece', NumberOfPoints=len(points),
                           NumberOfCells=len(offsets))

        # write points
        self._open_element('Points')
        self._write_array(points, 'points')
        self._close_element('Points')

        # write cells
        self._open_element('Cells')
        self._write_array(connectivity, 'connectivity')
        self._write_array(offsets, 'offsets')
        self._write_array(types, 'types')
        self._close_element('Cells')

    def start_pointdata(self):
        self._open_data('PointData')

    def start_celldata(self):
        self._open_data('CellData')

    def add_data(self, arr, name, compnames=None):
        """
        """
        if not self._pointdata_open and not self._celldata_open:
            raise VtkError('Call start_pointdata() or start_celldata() first')
        self._write_array(arr, name, compnames)

    def add_data_cellpoints(self, arr, name, compnames=None):
        """
        """
        _reorder_cellpointdata_vtk(
            arr, self._current_offsets, self._current_types)
        self.add_data(arr, name, compnames)


    def close(self):
        self._close_data()
        self._close_piece()
        self._close_element('UnstructuredGrid')
        self._close_element('VTKFile')
        self._file.close()

    def _write_array(self, arr, name, compnames=None):

        # accept 1d or 2d array
        if arr.ndim == 2:
            ncomps = arr.shape[1]
        elif arr.ndim == 1:
            ncomps = 1
        else:
            raise VtkError('array must be 1d or 2d')

        # create comp name kwargs
        if compnames:  # sequence of str
            compname_dict = {'ComponentName{}'.format(i): cn \
                             for i, cn in enumerate(compnames)}
        else:
            compname_dict = {}

        vtktype = numpy2vtk_datatype[arr.dtype.name]
        self._open_element('DataArray', type=vtktype, Name=name,
                           NumberOfComponents=ncomps, format='binary',
                           **compname_dict)

        # write indent
        self._file.write(b'  '*len(self._elements))

        # write header
        nbytes = arr.size * arr.dtype.itemsize
#        self._file.write(b64encode(pack('<Q', nbytes)))

        # write array
#        self._file.write(b64encode(arr.T.ravel().tostring(order='F')))

        header = pack('<Q', nbytes)
        arrstr = arr.ravel().tostring()
        self._file.write(b64encode(header + arrstr))

        # write newline
        self._file.write(b'\n')

        self._close_element('DataArray')

    def _close_piece(self):
        if self._piece_open:
            self._close_element('Piece')
            self._piece_open = False

    def _open_data(self, datatype):
        if not self._piece_open:
            raise VtkError('Call new_piece() first')
        self._close_data()
        self._open_element(datatype)
        if datatype == 'PointData':
            self._pointdata_open = True
        elif datatype == 'CellData':
            self._celldata_open = True

    def _close_data(self):
        if self._pointdata_open:
            self._close_element('PointData')
            self._pointdata_open = False
        elif self._celldata_open:
            self._close_element('CellData')
            self._celldata_open = False

    def _open_element(self, tag, selfclosing=False, **kwargs):
        s = b'  '*len(self._elements)  # indentation
        s += b'<%s' % tag.encode()
        if kwargs:
            s += b' '
            s += b' '.join([b'%s="%s"' % (str(k).encode(), str(v).encode())
                for k, v in list(kwargs.items())])
        if selfclosing:
            s += b' />\n'
        else:
            s += b'>\n'
            self._elements.append(tag)
        self._file.write(s)

    def _close_element(self, tag):
        if not tag == self._elements.pop():
            raise VtkError('{} has no corresponding opening tag'.format(tag))
        s = b'  '*len(self._elements)
        s += b'</%s>\n' % tag.encode()
        self._file.write(s)



