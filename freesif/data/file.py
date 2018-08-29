# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""Defines a class to operate on SIF data records stored on a HDF5 file.
"""

import tables as tb
from os import path
from .helpers import parse_sifname
from .strucdata import TopLevelData, InterLevelData, FirstLevelData
from .hydrodata import HydroData
from ..sequentialparser.sif2hdf5 import sif2hdf5
from ..exceptions import ClosedFileError


# keep track of open files
_open_files = set()

def open_hdf5(filename, mode='r'):
    """Returns a File object.
    """

    # currently this only return a File, passing the args directly to the
    # constructor. However in the future a factory pattern may be applied,
    # so this function should be used instead of creating a File directly.

    return File(filename, mode)

def open_sif(filename):
    """Convert a SIF file to an in-memory HDF5 file and return the SifData
    object representing the SIF file.

    Suitable for smaller SIF files, e.g. hydrodynamic results, where conversion
    time is short.
    """

    # the in-memory file will be in 'w' (write) mode !
    # this prevents files with the same name (at different paths) to be open
    # at the same time.

    tbfile = sif2hdf5(filename, in_memory=True)
    prefix, letter, idno, ext = parse_sifname(filename)
    groupname = prefix + letter + idno
    return File(tbfile)[groupname]

def close_all_files():
    global _open_files
    for f in _open_files:
        f.close()
    _open_files = set()


class File(object):
    """Class to operate on SIF data stored on a HDF5 file.
    """

    def __init__(self, filename, mode='r'):
        """
        """

        # mode is currently not relevant, but maybe in the future, e.g.:
        # f = fs.open_hdf5('data.h5', 'w')
        # t1 = f.new_strucdata('T1', level=1)  # create a new (empty) StrucData
        # t1.add_nodes(nodes)
        # t1.add_elements(elems)
        # ...
        # t1.write('T1.FEM')

        # accept filemane or tables.File instance
        if isinstance(filename, str):
            # serve the full path to tables to ensure files with the same name
            # at different locations are treated independently
            fullpath = path.abspath(filename)

            # put this in a try/except to clean up if tables raises an
            # exception? probably not necessary as no initialization is done
            self._tbfile = tb.open_file(fullpath, mode)
            self.filename = path.split(fullpath)[1]

        elif isinstance(filename, tb.File):
            self._tbfile = filename
            self.filename = self._tbfile.filename

        self._data = {}
        # loop over independent (top level or non-hierarhcy) groups
        for gr in self._tbfile.root._f_iter_nodes('Group'):

            datatype = gr._f_getattr('DATATYPE')

            if datatype == 'STRUCTURAL':
                level = gr._f_getattr('LEVEL')
                if level == 'TOPLEVEL':
                    data = TopLevelData(gr, self)

                elif level == 'FIRSTLEVEL':
                    data = FirstLevelData(gr, self)

            elif datatype == 'HYDRODYNAMIC':
                data = HydroData(gr, self)

            self._data[gr._v_name] = data

            # loop over subgroups (intermediate or 1st level groups)
            for subgr in gr._f_iter_nodes('Group'):
                level = subgr._f_getattr('LEVEL')
                key = gr._v_name + '/' + subgr._v_name
                if level == 'INTERLEVEL':
                    self._data[key] = InterLevelData(subgr, data)

                elif level == 'FIRSTLEVEL':
                    self._data[key] = FirstLevelData(subgr, data)

        _open_files.add(self)
        self._isopen = True

    def __getitem__(self, name):  # dict interface
        self._check_isopen()
        return self._data[name]

    def get(self, name, d=None):  # dict interface
        self._check_isopen()
        return self._data.get(name, d)

    def keys(self):  # dict interface
        self._check_isopen()
        return self._data.keys()

    def values(self):  # dict interface
        self._check_isopen()
        return self._data.values()

    def items(self):  # dict interface
        self._check_isopen()
        return self._data.items()

    def __contains__(self, name):  # dict interface
        self._check_isopen()
        return name in self._data

    def __str__(self):
        # brief description

        if not self._isopen:
            return '<Closed File object>'

        s = ''
        for k, v in sorted(self._data.items()):
            s += '{!r}: {!s}\n'.format(k, v)
        return s

    def __repr__(self):
        # detailed description

        if not self._isopen:
            return '<Closed File object>'

        s = ''
        for d in sorted(self._data.values()):
            s += '{!r}\n'.format(d)
        return s

    def __iter__(self):  # NOT dict behaviour! values instead of keys
        return self.values()

    def _check_isopen(self):
        if not self._isopen:
            raise ClosedFileError(
                'Trying to operate on closed File object')

    def close(self):
        """
        """
        if not self._isopen:
            return

        for data in self._data.values():
            data._close()

        self._tbfile.close()

        _open_files.remove(self)
        filename = self.filename
        self.__dict__.clear()
        self._isopen = False
        self.filename = filename

#    @classmethod
#    def from_sif(cls, sifname, usefileprefix=True, prefix='', only=None):
#        """
#        """
#        import time
#        hdf5name = 'sifdata_{}.h5'.format(time.strftime('%Y%m%d_%H%M%S'))
#        sif2hdf5(sifname, hdf5name, False, usefileprefix, prefix, only)
#        return cls(hdf5name)
