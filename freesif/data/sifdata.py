# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""Base class for all data classes.
"""

import tables as tb
import numpy as np

from ..sequentialparser.parsers import atoms
from ..exceptions import NoSuchRecordError, ClosedFileError

class SifData(object):
    """Base class for StrucData and HydroData
    """

    def __init__(self, tbgroup, fileinstance=None):
        """
        Parameters
        ----------
        group : tables.Group
            A ``tables.Group`` instance holding Sesam data.
        """

        # check that group actually is part of fileinstance._tbfile?

        self.name = tbgroup._v_name
        self.rec_names = tbgroup._f_getattr('RECORDS')
        self.file = fileinstance
        self._tbgroup = tbgroup
        self._seltyp = tbgroup._f_getattr('SELTYP')
        self._level = tbgroup._f_getattr('LEVEL')
        self._filetype = tbgroup._f_getattr('FILETYPE')
        self._datatype = tbgroup._f_getattr('DATATYPE')

        self._isopen = True

    def __str__(self):
        # brief description

        if not self._isopen:
            return '<Closed SifData object>'

        if self._datatype == 'HYDRODYNAMIC':
            return '{} data'.format(self._datatype.lower())

        else:
            return '{} {} data'.format(self._level.lower(),
                                       self._datatype.lower())

    def __repr__(self):
        # detailed description

        if not self._isopen:
            return '<Closed SifData object>'

        return 'SifData(name={})'.format(self.name)

    def __del__(self):
        pass

    def _check_isopen(self):
        if not self._isopen:
            raise ClosedFileError(
                'Trying to operate on closed SifData object')

    def _get_record(self, rec_name):
        """Get pytables dataset(s) for a record.

        Returns the Table with rec_name records. If arrays exist for this
        record type, a tuple is returned: (rec, arr1, ..., arrn)
        """

        # check if file is open
        self._check_isopen()

        # get main record
        try:
            rec = self._tbgroup._f_get_child(rec_name.lower())
        except tb.NoSuchNodeError:
            raise NoSuchRecordError(
                '"{}" does not have any "{}" records'.format(self.name,
                                                             rec_name.upper()))
        # get any associated arrays
        arr_atoms = atoms.get(rec_name.upper())
        if arr_atoms:
            arrs = []
            for name, atom, shape in arr_atoms:
                arrname = '{}_{}'.format(rec_name.lower(), name)
                arrs.append(self._tbgroup._f_get_child(arrname))
            return (rec,) + tuple(arrs)
        return rec

    def _get_varlendata(self, rec, arr, indices):
        """
        """

        colname = arr.name[len(rec.name)+1:]
        start_indices = rec.col(colname + '_start')[indices]
        stop_indices = rec.col(colname + '_stop')[indices]

        datalist = []
        offsetlist = []
        offset = 0
        for start, stop in zip(start_indices, stop_indices):
            datalist.append(arr[start:stop])
            offset += stop - start
            offsetlist.append(offset)
        return np.concatenate(datalist), np.array(offsetlist)

    def _has_record(self, rec_name):
        """Check if record type is present on this dataset"""
        return rec_name.upper() in self.rec_names

    def _record_info(self):
        infostr = ''
        for t in self._tbgroup._f_iter_nodes('Table'):
            infostr += '{:8}: {}\n'.format(t.name, len(t))
        return infostr

    def close(self):
        if not self._isopen:
            return

        if self.file:
            self.file.close()
        else:
            self._close()

    def _close(self):
        self.__dict__.clear()
        self._isopen = 0






