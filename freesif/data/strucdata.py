# -*- coding: utf-8 -*-
# Copyright (c) 2015-2018 Audun Gravdal Johansen
"""Classes to operate on SESAM structural data.
"""

from __future__ import division

import numpy as np
from .sifdata import SifData
from .helpers import getrow
from ..exceptions import HierarchyError, NoSuchRecordError, ResultError


class StrucData(SifData):
    """Base class for the StrucData classes.
    """

    def __init__(self, tbgroup, fileinstance, toplevel=None):
        """
        Parameters
        ----------
        toplevel : TopLevelData
            If this superelement is part of a superelement hierarchy and is not
            the top level superelement, provide the top level data object
        """

        super(StrucData, self).__init__(tbgroup, fileinstance)
        self._toplevel = toplevel

    def _get_setmember_indices(self, setnames, istype=2):
        """Get (zero-based) indices of set members.

        if istype=1: indices are into the gnode/gcoord tables (nodes),
        if istype=2 (default): indices are into the gelmnt1/gelref1 tables
        (elements)
        """

        tdsetnam, tdsetnam_text = self._get_record('tdsetnam')
        gsetmemb, irmemb = self._get_record('gsetmemb')

        setnames = (setnames,) if isinstance(setnames, str) else setnames
        setmembs = []
        for setname in setnames:
            cond = 'name==b{}'.format(repr(setname))
            isref = getrow(tdsetnam, cond)['idno']
            cond = '(isref=={})&(istype=={})'.format(isref, istype)
            for r in gsetmemb.where(cond):
                setmembs.append(irmemb[r['irmemb_start']:r['irmemb_stop']])
        return np.unique(np.concatenate(setmembs)) - 1

    def _get_elementkind_indices(self, kind, elemindices=None):
        """Get (zero-based) indices of elements of a certain *kind*.

        Legal alternatives for *kind* are 'beam' and 'shell'

        Optionally provide *elemindices* which can be a subset of the complete
        set of element indices.
        """

        gelmnt1, nodin = self._get_record('gelmnt1')

#        # This is slower:
#        if kind == 'beam':
#            cond = '(eltyp==15) | (eltyp==23)'
#        elif kind == 'shell':
#            cond = '(eltyp==24) | (eltyp==25) | (eltyp==26) | (eltyp==28)'
#        else:
#            raise ValueError("kind must be 'beam' or 'shell'")
#        return gelmnt1.get_where_list(cond)

        if elemindices is not None:
            arr = gelmnt1[:][elemindices]
        else:
            arr = gelmnt1[:]

        if kind == 'beam':
            return arr[(arr['eltyp']==15) | (arr['eltyp']==23)]['elno'] - 1
        elif kind == 'shell':
            return arr[(arr['eltyp']==24) | (arr['eltyp']==25) |
                       (arr['eltyp']==26) | (arr['eltyp']==28) ]['elno'] - 1
        else:
            raise ValueError("kind must be 'beam' or 'shell'")

    def _get_elementindices(self, sets=None, kind=None):
        """
        """

        if sets:
            elemindices = self._get_setmember_indices(sets)
        if sets and kind:
            elemindices = self._get_elementkind_indices(kind, elemindices)
        elif kind:
            elemindices = self._get_elementkind_indices(kind)
        elif not sets:
            gelmnt1, nodin = self._get_record('gelmnt1')
            elemindices = gelmnt1.col('elno') - 1

        return elemindices

    def _get_connectivity(self, elemindices):
        """Get connectivity for elements given by *elemindices*.

        Returns a tuple *(con, offset)* where *con* is a 1d array with
        connectivity (nodeindices) for each element given sequentially.
        *offset* is a 1d array with indices into *con* giving the start index
        for each element.
        """

        gelmnt1, nodin = self._get_record('gelmnt1')
        con, offset = self._get_varlendata(gelmnt1, nodin, elemindices)
        con -= 1  # change to zero-based
        return con, offset

    def _update_connectivity(self, con):
        """When elements are a subset (i.e. using sets or kind), the
        connectivity need to be updated to match the node array for the subset.
        """

        nodeindices = np.unique(con)
        indexmap = dict(zip(nodeindices, np.arange(len(nodeindices))))
        indexmapper = np.vectorize(lambda x: indexmap[x])
        return indexmapper(con)

    def _get_nodeindices_ofelems(self, elemindices, disconnected=False):
        """Get (zero-based) indices (into gnode/gcoord table) of nodes used by
        elements given by *elemindices*.
        """
        con, offset = self._get_connectivity(elemindices)
        if disconnected:
            return con
        else:
            return np.unique(con)

    def _get_nodeindices(self, sets=None, kind=None, disconnected=False):
        """Get (zero-based) indices into gnode/gcoord table
        """
        gcoord = self._get_record('gcoord')
        gelmnt1, nodin = self._get_record('gelmnt1')

        if sets or kind:
            elemindices = self._get_elementindices(sets, kind)
            con, _ = self._get_connectivity(elemindices)
#            con = self._update_connectivity(con)
            if disconnected:
                return con
            else:
                return np.unique(con)
        else:
            if disconnected:
                return nodin[:] - 1
            else:
                return np.arange(len(gcoord))

    def _get_resrefs(self, run=1, rescases=None):
        """Returns resref records for the given run and external result case
        number(s) as a numpy structured array. *rescases* can be None (default,
        returns all), int or sequence of ints.
        """

        rdresref, reftyps = self._get_record('rdresref')

        if isinstance(rescases, (list, tuple, np.ndarray)):
            rescases = set(rescases)
        elif isinstance(rescases, int):
            rescases = set([rescases])
        elif rescases is None:
            indices = rdresref.get_where_list('irno=={}'.format(run))
            return rdresref[indices]

        resrefs = []
        for row in rdresref.iterrows():
            if row['irno'] == run and row['ieres'] in rescases:
                resrefs.append(row.fetch_all_fields())
        return np.array(resrefs)

    def _check_complex_and_get_flag(self, complexflags):
        complexflag = complexflags[0]
        if not np.all(complexflags == complexflag):
            raise ResultError(
                'Result cases must all be either real or complex')
        return complexflag

    def get_setnames(self):
        """Get names of all sets in superelement.Returns None if no sets exist
        """
        try:
            return list(self._get_record('tdsetnam')[0].col('name'))
        except NoSuchRecordError:
            pass


class HigherLevelData(StrucData):
    """Base class for higher level structural data
    """
    # this class should contain many of the same methods as FirstLevelData,
    # but with the effect that it calls the method on all children and returns
    # a dict of data..
    # should this dict be nested? or skip intermediate levels
    # is there a difference between input data and results data?
    # should call methods directly on FirstLevelData and use hierarchy data to
    # resolve which objects to call. This to have uniform behaviour for input
    # and results data (as intermediate levels are not physically present on
    # results data)

    # an optional 'only' argument to filter which superelements to get data
    # from
    # the 'sets' argument probably dont make sense.. or? lets keep it: if any
    # of the sets are present on a superelement, data is returned, if not,
    # no data is returned for that superelement

    def get_nodes(self, only=None, sets=None, index=1, trans=None,
                  concepts=False, disconnected=False):
        pass


class TopLevelData(HigherLevelData):
    """
    """
    # hierarchy and transformation data should be kept/gotten here..

    def get_transformation(self, seltyp, index, level='top'):
        """Get transformation matrix from 1stlevel to level(int or 'top')
        """
        hierarch, ihsref = self._get_record('hierarch')
        hsuptran = self._get_record('hsuptran')
        toplevel = hierarch[0]['islevl']

        if level == 'top':
            level = toplevel
        elif level == 1:
            return np.diag(np.ones(4, 'f'))
        elif level > toplevel:
            raise HierarchyError('level {} does not exist, toplevel is '
                                 '{}'.format(level, toplevel))
        elif level < 1:
            raise HierarchyError('level must be 1 or above, got '
                                 '{}'.format(level))

        cond = '(iselty=={})&(indsel=={})'.format(seltyp, index)
        r = getrow(hierarch, cond)
        t = getrow(hsuptran, 'itref=={}'.format(r['itref']))['t']
        rp = getrow(hierarch, 'ihref=={}'.format(r['ihpref']))
        if rp['islevl'] > level:
            raise HierarchyError('Superelement {} with index {} not present on'
                                 ' level {}'.format(seltyp, index, level))
        while rp['islevl'] < level:
            t2 = getrow(hsuptran, 'itref=={}'.format(rp['itref']))['t']
            t = np.dot(t, t2)
            if rp['ihpref'] != 0:
                rp = getrow(hierarch, 'ihref=={}'.format(rp['ihpref']))
                if rp['islevl'] > level:
                    raise HierarchyError('Superelement {} with index {} not '
                                         'present on level {}'.format(seltyp,
                                                                      index,
                                                                      level))
        return t


class InterLevelData(HigherLevelData):
    """
    """
    pass


class FirstLevelData(StrucData):
    """
    """

    # extrapolation matrices for 2nd order elements
    s3 = np.sqrt(3)

    # 8-node quads
    xquad = np.array(
        [[(1+s3)*(1+s3), (1-s3)*(1+s3), (1+s3)*(1-s3), (1-s3)*(1-s3)],
         [(1+s3),        (1+s3),        (1-s3),        (1-s3)       ],
         [(1-s3)*(1+s3), (1+s3)*(1+s3), (1-s3)*(1-s3), (1+s3)*(1-s3)],
         [(1-s3),        (1+s3),        (1-s3),        (1+s3)       ],
         [(1-s3)*(1-s3), (1+s3)*(1-s3), (1-s3)*(1+s3), (1+s3)*(1+s3)],
         [(1-s3),        (1-s3),        (1+s3),        (1+s3)       ],
         [(1+s3)*(1-s3), (1-s3)*(1-s3), (1+s3)*(1+s3), (1-s3)*(1+s3)],
         [(1+s3),        (1-s3),        (1+s3),        (1-s3)       ]]) * .25

    # 6-node triangles
    xtri = np.array([[ 5./3, -1./3, -1./3],
                     [ 2./3,  2./3, -1./3],
                     [-1./3,  5./3, -1./3],
                     [-1./3,  2./3,  2./3],
                     [-1./3, -1./3,  5./3],
                     [ 2./3, -1./3,  2./3]])

    # 3-node beam
    xbm3 = np.array([[1. + s3, 1. - s3],
                     [1.     , 1.     ],
                     [1. - s3, 1. + s3]]) * .5

    # index arrays for 1st order shells
    # (indices of result points located at nodes)
    idxtri = np.array([0, 4, 1, 5, 3, 7])
    idxquad = np.array([0, 5, 1, 6, 4, 9, 3, 8])

    # transformation matrices for decomposision of stresses

    # 1st order plates
    x_decomp_thin = np.array([
            [ 0.5,  0. ,  0. ,  0.5,  0. ,  0. ],
            [ 0. ,  0.5,  0. ,  0. ,  0.5,  0. ],
            [ 0. ,  0. ,  0.5,  0. ,  0. ,  0.5],
            [-0.5,  0. ,  0. ,  0.5,  0. ,  0. ],
            [ 0. , -0.5,  0. ,  0. ,  0.5,  0. ],
            [ 0. ,  0. , -0.5,  0. ,  0. ,  0.5]])

    # 2nd order shells
    x_decomp_thick = np.array([
            [ 0.5,  0. ,  0. ,  0. ,  0. ,  0.5,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0.5,  0. ,  0. ,  0. ,  0. ,  0.5,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0.5,  0. ,  0. ,  0. ,  0. ,  0.5,  0. ,  0. ],
            [-0.5,  0. ,  0. ,  0. ,  0. ,  0.5,  0. ,  0. ,  0. ,  0. ],
            [ 0. , -0.5,  0. ,  0. ,  0. ,  0. ,  0.5,  0. ,  0. ,  0. ],
            [ 0. ,  0. , -0.5,  0. ,  0. ,  0. ,  0. ,  0.5,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  1.5,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  1.5,  0. ,  0. ,  0. ,  0. ,  0. ]])


    def __init__(self, tbgroup, toplevel=None):
        super(FirstLevelData, self).__init__(tbgroup, toplevel)
        if self._filetype == 'R':
            self._index = tbgroup._f_getattr('INDEX')

        # setup dict of element result processing functions
        self._process_eltyp = {15 : self._process_eltyp15,
                               23 : self._process_eltyp23,
                               24 : self._process_eltyp24,
                               25 : self._process_eltyp25,
                               26 : self._process_eltyp26,
                               28 : self._process_eltyp28}

    def _process_element_results(self, eltyp, restype, pos, d):
        """
        """
        return self._process_eltyp[eltyp](restype, pos, d)

    def _process_eltyp15(self, restype, pos, d):
        """BEAS, 3D Beam (2-node)"""
        # it seems this element has 3 result points.
        # return first and last...
        return d.reshape(3, 6)[[0,2],:]

    def _process_eltyp23(self, restype, pos, d):
        """BTSS, General Curved Beam (3-node)"""
        d.shape = (2, 6)
        if pos == 'nodes':
            return np.dot(self.xbm3, d)
        return d

    def _process_eltyp24(self, restype, pos, d):
        """FQUS, Flat Quadrilateral Thin Shell (4-node)"""
        d = d.reshape(10, 3)[self.idxquad,:].reshape(4, 2, 3)
        if restype == 'decomposedstress':
            d = self._decompose_stresses_thin(d)  # shape is (4, 6)
        return d

    def _process_eltyp25(self, restype, pos, d):
        """FTRS, Flat Triangular Thin Shell (3-node)"""
        d = d.reshape(8, 3)[self.idxtri,:].reshape(3, 2, 3)
        if restype == 'decomposedstress':
            d = self._decompose_stresses_thin(d)  # shape is (4, 6)
        return d

    def _process_eltyp26(self, restype, pos, d):
        """SCTS, Subparametric Curved Triangular Thick Shell (6-node)"""
        d.shape = (2, 3, 5)  # (side, respt, comp)
        if pos == 'respts':
            return np.rollaxis(d, 1)
        elif pos == 'nodes':
            d[...,:3] = self._extrapolate_to_surface(d[...,:3])
            if restype == 'generalstress':
                d[...,3:] = 0.
            elif restype == 'decomposedstress':
                d = self._decompose_stresses_thick(d)  # shape is (4, 8)
            # extrapolate to nodes:
            return np.dot(self.xtri, d)  # shape is (6, 2, 5)

    def _process_eltyp28(self, restype, pos, d):
        """SCQS, Subparametric Curved Quadrilateral Thick Shell (8-node)"""
        d.shape = (2, 4, 5)  # (side, respt, comp)
        if pos == 'respts':
            return np.rollaxis(d, 1)
        elif pos == 'nodes':
            d[...,:3] = self._extrapolate_to_surface(d[...,:3])
            if restype == 'generalstress':
                d[...,3:] = 0.
            elif restype == 'decomposedstress':
                d = self._decompose_stresses_thick(d)  # shape is (4, 8)
            # extrapolate to nodes:
            return np.dot(self.xquad, d)  # shape is (8, 2, 5)

    def _extrapolate_to_surface(self, d):
        lp, up = d[0], d[1]
        a = self.s3 / 3
        ls = lp + (-1 + a) / (2*a) * (up - lp)
        us = lp + (1 + a) / (2*a) * (up - lp)
        return np.array((ls, us))

    def _decompose_stresses_thick(self, d):
        return np.dot(self.x_decomp_thick, np.rollaxis(d,1).reshape(-1,10).T).T

    def _decompose_stresses_thin(self, d):
        return np.dot(self.x_decomp_thin, d.reshape(-1,6).T).T


    def get_nodes(self, sets=None, kind=None, disconnected=False, trans=None,
                  index=1):
        """Get node coordinates.

        Parameters
        ----------
        sets : str or sequence of str
            set name or sequence of set names. If sets=None (default), all
            nodes in the superelement is returned.
        kind : str
            'beam', 'shell' or None (None is default, returns all kinds)
        disconnected : bool
            if True, elements do not share nodes. Instead each element has its
            own set of nodes. This results in more nodes than in the original
            data.
        trans : numpy.ndarray, int or 'top'
            | four calling patterns are supported:
            | 1: provide a 4x4 transformation matrix
            | 2: provide the hierarchy level (int)
            | 3: provide 'top' to specify the top level
            | 4: None (default): No transformation is applied
        index : int
            superelement index number. Relevant when applying hierarchy
            transformations (T-files only). Default is 1.

        Returns
        -------
        coords : numpy.ndarray
            A 2d array with shape (nnodes, 3), where the second axis represent
            the x, y and z coordinates.
        """

        # index is redundant for resultsdata and non-hierarchy inputdata

        gcoord = self._get_record('gcoord')
        nodeindices = self._get_nodeindices(sets, kind, disconnected)
        recs = gcoord[:].take(nodeindices, axis=0)  # is allegedly faster

        # add a column of ones in order to perform affine transformation
        coords = np.column_stack((recs['xcoord'],
                                  recs['ycoord'],
                                  recs['zcoord'],
                                  np.ones(len(recs), 'f')))

        # transformation
        if trans is not None:
            if self._filetype == 'R':
                index = self._index
            if isinstance(trans, np.ndarray):
                if trans.shape != (4,4):
                    raise ValueError(
                        'Transformation matrix must have shape (4,4), got '
                        '({},{})'.format(*trans.shape))
            elif self._toplevel:
                if isinstance(trans, str):
                    if trans == 'top':
                        trans = self._toplevel.get_transformation(self._seltyp,
                                                                  index, 'top')
                    else:
                        raise ValueError('legal string argument: "top"')
                elif isinstance(trans, int):
                    trans = self._toplevel.get_transformation(self._seltyp,
                                                              index, trans)
            elif self._has_record('rsuptran') and trans=='top':
                trans = self._get_record('rsuptran')[0]['t']
            else:
                raise HierarchyError(
                    '{0}: No transformation data for trans={1}. {0} has no '
                    'toplevel reference.'.format(self.name, trans))
            coords[:] = np.dot(coords, trans)


        return coords[:,:3]

    def get_nodenumbers(self, sets=None, kind=None, disconnected=False,
                        numbertype='external'):
        """Get internal or external (default) node numbers.
        """
        gnode = self._get_record('gnode')
        if numbertype == 'external':
            col = 'nodex'
        elif numbertype == 'internal':
            col = 'nodeno'
        else:
             raise ValueError("numbertype must be 'internal' or 'external'")
        if any((sets, kind, disconnected)):
            nodeindices = self._get_nodeindices(sets, kind, disconnected)
            return gnode.col(col).take(nodeindices)
        else:
            return gnode.col(col)

    def get_noderesults(self, restype, run=1, rescases=None,
                        sets=None, kind=None, disconnected=False):
        """Get node result data.

        Parameters
        ----------
        restype : str
            'displacement', 'velocity' or 'acceleration'
        run : int
            Analysis run number (default is 1)
        rescases : int, sequence or None
            External result case number(s)
        sets : str or sequence of str
            set name or sequence of set names. If sets=None (default), all
            elements in the current dataset is returned.
        kind : str
            'beam', 'shell' or None (None is default, returns all kinds)
        disconnected : bool
            if True, elements do not share nodes. Instead each element has its
            own set of nodes. This results in more nodes than in the original
            data.

        Returns
        -------
        result : numpy.ndarray
            A 3d array with shape (nrescases, nnodes, 6), where the last axis
            represent the six degrees of freedom for the result type.
        """

        if restype == 'displacement':
            rvnodres, res = self._get_record('rvnoddis')
            start, stop = 'dis_start', 'dis_stop'
        elif restype == 'velocity':
            rvnodres, res = self._get_record('rvnodvel')
            start, stop = 'vel_start', 'vel_stop'
        elif restype == 'acceleration':
            rvnodres, res = self._get_record('rvnodacc')
            start, stop = 'acc_start', 'acc_stop'
        else:
            raise ValueError(
                "restype must be 'displacement', 'velocity' or 'acceleration'")

        resrefs = self._get_resrefs(run, rescases)
        complexflag = self._check_complex_and_get_flag(resrefs['icompl'])

        # get data
        nodenumbs = set(
            self.get_nodenumbers(sets, kind, numbertype='internal'))
        ires = set(resrefs['ires'])

        data = []
        for r in rvnodres.iterrows():
            if r['ires'] in ires and r['iinod'] in nodenumbs:
                data.append(res[r[start]:r[stop]])  # lookup in res is slow..
        data = np.concatenate(data)

        # change dtype if complex
        if complexflag:
            data.dtype = np.complex64

        # reshape
        ncomps = 6  # for now assume allways 6 comps for node results..
        nres = len(ires)  # TODO: need an axis for rescase!
        if nres > 1:
            data.shape = (nres, -1, ncomps)
        else:
            data.shape = (-1, ncomps)

        # disconnect
        if disconnected:
            indices = self._get_nodeindices(sets, kind, disconnected)
            indices = self._update_connectivity(indices)
            axis = 0 if nres == 1 else 1
            return data.take(indices, axis=axis)

        return data

        #   - get internal rescase number(s) from external and run no.
        #     e.g. self._get_ires(run, rescases)
        #   - check that results are either real or complex
        #   - get internal nodenumbers
        #   - iterate rvnodres:
        #       data = []
        #       if r[ires] in ires_set and r[iinod] in nodenumbers:
        #           data.append(res[r[start]:r[stop]])
        #       data = np.concatenate(data)
        #   - reshape according to number of components
        #   - if disconnected:
        #       indices = self._get_nodeindices(sets, kind, disconnected)
        #       data.take(indices, axis=0)


        # TODO:
        # sets, kind and disconnected are used everywhere, could it
        # optionally be set global defaults for these ?
        # fs.set_defaults(sets=[...], kind='shell', disconnected=True)
        # Then at the beginning in the functions who have these parameters:
        # sets, kind, disconnected = self._get_global_defaults(
        #   sets, kind, disconnected)
        # if a value is passed to the function (other than the local default)
        # this will override the global default and the argument will be
        # returned unchanged. Same if a global default have not
        # been set.


        # result types:
        #   displacement
        #   velocity
        #   acceleration
        #   reactions ?
        #   nodal average of element results (coplanar elements only)

        # rescase patterns:
        #   None (default) returns all
        #   rescase (int): choose a single resultcase
        #   sequence of rescases (ints): (1,4,5,6)

    def get_elements(self, sets=None, kind=None, disconnected=False):
        """Get element connectivity.

        Parameters
        ----------
        sets : str or sequence of str
            set name or sequence of set names. If sets=None (default), all
            elements in the current dataset is returned.
        disconnected : bool
            if True, elements do not share nodes. Instead each element has its
            own set of nodes. This results in more nodes than in the original
            data.
        kind : str
            'beam', 'shell' or None (None is default, returns all kinds)

        Returns
        -------
        Returns three 1d arrays:

        connectivity : numpy.ndarray
            Element definitions are given sequentially as indices into a
            corresponding *coords* array, ref *get_nodes* method. the
            *get_nodes* method must have been called with the same values for
            *sets*, *disconnected*, and *kind*.
        offset : numpy.ndarray
            indices into *connectivity* representing the end of an element
            definition. This array has length equal to the number of elements
            returned.
        eltyp : numpy.ndarray
            element type id. corresponds to the *offset* array. See table below
            for the supported element types.

        Supported element types
        -----------------------
        ==== ====== ============ ==============================================
        Name   Id   No. of nodes Description
        ==== ====== ============ ==============================================
        BEAS  15         2       3D 2 Node Beam
        BTSS  23         3       General Curved Beam
        FQUS  24         4       Flat Quadrilateral Thin Shell
        FTRS  25         3       Flat Triangular Thin Shell
        SCTS  26         6       Subparametric Curved Triangular Thick Shell
        SCQS  28         8       Subparametric Curved Quadrilateral Thick Shell
        ==== ====== ============ ==============================================

        """

        # TODO: should disconnected be True as default? If this turns out to be
        #       the most common use. Will this always be True when working with
        #       element results?

        gelmnt1, nodin = self._get_record('gelmnt1')

        if sets or kind:
            elemindices = self._get_elementindices(sets, kind)
            con, offset = self._get_connectivity(elemindices)
            eltyp = gelmnt1.col('eltyp')[elemindices]
            if disconnected:
                con = np.arange(len(con))
            else:
                con = self._update_connectivity(con)
        else:
            if disconnected:
                con = np.arange(len(nodin))
            else:
                con = nodin[:] - 1
            offset = gelmnt1.col('nodin_stop')
            eltyp = gelmnt1.col('eltyp')

        return con, offset, eltyp

    def get_elementnumbers(self, sets=None, kind=None, numbertype='external'):
        """Get internal or external (default) element numbers.
        """

        gelmnt1, nodin = self._get_record('gelmnt1')
        if numbertype == 'external':
            col = 'elnox'
        elif numbertype == 'internal':
            col = 'elno'
        else:
             raise ValueError("numbertype must be 'internal' or 'external'")
        if sets or kind:
            elemindices = self._get_elementindices(sets, kind)
            return gelmnt1.col(col).take(elemindices)
        else:
            return gelmnt1.col(col)

    def get_elementresults(self, restype, pos='nodes', run=1, rescases=None,
                           sets=None):
        """Get element result data.

        Parameters
        ----------
        restype : str
            'beamforce', 'generalstress' or 'decomposedstress'
        pos : str
            'nodes' (default), 'respts' (not supported yet) or average
            (not supported yet)
        run : int
            Analysis run number (default is 1)
        rescases : int, sequence or None
            External result case number(s)
        sets : str or sequence of str
            optionally limit the data by set(s)

        Returns
        -------
        result : numpy.ndarray
            The shape of the returned array will depend on the provided values
            for the *restype* and *pos* parameters. See below for description
            of the different types.

        Result type: Beam Force
        -----------------------
        ...

        ===== ===== ========== ==========================================
        Pos.  Elem.  No. of    No. of comps.
              Type   Res. pts.
        ===== ===== ========== ==========================================
        nodes BEAS  2          6 (nxx, nxy, nxz, mxx, mxy, mxz)
        nodes BTSS  3          6 (nxx, nxy, nxz, mxx, mxy, mxz)
        ===== ===== ========== ==========================================



        Result type: General Stress
        ---------------------------
        ...

        ===== ===== ========== ==========================================
        Pos.  Elem.  No. of    No. of comps.
              Type   Res. pts.
        ===== ===== ========== ==========================================
        nodes FQUS   4            3 (sigxx, sigyy, tauxy)
        nodes FTRS   3            3 (sigxx, sigyy, tauxy)
        nodes SCTS   6            5 (sigxx, sigyy, tauxy, tauxz, tauyz)
        nodes SCQS   8            5 (sigxx, sigyy, tauxy, tauxz, tauyz)
        ===== ===== ========== ==========================================

        """

        if restype == 'beamforce':
            rvres, res = self._get_record('rvforces')
            start, stop = 'force_start', 'force_stop'
            kind = 'beam'
        elif restype in ('generalstress', 'decomposedstress'):
            rvres, res = self._get_record('rvstress')
            start, stop = 'stress_start', 'stress_stop'
            kind = 'shell'
        else:
            raise ValueError(
                'restype={!r} not supported'.format(restype))

        resrefs = self._get_resrefs(run, rescases)
        complexflag = self._check_complex_and_get_flag(resrefs['icompl'])

        # get data
        elemnumbs = set(
            self.get_elementnumbers(sets, kind, numbertype='internal'))
        ires = set(resrefs['ires'])

        data = []
        for r in rvres.iterrows():
            if r['ires'] in ires and r['iielno'] in elemnumbs:

                eltyp = r['ir']
                d = res[r[start]:r[stop]]

                if complexflag:
                    d.dtype = np.complex64

                data.append(
                    self._process_element_results(eltyp, restype, pos, d))

        data = np.concatenate(data)  # shape is (nres*npos, comp)
                                     # or (nres*npos, surface, comp)

        # reshape
        if len(ires) > 1:
            if restype == 'generalstress':
                ncomps = data.shape[2]
                data.shape = (len(ires), -1, 2, ncomps)
            else:
                ncomps = data.shape[1]
                data.shape = (len(ires), -1, ncomps)

        return data

        # this function could also support 'kind':
        # The purpose would be (for kind=None (all)) to have same shape datasets
        # for both beam and shell results which then could be used with a single
        # dataset for nodes and elements. irrelevant element types would be
        # given zero or nan values.

        #   - get internal rescase number(s) from external and run no.
        #     e.g. self._get_ires(run, rescases)
        #   - check that results are either real or complex
        #   - get internal elementnumbers
        #   - get element types: eltyps = np.unique(rvstress.col('ir'))

        #   - iterate rvres:
        #   - Alt 1 separate elementtypes
        #       data = {}
        #       for eltyp in eltyps:
        #           data[eltyp] = []
        #       if r[ires] in ires_set and r[iinod] in nodenumbers:
        #           data[r[ir]].append(res[r[start]:r[stop]])

        #       for k,v in data.iteritems:
        #           data[k] = np.concatenate(v)

        #   - we now have data for each element type, which is needed to
        #     perform extrapolation etc..
        #   - but have we lost the original order??
        #   - data must be put back in original order berfore returning...
        #   - but reshaped first ?

        #   - Alt 2 perform manipulations inside loop
        #       data = []
        #       if r[ires] in ires_set and r[iinod] in nodenumbers:
        #
        #           d = res[r[start]:r[stop]]
        #
        #           perform processing:
        #           extrapolation etc..
        #
        #           data.append(d)

        #       data = np.concatenate(data)
        #       change dtype if complex
        #       reshape data according to no. of components
        #       must ensure that all element types have the same number of
        #       components...

        # Sequence of stress computations:
        #   1. stresses at result points (as on resultfile)
        #   2. surface stresses at result points (extrapolation for 2nd order)
        #   3. element (node upper/lower) stresses (extrapolation to nodes)
        #   4. element average (at upper/lower) by averaging elem corner vals
        #   5. nodal average
        #   6. decomposed stresses (bending/membrane)
        #   7. Combinations (complex results must be evaluated/expanded)
        #   8. principal and von mises stress calc


        # kind:
        #   beam
        #   shell
        #   solid

        # restypes:
        #   beamforce (beams)
        #   beamstress (beams)
        #   generalstress (shell and solid)
        #   decomposedstress (i.e. membrane/bending) (shell only)
        #   surfaceloads (shell and solid)

        # pos:
        #   node (default)
        #   respt
        #   average

        # layer: (consider returning both?)
        #   upper
        #   lower

        # assume all elements have same (number of) components

        # TODO
        # provide functions to operate on result data from get_element results?
        # freesif.calc.vonmises(arr, ...)
        # freesif.calc.principal(arr, ...)
        # freesif.calc.complex_eval(arr, phase)
        # freesif.calc.complex_expand(arr, phase_step=10)


    def get_concepts(self, sets=None):
        pass

    def get_conceptresults(self, sets=None):
        pass

    def get_conceptnames(self, sets=None):
        pass

    def get_section_properties(self, sets=None):
        """Get beam section properties
        """
        gelref1, geono, fixno, eccno, transno = self._get_record('gelref1')
        gbeamg = self._get_record('gbeamg')
        gbeamg_arr = gbeamg[:]
        elemindices = self._get_elementindices(sets, 'beam')

        geono_bm = gelref1.col('geono').take(elemindices)
        gbeamg_indices = dict(zip(gbeamg.col('geono'), range((len(gbeamg)))))
        return np.array([gbeamg_arr[gbeamg_indices[gbm]] for gbm in geono_bm])

    def get_shell_thicknesses(self, sets=None):
        """Get shell thicknesses
        """
        gelref1, geono, fixno, eccno, transno = self._get_record('gelref1')
        gelth = self._get_record('gelth')
        gelth_arr = gelth.col('th')
        elemindices = self._get_elementindices(sets, 'shell')

        geono_sh = gelref1.col('geono').take(elemindices)
        gelth_indices = dict(zip(gelth.col('geono'), range((len(gelth)))))
        return np.array([gelth_arr[gelth_indices[gsh]] for gsh in geono_sh])


    # to be memory efficient:
    # iter_noderesults(), iter_elementresults(), iter_conceptresults()  ??

    def get_properties(self, prop_type, sets=None, concepts=False):
        # choose per element or per concept
        # one structured array with properties, and another array of ints
        # corresponding to elements/concepts array pointing into the properties
        # array

        # only relevant for multi field properties like section, material etc?

        pass

    def get_resultcase_info(self, run=None, rescases=None):
        pass

    def get_result_component_names(self, restype):
        pass

    def get_resultcase_names(self, rescases=None):
        """Get loadcase names.

        Parameters
        ----------
        rescases : int, sequence or None
            External result case number(s)

        Returns
        -------
        result : dict
            Result case numbers as keys, result case names as values.
        """
        loadtable, _ = self._get_record('tdload')
        idnos = loadtable.col('idno')
        names = loadtable.col('name')
        if rescases is None:
            pass
        else:
            extract_elements = np.isin(idnos, rescases)
            idnos = idnos[extract_elements]
            names = names[extract_elements]
        return dict(zip(idnos, names))
