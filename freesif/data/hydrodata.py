# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""Defines class to operate on SESAM hydrodynamic data.
"""

from __future__ import division

import numpy as np
from .sifdata import SifData
from .helpers import getrow
from ..exceptions import ResultError, NoSuchRecordError


class HydroData(SifData):
    """Class to operate on SESAM hydrodynamic data.
    """

    # mulitbody and timedomain results not tested/implemented yet !!

    # resptypes:
    #
    # 1st order:
    #   W1EXFORC - implemented
    #   W1MOTION - implemented
    #   W1SFORCE - implemented
    #   W1MATRIX - implemented
    #   W1PANPRE - not documented!
    #
    # 2nd order:
    #   W2EXFDIF
    #   W2EXFSUM
    #   W2FLUDIF
    #   W2FLUSUM
    #   W2HDRIFT - implemented
    #   W2MDRIFT - implemented
    #
    # other:
    #   WFLUIDKN - implemented
    #   WGRESPON
    #   WSURFACE

    # ASSUMPTIONS:
    #  -- all result records are ordered by increasing iwres
    #  -- wdresref recs are ordered by direction, then angular frequencies
    #     i.e. result data can be reshaped to (ndirs, nfreqs)
    #  -- all results on one .SIF file are associated with the same set of
    #     directions and frequencies
    #  -- only one ibcond for each ibody (not necessary the case for time
    #     domain? other situations?)

    def __init__(self, tbgroup, fileinstance):
        """
        """
        super(HydroData, self).__init__(tbgroup, fileinstance)

        wdresref, dirs, freqs, times = self._get_record('wdresref')

        # get unique frequencies without sorting
        _, idx = np.unique(freqs[:,1], return_index=True)
        self._angfreqs = freqs[:,1][np.sort(idx)]
        self._periods = 2*np.pi / self._angfreqs

        # get unique directions without sorting
        _, idx = np.unique(dirs[:,1], return_index=True)
        self._directions = dirs[:,1][np.sort(idx)]

        self._timesteps = np.unique(times[:,1])

        self._nfreqs = len(self._angfreqs)
        self._ndirs = len(self._directions)
        self._ntimesteps = len(self._timesteps)

        try:
            self._nsections = len(self._get_record('wsection'))
        except NoSuchRecordError:
            self._nsections = 0

        try:
            self._npoints = len(self._get_record('wfkpoint'))
        except NoSuchRecordError:
            self._npoints = 0

        # set g, ro and depth as global instance variables
        wglobdef = self._get_record('wglobdef')[0]  # allways 1 record only
        self.g = wglobdef['g']
        self.ro = wglobdef['ro']
        self.depth = wglobdef['depth']

        self.domaintype = 'frequency' if self._nfreqs else 'time'

    def _get_ibcond(self, bodyid=1):
            # get internal body and condition reference number for bodyid
            # assume only one ibcond per bodyid (not the case if fwd speed?)
            wbodcon, refcond = self._get_record('wbodcon')
            return getrow(wbodcon, 'ibody=={}'.format(bodyid))['ibcond']

    def _get_results(self, rec_name, bodyid, *fields):
        """Get a set of columns (given by `fields`) from result records
        (`rec_name`). `bodyid` should be None if not relevant (e.g. wfluidkn)
        Returns an array for each field
        """

        res = self._get_record(rec_name)
        if bodyid:
            ibcond = self._get_ibcond(bodyid)
            # get all records for this ibcond
            res_recs = res.read_where('ibcond=={}'.format(ibcond))
        else:
            res_recs = res[:]
        # return res_recs[list(fields)]  # structured array
        if len(fields) > 1:
            return tuple(res_recs[field] for field in fields)
        return res_recs[fields[0]]

    def _assert_domaintype(self, domaintype):
        self._check_isopen()
        if not domaintype == self.domaintype:
            raise ResultError('No {} domain results on {}'.format(domaintype,
                                                                  self.name))

    def get_angular_freqs(self):
        """Get angular frequencies.
        """
        self._assert_domaintype('frequency')
        return self._angfreqs

    def get_periods(self):
        """Get periods.
        """
        self._assert_domaintype('frequency')
        return self._periods

    def get_directions(self, unit='degrees'):
        """Get wave directions in 'degrees' (default) or 'radians'
        """
        self._check_isopen()

        if unit == 'degrees':
            return self._directions * 180. / np.pi
        elif unit == 'radians':
            return self._directions
        else:
            raise ValueError("unit must be 'degrees' or 'radians'")

    def get_timesteps(self):
        """Get time instants.
        """
        self._assert_domaintype('time')
        return self._timesteps

    def get_sections(self):
        """Get section definitions

        Returns
        -------
        sections : numpy structured array
            a numpy structured array with fields *secusr*, *isecty*, *p1*, *p2*
            and *p3*. See the Sesam Interface File documentation for an
            explanation of these fields.
        """
        wsection = self._get_record('wsection')
        return wsection[:][['secusr','isecty','p1','p2','p3']]

    def get_points(self):
        """Get fluid kinematics reference points

        Returns
        -------
        points : numpy.ndarray
            A 2d array with shape (npts, 3), where the three positions along
            the second axis represent the x, y and z coordinates.
        """


        wfkpoint = self._get_record('wfkpoint')
        return wfkpoint[:]['fkpnt']

    def get_bodyproperties(self, bodyid=1):
        """Get body characteristica as a dict
        """
        wbody = self._get_record('wbody')
        ibcond = self._get_ibcond(bodyid)
        row = getrow(wbody, 'ibcond=={}'.format(ibcond))
        bodyprops = dict(zip(row.dtype.names, row))
        del bodyprops['nfield']
        del bodyprops['ibcond']
        return bodyprops

    def get_motion_raos(self, bodyid=1):

        """Get rigid body motion RAOs.

        Parameters
        ----------
        bodyid : int, optional
            External body identification number. Default is 1

        Returns
        -------
        data : 3d numpy.ndarray (complex64)
            motion data with shape (6, len(dirs), len(freqs)) where 1st
            axis corresponds to the motion degree of freedom, 2nd to direction
            and 3rd to frequency.
        """

        # should include an 'uncoupled' option ?

        self._assert_domaintype('frequency')
        data = self._get_results('w1motion', bodyid, 'data')

        # reshape data to (dof, direction, frequency)
        data.shape = (self._ndirs, self._nfreqs, 6)
        return np.rollaxis(data, 2)

    def get_excitationforce_raos(self, bodyid=1):
        """Get excitation force RAOs.

        Parameters
        ----------
        bodyid : int, optional
            External body identification number. Default is 1

        Returns
        -------
        data : 3d numpy.ndarray (complex64)
            force data with shape (6, len(dirs), len(freqs)) where 1st
            axis corresponds to the force degree of freedom, 2nd to direction
            and 3rd to frequency.
        """

        self._assert_domaintype('frequency')
        data = self._get_results('w1exforc', bodyid, 'data')

        # reshape data to (dof, direction, frequency)
        data.shape = (self._ndirs, self._nfreqs, 6)
        return np.rollaxis(data, 2)

    def get_sectionforce_raos(self, bodyid=1):
        """Get sectional force RAOs.

        Parameters
        ----------
        bodyid : int, optional
            External body identification number. Default is 1

        Returns
        -------
        data : 4d numpy.ndarray (complex64)
            force data with shape (len(sections), 6, len(dirs), len(freqs))
            where 1st axis corresponds to section, 2nd axis to the force degree
            of freedom, 3rd to direction and 4th to frequency.
        """

        self._assert_domaintype('frequency')
        data = self._get_results('w1sforce', bodyid, 'data')

        # reshape data to (section, dof, direction, frequency)
        data.shape = (self._nsections, self._ndirs, self._nfreqs, 6)
        return np.rollaxis(data, 3, 1)

    def get_fluidkinematics_raos(self, kind='elevation'):
        """Get fluid kinematics.

        Parameters
        ----------
        kind : 'elevation', 'pressure' or 'velocity', optional
            Used to specify result type. Default is 'elevation'
        Returns
        -------
        data : numpy.ndarray
            if kind = 'elevation' or 'pressure':
                shape = (len(points), len(dirs), len(freqs))

            if kind = 'velocity':
                shape = (len(points), 3, len(dirs), len(freqs)). The second
                dimension represent the x, y and z components of the paricle
                velocity
        """
        self._assert_domaintype('frequency')
        p, v = self._get_results('wfluidkn', None, 'p', 'v')
        if kind in ('elevation', 'pressure'):
            p.shape = (self._ndirs, self._nfreqs, self._npoints)
            p = np.rollaxis(p, 2)
            if kind == 'elevation':
                return p / self.ro / self.g
            return p
        elif kind == 'velocity':
            v.shape = (self._ndirs, self._nfreqs, self._npoints, 3)
            v = np.rollaxis(v, 3)
            return np.rollaxis(v, 3)
        else:
            raise ValueError("kind must be 'elevation', 'pressure' or "
                             "'velocity'")

    def get_matrix(self, mtyp, bodyid=1, bodyid2=1):
        """Get matrix.

        Parameters
        ----------
        mtyp : int
            Used to specify matrix type:
                - 11 = First order body mass matrix
                - 12 = First order added mass matrix (freq. dependant)
                - 21 = First order potential damping matrix (freq. dependant)
                - 22 = Viscous damping matrix
                - 31 = Hydrostatic restoring matrix
                - 32 = Mooring restoring matrix

        bodyid : int, optional
            External body identification number. Default is 1

        bodyid2 : int, optional
            External body identification number for coupling body. Default is 1

        Returns
        -------
        mat : numpy.ndarray
            | shape = (6, 6)
            | shape = (6, 6, len(freqs)) if frequency dependant
        """

        w1matrix = self._get_record('w1matrix')
        ibconi = self._get_ibcond(bodyid)
        ibconj = self._get_ibcond(bodyid2)
        cond = '(ibconi=={}) & (ibconj=={}) & (imtyp=={})'.format(ibconi,
                                                                  ibconj,
                                                                  mtyp)
        if mtyp in (11, 22, 31, 32):  # freq indep
            return getrow(w1matrix, cond)['hmat']
        elif mtyp in (12, 21):  # freq dep
            return np.rollaxis(w1matrix.read_where(cond)['hmat'], 0, 3)
        else:
            raise ValueError('mtyp must be 11, 12, 21, 22, 31 or 32, '
                             'got {}'.format(mtyp))

    def get_bodymass(self, bodyid=1, bodyid2=1):
        """Get first order body mass matrix
        """
        return self.get_matrix(11, bodyid, bodyid2)

    def get_addedmass(self, bodyid=1, bodyid2=1):
        """Get frequency dependant first order added mass matrices
        """
        return self.get_matrix(12, bodyid, bodyid2)

    def get_potentialdamping(self, bodyid=1, bodyid2=1):
        """Get frequency dependant first order potential damping matrices
        """
        return self.get_matrix(21, bodyid, bodyid2)

    def get_viscousdamping(self, bodyid=1, bodyid2=1):
        """Get viscous damping matrix, frequency independent
        """
        return self.get_matrix(22, bodyid, bodyid2)

    def get_hydrostatic_restoring(self, bodyid=1, bodyid2=1):
        """Get hydrostatic restoring matrix
        """
        return self.get_matrix(31, bodyid, bodyid2)

    def get_mooring_restoring(self, bodyid=1, bodyid2=1):
        """Get mooring restoring matrix
        """
        return self.get_matrix(32, bodyid, bodyid2)

    def get_meandrift(self, bodyid=1):
        """Get second order mean drift forces

        Parameters
        ----------
        bodyid : int, optional
            External body identification number. Default is 1

        Returns
        -------
        forces : 3d numpy.ndarray (complex64)
            force data with shape (6, len(dirs), len(freqs)) where 1st axis
            corresponds to force component, 2nd to direction and 3rd to
            frequency.
        """

        self._assert_domaintype('frequency')
        forces = self._get_results('w2mdrift', bodyid, 'smfor')
        forces.shape = (self._ndirs, self._nfreqs, 6)
        return np.rollaxis(forces, 2)

    def get_horiz_meandrift(self, bodyid=1):
        """Get second order horizontal mean drift forces

        Parameters
        ----------
        bodyid : int, optional
            External body identification number. Default is 1

        Returns
        -------
        forces : 3d numpy.ndarray (complex64)
            force data with shape (3, len(dirs), len(freqs)) where 1st axis
            corresponds to force component, 2nd to direction and 3rd to
            frequency.
        """

        self._assert_domaintype('frequency')
        forces = self._get_results('w2hdrift', bodyid, 'shfor')
        forces.shape = (self._ndirs, self._nfreqs, 3)
        return np.rollaxis(forces, 2)

    def _n_pressure_panels(self, bodyid=1):
        try:
            return len(np.unique(self._get_results('w1panpre', bodyid, 'ipan')))
        except NoSuchRecordError:
            return 0

    def _n_symmetry_parts(self, bodyid=1):
        try:
            return len(np.unique(self._get_results('w1panpre', bodyid, 'isymm1')))
        except NoSuchRecordError:
            return 0

    def get_pressure_panel_ids(self, bodyid=1):
        """Get panels exposed to panel pressure

        Parameters
        ----------
        bodyid : int, optional
            External body identification number. Default is 1

        Returns
        -------
        data : 3d numpy.ndarray (complex64)
            panel overview (panel number, symmetry) with shape (len(panels), len(symmetry_part))
            where 1st axis corresponds to the panel number,  2nd to the corresponding.
        """
        self._assert_domaintype('frequency')
        _npanels = self._n_pressure_panels(bodyid=bodyid)
        _nsymm = self._n_symmetry_parts(bodyid=bodyid)

        if _npanels != 0:

            panels = self._get_results('w1panpre', bodyid, 'ipan')
            symms = self._get_results('w1panpre', bodyid, 'isymm1')

            panel_ids = []
            for symm in np.unique(symms):
                for panel in np.unique(panels):
                    panel_ids.append((panel, symm))

            for panel, symm in zip(panels, symms):
                assert (panel, symm) in panel_ids

            return np.array(panel_ids)
        else:
            return None

    def get_panel_pressures(self, bodyid=1):
        """Get panel pressure

        Parameters
        ----------
        bodyid : int, optional
            External body identification number. Default is 1

        Returns
        -------
        data : 3d numpy.ndarray (complex64)
            force data with shape (len(panels), len(symmetry_part), len(dirs), len(freqs)) where 1st
            axis corresponds to the panel number,  2nd to the symmertry part, 3nd to direction
            and 4th to frequency.
        """

        self._assert_domaintype('frequency')
        _npanels = self._n_pressure_panels(bodyid=bodyid)
        _nsymm = self._n_symmetry_parts(bodyid=bodyid)

        if _npanels != 0:
            pressures = self._get_results('w1panpre', bodyid, 'p')
            pressures.shape = (self._ndirs, self._nfreqs,  _nsymm, _npanels)
            return np.moveaxis(pressures, [2, 3], [1, 0])
        else:
            return None
