# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 22:25:38 2018

@author: audun
"""

import numpy as np


def absmax(a, axis=None):
    """Signed abs max
    """
    amax = a.max(axis)
    amin = a.min(axis)
    return np.where(-amin > amax, amin, amax)


def find_elemres_indices(pairs, elems, elemnumb, nodenumb):
    """Get indices into array axis representing element node result points
    for a list of external (elementnumber, nodenumber) pairs.

    This is typically used to get results at specific element nodes.
    """

    con, offset, eltyp = elems
    indices = []
    for elno, nodeno in pairs:
        nodeind_arr, = np.where(nodenumb==nodeno)
        elind, = np.where(elemnumb==elno)[0]
        tmp_arr = offset[elind] - nodeind_arr
        tmp_ind = tmp_arr[tmp_arr>0].argmin()
        indices.append(nodeind_arr[tmp_ind])
    return np.array(indices)
