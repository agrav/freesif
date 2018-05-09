# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 20:24:59 2018

@author: audun
"""

from __future__ import division

import numpy as np


def principal2d(arr):
    """2d principal stress calculation.

    *arr* must be at least a (...,3) shaped array where the three first
    elements on the last axis are SIGXX, SIGYY and TAUXY.

    Returns a (...,2) shape array with P1 and P2 on the last axis.
    """

    sxx = arr[...,0]
    syy = arr[...,1]
    txy = arr[...,2]
    part1 = (sxx + syy) / 2
    part2 = 0.5*np.sqrt((sxx - syy)**2 + 4*txy**2)
    p1 = part1 + part2
    p2 = part1 - part2
    return np.stack([p1,p2], axis=-1)


def principal_thickshell(arr):
    """Principal stress calculation for 6-node and 8-node shell elements.

    *arr* must be at least a (...,5) shaped array where the five first
    elements on the last axis are SIGXX, SIGYY, TAUXY, TAUXZ and TAUYZ.

    Returns a (...,2) shape array with P1 and P2 on the last axis.
    """

    # stress components
    sxx = arr[...,0]
    syy = arr[...,1]
    szz = np.zeros_like(syy)
    txy = arr[...,2]
    txz = arr[...,3]
    tyz = arr[...,4]

    # build stress tensor
    s = np.stack([sxx,txy,txz,txy,syy,tyz,txz,tyz,szz], axis=-1)
    s.shape = s.shape[:-1] + (3,3)

    # solve eigenproblem
    p,v = np.linalg.eig(s)

    # remove mid value
    arem = np.argmin(np.abs(p), axis=-1)
    p_mod = np.empty(p.shape[:-1] + (2,))
    for ind in np.ndindex(p.shape[:-1]):
        p_mod[ind] = np.delete(p[ind], arem[ind])

    # get P1 and P2
    p1 = p_mod.max(axis=-1)
    p2 = p_mod.min(axis=-1)

    return np.stack([p1,p2], axis=-1)











