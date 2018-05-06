# -*- coding: utf-8 -*-
"""
Created on Wed May 02 21:08:48 2018

@author: audun
"""

import numpy as np


def vonmises2d(arr):
    """2d von Mises stress calculation.

    *arr* must be at least a (...,3) shaped array where the three first
    elements on the last axis are SIGXX, SIGYY and TAUXY.

    Returns a arr.shape[:-1] shape array with the calculated von Mises stress.
    """
    sxx = arr[...,0]
    syy = arr[...,1]
    txy = arr[...,2]
    return np.sqrt(sxx**2 + syy**2 - sxx*syy + 3*txy**2)


def vonmises_thickshell(arr):
    raise NotImplementedError
