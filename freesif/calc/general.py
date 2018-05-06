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
