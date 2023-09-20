# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""
"""

from .sequentialparser.sif2hdf5 import sif2hdf5
from .data.file import File, open_hdf5, open_sif
from . import utils
from . import calc

__version__ = '0.1.2'
__all__ = ['sif2hdf5', 'File', 'utils', 'open_hdf5', 'open_sif', 'calc']
