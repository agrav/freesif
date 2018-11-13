# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""
"""

import numpy as np
from os import path
import re

from ..exceptions import ConditionError

def getrow(table, condition):
    """Get first row in a tables.Table fulfilling condition
    """
    for r in table.where(condition):
        return r.fetch_all_fields()
    raise ConditionError('Table "{}" has no row fulfilling '
                         'condition'.format(table.name))


def getvldata(arr, start_indices, stop_indices):
    """
    """
    data = []
    for start, stop in zip(start_indices, stop_indices):
        data.append(arr[start:stop])
    return np.concatenate(data)

def parse_sifname(filename):
    """returns prefix, letter, selno, extension
    """
    _, fname = path.split(filename)

    name_pattern = '(?i)(?P<prefix>[\w-]*?)' + \
                   '(?P<letter>[rtgl])' + \
                   '(?P<selno>[1-9]\d*)' + \
                   '(H[1-9]\d*)?.' + \
                   '(?P<ext>(si[fun]|fem))'

    match = re.match(name_pattern, fname)

    if not match:
        raise ValueError('{} is not a valid SESAM file name.'.format(fname))

    d = match.groupdict()
    return d['prefix'], d['letter'], d['selno'], d['ext']
