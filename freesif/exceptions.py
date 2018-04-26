# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""
"""


class NoSuchRecordError(Exception):
    pass


class ConditionError(Exception):
    pass


class HierarchyError(Exception):
    pass


class ResultError(Exception):
    pass


class ClosedFileError(Exception):
    pass


class VtfError(Exception):
    pass


class VtkError(Exception):
    pass
