# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""Convert SESAM Interface files into HDF5 format.
"""
from __future__ import print_function

from os import path
import tables as tb
import re
from . import parsers
from .sequentialfiles import FormattedFile, UnformattedFile


def sif2hdf5(sifname, hdf5name=None, append=False, usefileprefix=True,
            prefix='', only=None, in_memory=False):

    """Convert  data from sesam sequential file(s) into HDF5 format

    Parameters
    ==========

    sifname : str or sequence
        Name of the (top level) Sesam file to be converted. Lower level
        files are found automatically. A sequence of file names is also
        accepted.

    hdf5name : str
        Optional. Name of the HFD5 file for which data are to be written.
        Defaults to `sifname` with .h5 extension. If `sifname` is a sequence of
        filenames and `hdf5name` is provided, all SIF files are read into the
        same file.

    append : bool
        Optional.
        If True: writes data next to any existing data on an existing .h5 file.
        If False (default): creates a new file overwriting any existing data.

    usefileprefix : bool
        Optional. If True (default): filename prefix will be added to hdf5
        group name. If False: any filename prefix will not be added. Applies to
        top level groupnames only. prefixes will never be added to groups below
        the top level

    prefix : str
        Optionally add a prefix to the hdf5 group name. If `sifname` is a
        sequence of filenames, `prefix` must be a sequence of (unique)
        prefixes.

    only : sequence
        Optionally provide a list of record names. If the list is provided,
        only records with those names will be parsed. **Warning**: Use with
        care, some records are required in order to read other records
        correctly.

    in-memory : bool
        if this parameter is True a tables.File instance is returned with a
        pure in-memory representation. the *hdf5name* and *append* parameters
        are ignored in this case.

    General
    =======

    There are three types of SESAM sequential files:

    1. Input Interface File (T#.FEM, L#.FEM)
    2. Structural Results Interface File (R#.SIF, R#.SIU)
    3. Hydrodynamic Results Interface File (G1.SIF)

    Sequential files are either formatted of unformatted. Both types are
    supported. Super element hierarchies are supported. Super element
    hierarchies (multiple Input- or Result-files) are read into a single .h5
    file. In addition, several independent models/datasets may reside side by
    side in the same .h5 file (with append=True).

    How the different file types are treated
    ========================================

    Input Interface files (T#.FEM)
    ------------------------------

    Both intermediate and 1st level superelements will be stored directly under
    the toplevel group.


    Structural Results Interface File (R#.SIF, R#.SIU)
    --------------------------------------------------

    Intermediate level data is not present on Structural Results files, only
    top and 1st level. 1st level group names will include the index number to
    separate repeating superelements.

    Hydrodynamic Results Interface File (G1.SIF)
    --------------------------------------------

    All records are stored directly under a single (top level) group.

    Typical h5 file layout (super element hierarichy)
    -------------------------------------------------

    - root

      - toplevel1 (group)

        - hierarch (table)
        - hierarch_ihsref (array)
        - higher_level_data1
        - higher_level_data2
        - etc..
        - sub_element_1 (group) (higher level)

          - higher_level_data1
          - higher_level_data2
          - etc..

        - sub_element_2 (group) (1st level)

          - 1st_level_data
          - gnode (table)
          - gcoord (table)
          - gelmnt1 (table)
          - gelmnt1_nodin (array)
          - etc..

        - sub_element_n (higher level or 1st level)

      - toplevel2 (group)

        - etc..

    Examples
    ========

    Convert a single SIF file::

        import freesif as fs
        fs.sif2hdf5('R1.SIU')

    Other notes
    ===========

    The *File* class provides high level methods to operate on the data
    records on a HDF5 file (via the *SifData* sub-classes).
    """

    # handle a sequence of file names
    if not isinstance(sifname, str):
        if not isinstance(sifname, (list, tuple)):
            raise TypeError('sifname must be a file name (str) or sequence '
                             'of file names')
        if in_memory:
            raise TypeError('A sequence of filenames is currently not '
                            'supported for the in-memory option')

        if not prefix:
            prefix = ['' for _ in sifname]

        elif isinstance(prefix, str):
            prefix = [prefix for _ in sifname]

        if hdf5name:
            sif2hdf5(sifname[0], hdf5name, append, usefileprefix, prefix[0],
                     only)
            for name, prefix in zip(sifname[1:], prefix[1:]):
                sif2hdf5(name, hdf5name, True, usefileprefix, prefix, only)
        else:
            for name, prefix in zip(sifname, prefix):
                sif2hdf5(name, None, append, usefileprefix, prefix, only)
        return

    # check append
    if not hdf5name and append:
        raise ValueError('Cannot append unless a hdf5name is provided.')


    # determine filename info
    real_path = path.realpath(sifname)
    fdir, fname = path.split(real_path)

    name_pattern = '(?i)(?P<prefix>[\w-]*?)' + \
                   '(?P<letter>[rtgl])' + \
                   '[1-9]\d*(H[1-9]\d*)?.' + \
                   '(?P<ext>(si[fun]|fem))'

    match = re.match(name_pattern, fname)

    if not match:
        raise ValueError('{} is not a valid SESAM file name.'.format(fname))

    groupdict = match.groupdict()
    fileprefix = groupdict['prefix']
    letter = groupdict['letter']
    extension = groupdict['ext']

    if extension.upper() == 'SIN':
        raise ValueError('.SIN files not supported. Convert to .SIU')

    # reset in case several files are read in the same session
    parsers.reset()

    # filetype identifier (T, L, R or G)
    filetypeid = letter

    # set group name prefix
    if usefileprefix:
        parsers.set_groupname_prefix(prefix + fileprefix)
    else:
        parsers.set_groupname_prefix(prefix)

    # create hdf5 filename if not given
    if not hdf5name or in_memory:
        hdf5name = fname[:-3] + 'h5'

    # determine formatting
    with open(sifname, 'rb') as f:
        if f.read(1) == b'K':  # file is unformatted
            is_formatted = False
        else:
            is_formatted = True

    # open files
    if is_formatted:
        in_file = FormattedFile(sifname)
    else:
        in_file = UnformattedFile(sifname)

    if in_memory:
        out_file = tb.open_file(hdf5name, 'w', driver='H5FD_CORE',
                                driver_core_backing_store=0)
    else:
        out_file = tb.open_file(hdf5name, 'a' if append else 'w')

    # TODO: check if group name allready exist

    # set filetypeid
    parsers.set_filetypeid(filetypeid)

    # parse first file
    parsers.parse_file(in_file, out_file, only)
    in_file.close()

    # parse additional lower level files

    # if model file(FEM): intermediate level files exist
    #                     only one file per superelement (not each repetition)

    # if results file(SIF/SIU): only 1st level files exist
    #                           one file per superelement repetition (index)

    # need to look for valid file names in directory of top level file..
    # generate valid file names based on hierarchy data..

    if parsers.has_hierarchy():

        hierarch_table = parsers.get_hierarch_table()

        # generate candidate file names
        fnames = []
        ihref = []
        if extension in ('SIF', 'SIU'):  # structural results files
            for r in hierarch_table[1:]:
                if r['islevl'] == 1:  # 1st level
                    fnames.append(fileprefix +
                                  'R' + str(r['iselty']) +
                                  'H' + str(r['ihref']) +
                                  '.' + extension)
                    ihref.append(r['ihref'])

        else:  # structural input files
            for r in hierarch_table[1:]:
                if r['indsel'] == 1:  # skip additional repetitions
                    fnames.append(fileprefix +
                                  'T' + str(r['iselty']) +
                                  '.' + extension)
                    ihref.append(r['ihref'])

            # TODO: look for load files as well ? (L-files)
            # Should records from L-files be put into corresponding T-file
            # group, or be put in a separate L-group? Let the user choose ?

        # read files
        for i, fname in enumerate(fnames):
            fname_fullpath = path.join(fdir, fname)
            if path.isfile(fname_fullpath):
                if is_formatted:
                    in_file = FormattedFile(fname_fullpath)
                else:
                    in_file = UnformattedFile(fname_fullpath)

                parsers.reset(reset_toplevel=False)
                parsers.set_current_ihref(ihref[i])
                parsers.parse_file(in_file, out_file, only)
                in_file.close()

            else:
                print('file not found: ' + fname)

    if in_memory:
        # TODO: is it possible to change mode to 'r' here ?
        return out_file

    out_file.close()

    # return hdf5 file name (full path)
    return path.abspath(hdf5name)
