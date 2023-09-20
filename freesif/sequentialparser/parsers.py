# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""'Low level' functions for reading/writing SESAM data records.
"""

# This module could probably be implemented more elegantly as a SifParser
# class. with all record defititions (readers, descriptors, atoms) imported
# from a separate module.

from __future__ import division

import numpy as np
import tables as tb


# ---Global variables-------------------------------------------------------- #
# some global variables are updated from reader functions when necessary

# size of continuation data types on an unformatted file
rec_size = 1024  # max allowable, will be set whenever an IRECSIZE is read
                 # currently not used anywhere..

# current superelement (should be reset for each in_file)
current_ihref = None  # need to be set from outside module when reading lower
                      # level data in a superelement hierarchy
                      # (by set_current_ihref())

current_indsel = None # used for structural results
                      # (should be reset for each in_file)

current_seltyp = None  # set by ident_reader

filetypeid = 'sel'  # typ set to T, L, R or G
                    # (should be reset for each top-level in_file)

# need to know if we are at top level in get_path
# (should be reset for each in_file)
attoplevel = False

# current root group
root_group = '/'

# group name prefix
groupname_prefix = ''

# hierarchy data exists (should be reset for each top-level in_file)
has_hierarch = False

# reference to pytables table of HIERARCH records (will be set from
# hierarch_reader, should be reset for each top-level in_file)
hierarch_table = None
hierarch_ihsref = None

# The toplevel HIERARCH record
# (should be reset for each top-level in_file)
toplevel_hrec = None

# The current HIERARCH record
# (should be reset for each in_file)
current_hrec = None

# size of temporary array used when reading records
tmp_arr_size = 100000

# ---Containers-------------------------------------------------------------- #

# Dict with info about element types
# elem_types[ELTYP] = (name, no_of_nodes, description)
# (by far complete)
elem_types = {15: ('BEAS', 2, '3D 2 Node Beam'),
              16: ('AXIS', 2, 'Axial Spring'),
              18: ('GSPR', 1, 'Spring to Ground'),
              23: ('BTSS', 3, 'General Curved Beam'),
              24: ('FQUS', 4, 'Flat Quadrilateral Thin Shell'),
              25: ('FTRS', 3, 'Flat Triangular Thin Shell'),
              26: ('SCTS', 6, 'Subparametric Curved'
                   'Triangular Thick Shell'),
              28: ('SCQS', 8, 'Subparametric Curved'
                   'Quadrilateral Thick Shell')}

# Dict to hold parser functions for each record type
readers = {}

# Dict to hold record descriptors
descriptors = {}

# Dict to hold atoms used for EArrays
atoms = {}

# Dict to hold tables refs (should be reset for each in_file)
tables = {}

# List to hold names of read records (should be reset for each in_file)
# added as an attribute ('RECORDS') to the tables group holding data from
# in_file
record_names = []

# Use this dict to provide info about the number of records for a given record
# type (to optimize pytables I/O speed).
record_count = {}

###
# Below is dicts used to temporary hold data read in one record, which is
# needed in another record. These data could also be looked up from the
# tables in the .h5 file, but this approach is believed to be simpler/faster.
###

# Dict to hold element type info
# element_data[elno] = eltyp
elem_data = {}

# Dict to (temporary) hold nlay, nok, noj, noi from the RDPOINTS record.
# populated from inside rdpoints_reader.
# rdpoints_data[icoref] = (nlay, nok, noj, noi)
# (should be reset for each in_file)
rdpoints_data = {}

# Dict to (temporary) hold nsp from the RDPOINTS record.
# populated from inside rdpoints_reader.
# nsp_data[(ispalt, iielno)] = nsp
# (should be reset for each in_file)
nsp_data = {}

# Dict ot hold complex flag (icompl) from RDRESREF records
# populated from inside rdresref_reader.
# rdresref_data[ires] = icompl
# (should be reset for each in_file)
rdresref_data = {}

# Dict ot hold lenrec from RDFORCES, RDNODRES, RDSTRAIN, RDSTRESS, RDELNFOR.
# lenrec_data[(rec_name, ir)] = lenrec
# ir = irforc for RDFORCES and so on..
# (should be reset for each in_file)
lenrec_data = {}


def reset(reset_toplevel=True):

    global rec_size, current_ihref, current_seltyp, tables, record_names, \
        rdpoints_data, nsp_data, rdresref_data, lenrec_data, has_hierarch, \
        hierarch_table, record_count, attoplevel, filetypeid, current_indsel, \
        current_hrec, hierarch_ihsref, toplevel_hrec

    rec_size = 1024
    current_ihref = None
    current_seltyp = None
    current_indsel = None
    current_hrec = None
    attoplevel = False
    tables = {}
    record_names = []
    rdpoints_data = {}
    nsp_data = {}
    rdresref_data = {}
    lenrec_data = {}
    record_count = {}
    if reset_toplevel:
        has_hierarch = False
        hierarch_table = None
        hierarch_ihsref = None
        toplevel_hrec = None
        filetypeid = 'sel'



# ---Helper classes---------------------------------------------------------- #

# classes to factor and encapsulate repeating code related to reading data into
# memory and flushing data to file.

class recdata(object):
    def __init__(self, table, *vldata):
        self.table = table
        self.vldata = vldata
        self.counter = 0
        self.recs = np.empty(tmp_arr_size, self.table.dtype)

    def add(self, rec):
        for vl in self.vldata:
            rec += (vl.start, vl.stop)
        self.recs[self.counter] = rec
        self.counter += 1
        if self.counter == tmp_arr_size:  # array is full
            self.table.append(self.recs)
            self.table.flush()
            for vl in self.vldata:
                vl.write()
            self.counter = 0

    def write(self):
        self.table.append(self.recs[:self.counter])
        self.table.flush()
        # write any connected vldata as well
        for vl in self.vldata:
            vl.write()


class vldata(object):
    def __init__(self, earr):
        self.arr = earr
        self.data = []
        self.startindex = len(self.arr)
        self.start = 0
        self.stop = 0

    def add(self, seq):
        if seq:
            self.start = self.startindex
            self.data += seq
            self.startindex += len(seq)
            self.stop = self.startindex
        else:
            self.start = 0
            self.stop = 0

    def write(self):
        if self.data:
            self.arr.append(np.array(self.data))
            self.data = []


# ---Helper functions-------------------------------------------------------- #


def has_hierarchy():
    return has_hierarch


def get_hierarch_table():
    return hierarch_table


def set_current_ihref(ihref):
    global current_ihref, current_hrec, current_indsel
    current_ihref = ihref
    if has_hierarch:
        current_hrec = hierarch_table[ihref-1]
        if filetypeid == 'R':
            current_indsel = current_hrec['indsel']


def set_groupname_prefix(prefix):
    global groupname_prefix
    groupname_prefix = prefix


def set_filetypeid(ftid):
    global filetypeid
    filetypeid = ftid


def get_path():

    if has_hierarch:  # hierarchy data exist for the data currently being read

        if current_indsel:  # 1st level structural results

            return groupname_prefix + filetypeid + \
                   str(hierarch_table[0]['iselty']) + \
                   '/' + filetypeid + str(current_seltyp) + \
                   '_' + str(current_indsel)

        elif attoplevel:
            return groupname_prefix + filetypeid + str(current_seltyp)

        else:  # structural model
            return groupname_prefix + filetypeid + \
                   str(hierarch_table[0]['iselty']) + \
                   '/' + filetypeid + str(current_seltyp)

    elif current_seltyp:  # no hierarchy data exist
        return groupname_prefix + filetypeid + str(current_seltyp)

    else:
        raise NameError('Cannot create a valid group name. Neither '
                        'superelement id (IDENT) or hierarchy data (HIERARCH)'
                        'has been read')

# function 'get_table' to provide correct table to store records in:
# if the record table is accompanied by EArray(s) this function
# should return a tuple (table, earr1, ... , earrn)


def get_table(rec_name, out_file):

    # check if table exist
    table = tables.get(rec_name, None)

    if table:
        return table

    else:
        where = root_group + get_path()
        descriptor = descriptors[rec_name]
        atom_data = atoms.get(rec_name, None)
        nrecs = record_count.get(rec_name, 10000)

        if rec_name == 'RDPOINTS':
            pass

        if atom_data:

            # subclass descriptor to add start and stop fields
            # (need to do this as the same descriptor may be used for several
            # record types)
            class newdescriptor(descriptor):
                pass

            vlarrs = []
            pos = len(descriptor.columns)
            for name, atom, shape in atom_data:

                # create EArray
                shape = (0,) + shape
                arrname = '{}_{}'.format(rec_name.lower(), name)
                vlarrs.append(out_file.create_earray(where, arrname, atom,
                                                     shape=shape,
                                                     createparents=True,
                                                     expectedrows=nrecs))

                # add start ond stop columns to table descriptor
                newdescriptor.columns[name + '_start'] = tb.Int32Col(pos=pos)
                pos += 1
                newdescriptor.columns[name + '_stop'] = tb.Int32Col(pos=pos)
                pos += 1

            table = out_file.create_table(where, rec_name.lower(),
                                          newdescriptor,
                                          createparents=True,
                                          expectedrows=nrecs)

            tables[rec_name] = (table,) + tuple(vlarrs)

        else:
            tables[rec_name] = out_file.create_table(where, rec_name.lower(),
                                                     descriptor,
                                                     createparents=True,
                                                     expectedrows=nrecs)

        return tables[rec_name]


def standard_record_reader(record_name, descriptor, skipfirst=0):
    """
    """
    ncols = len(descriptor.columns)
    n_add_fl = ncols - 4 if ncols > 4 else 0

    def record_reader(rec_name, rec, in_file, out_file):

        table = get_table(rec_name, out_file)
        recs = recdata(table)

        while rec_name == record_name:

            if n_add_fl:
                rec += in_file.read_floats(n_add_fl, skipfirst)

            if len(rec) != ncols:  # pad or remove values to match ncols
                if len(rec) < ncols:
                    rec += tuple(0. for x in range(ncols-len(rec)))
                else:
                    rec = rec[:ncols-len(rec)]

            recs.add(rec)

            rec_name, rec = in_file.read_headerrec()

        recs.write()
        return rec_name, rec

    return record_reader


def reader_optionalfields(record_name, descriptor, skipfirst=0, noptf=0):

    # records describing section type have two optional fields at the end,
    # 'nloby' and 'nlobz'. Need to check if present or not

    ncols = len(descriptor.columns)
    n_add_fl = ncols - 4 - noptf  # no. of fields after header rec and before
                                  # optional fields

    def record_reader(rec_name, rec, in_file, out_file):

        table = get_table(rec_name, out_file)
        recs = recdata(table)

        while rec_name == record_name:

            rec += in_file.read_floats(n_add_fl, skipfirst)

            if not in_file.at_headerrec():
                rec += in_file.read_floats(noptf, skipfirst)

            # pad or remove values to match ncols
            if len(rec) != ncols:
                if len(rec) < ncols:
                    rec += tuple(0. for x in range(ncols-len(rec)))
                else:
                    rec = rec[:ncols-len(rec)]

            recs.add(rec)

            rec_name, rec = in_file.read_headerrec()

        recs.write()
        return rec_name, rec

    return record_reader


# ---Parser ----------------------------------------------------------------- #


def general_parser(rec_name, rec, in_file, out_file):
    """Responsible for:
    - delegating to the appropriate record parser
    - handling not implemented records
    - maintaining a list of records read
    """
    # skip record if first float is negative:
    if rec[0] < 0.:
        if rec[0] == -5.:
            in_file.read_floatrec()
        return in_file.read_headerrec()

    # call record parser
    if rec_name in readers:
        record_names.append(rec_name)
        return readers[rec_name](rec_name, rec, in_file, out_file)

    # handle not implemented record
    nrecs = 0
    current_name = rec_name
    while rec_name == current_name:
        nrecs += 1
        rec_name, rec = in_file.skip_restofrecord(rec_name, rec)

#    print 'not implemented: {:<10} no. of records: {:>6}'.format(current_name,
#                                                                 nrecs)

    return rec_name, rec


def parse_file(in_file, out_file, only=None):

    rec_name, rec = in_file.read_headerrec()

    if only:
        while rec_name:
            if rec_name in only:
                rec_name, rec = general_parser(rec_name, rec, in_file,
                                                  out_file)
            else:
                skip_rec_name = rec_name
                while rec_name == skip_rec_name:
                    rec_name, rec = in_file.skip_restofrecord(rec_name, rec)
    else:
        while rec_name:
            rec_name, rec = general_parser(rec_name, rec, in_file, out_file)

    # set attributes on superelement group
    # ------------------------------------

    where = root_group + get_path()
    group = out_file.get_node(where)

    def remove_dups_preseve_order(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if x not in seen and not seen_add(x)]

    group._f_setattr('RECORDS', remove_dups_preseve_order(record_names))

    if filetypeid in ('T', 'L', 'R'):
        datattype = 'STRUCTURAL'
        if filetypeid == 'R':
            if attoplevel:
                if hierarch_table[0]['islevl'] == 1:
                    level = 'FIRSTLEVEL'
                    group._f_setattr('INDEX', current_indsel)
                else:
                    level = 'TOPLEVEL'

            else:
                level = 'FIRSTLEVEL'
                group._f_setattr('INDEX', current_indsel)

        elif filetypeid in ('T', 'L'):
            if attoplevel:
                level = 'TOPLEVEL'

            elif has_hierarch:  # intermediate level
                if current_hrec['islevl'] > 1:
                    level = 'INTERLEVEL'
                else:
                    level = 'FIRSTLEVEL'

            else:
                level = 'FIRSTLEVEL'



    elif filetypeid == 'G':
        datattype = 'HYDRODYNAMIC'
        level = ''
    else:
        datattype = 'UNKNOWN'
        level = ''

    group._f_setattr('LEVEL', level)
    group._f_setattr('DATATYPE', datattype)
    group._f_setattr('SELTYP', current_seltyp)
    group._f_setattr('FILETYPE', filetypeid)


### Record definitions ###


# ---IRECSIZE --------------------------------------------------------------- #

# record type is not written to hdf5 file
def irecsize(rec_name, rec, in_file, out_file):
    global rec_size
    rec_size = int(rec[1])
    return in_file.read_headerrec()

readers['IRECSIZE'] = irecsize


# ---HIERARCH---------------------------------------------------------------- #


# descriptor
class hierarch_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ihref = tb.Int32Col(pos=1)
    iselty = tb.Int32Col(pos=2)
    indsel = tb.Int32Col(pos=3)
    islevl = tb.Int32Col(pos=4)
    itref = tb.Int32Col(pos=5)
    ihpref = tb.Int32Col(pos=6)
    nsub = tb.Int32Col(pos=7)

descriptors['HIERARCH'] = hierarch_rec

atoms['HIERARCH'] = [('ihsref', tb.Int32Atom(), ())]

def hierarch_reader(rec_name, rec, in_file, out_file):

#    On T-files  (T#.FEM), all HIERARCH records are located on the highest
#    level (top-level) T-file only. Not on first or intermediate level files
#
#    On R-files (R#.SIF, R#.SIU) the HIERARCH records are located on the
#    top-level file if superelement analysis. If direct analysis (one super
#    element only), the single .SIU file includes one HIERARCH record.

    rec2 = in_file.read_floats(4)

    global current_seltyp
    current_seltyp = int(rec[2])

    table, ihsref_arr = get_table(rec_name, out_file)
    ihsrefs = vldata(ihsref_arr)
    recs = recdata(table, ihsrefs)

    nsub = int(rec2[3])
    ihsrefs.add(in_file.read_floats(nsub))
    recs.add(rec + rec2)

    rec_name, rec = in_file.read_headerrec()

    # read remaining HIERARCH records (if any)
    while rec_name == 'HIERARCH':
        rec2 = in_file.read_floats(4)
        nsub = int(rec2[3])
        ihsrefs.add(in_file.read_floats(nsub))
        recs.add(rec + rec2)

        rec_name, rec = in_file.read_headerrec()

    recs.write()

    # set module level swithches
    global has_hierarch, attoplevel
    has_hierarch = True
    attoplevel = True

    # keep a reference to HIERARCH table for subsequent lower level files
    global hierarch_table, hierarch_ihsref, toplevel_hrec, current_hrec
    hierarch_table = table
    hierarch_ihsref = ihsref_arr
    toplevel_hrec = current_hrec = hierarch_table[0]

    return rec_name, rec

readers['HIERARCH'] = hierarch_reader


# ---DATE, TEXT-------------------------------------------------------------- #


class text_rec(tb.IsDescription):
    type = tb.Int32Col(pos=0)
    subtype = tb.Int32Col(pos=1)
    nrecs = tb.Int32Col(pos=2)
    nbyte = tb.Int32Col(pos=3)

descriptors['DATE'] = text_rec
descriptors['TEXT'] = text_rec

atoms['DATE'] = [('recs', tb.StringAtom(72), ())]
atoms['TEXT'] = [('recs', tb.StringAtom(72), ())]


def text_reader(rec_name, rec, in_file, out_file):

    table, recs_arr = get_table(rec_name, out_file)
    textrecs = vldata(recs_arr)
    recs = recdata(table, textrecs)

    text_rec_name = rec_name
    while rec_name == text_rec_name:

        text = []
        for _ in range(int(rec[2])):  # nrecs
            text.append(in_file.read_stringrec())

        textrecs.add(text)
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

readers['DATE'] = text_reader
readers['TEXT'] = text_reader


# ---IDENT------------------------------------------------------------------- #


class ident_rec(tb.IsDescription):
    slevel = tb.Int32Col(pos=0)
    seltyp = tb.Int32Col(pos=1)
    selmod = tb.Int32Col(pos=2)

descriptors['IDENT'] = ident_rec

def ident_reader(rec_name, rec, in_file, out_file):

    global current_seltyp
    if not current_seltyp:
        current_seltyp = int(rec[1])

    table = get_table(rec_name, out_file)

    table.row['slevel'] = int(rec[0])
    table.row['seltyp'] = int(rec[1])
    table.row['selmod'] = int(rec[2])

    table.row.append()
    table.flush()

    return in_file.read_headerrec()

readers['IDENT'] = ident_reader


# ---IEND-------------------------------------------------------------------- #

# record type is not written to hdf5 file
def iend(rec_name, rec, in_file, out_file):
    if rec[0] == 0.:  # end of file
        return None, None
    else:
        return in_file.read_headerrec()

readers['IEND'] = iend


# ---TDxxxxx----------------------------------------------------------------- #


class td_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    idno = tb.Int32Col(pos=1)
    codnam = tb.Int32Col(pos=2)
    codtxt = tb.Int32Col(pos=3)
    name = tb.StringCol(64, pos=4)

rec_names = ('TDMATER', 'TDSECT', 'TDSETNAM', 'TDSUPNAM', 'TSLAYER',
             'TDSCONC', 'TDLOAD', 'TDRESREF', 'TDBODNAM', 'TDRSNAM',
             'TDSCATTER')

for rec_name in rec_names:
    descriptors[rec_name] = td_rec
    atoms[rec_name] = [('text', tb.StringAtom(64), ())]


def td_reader(rec_name, rec, in_file, out_file):

    table, text_arr = get_table(rec_name, out_file)
    textdata = vldata(text_arr)
    recs = recdata(table, textdata)

    td_rec_name = rec_name
    while rec_name == td_rec_name:

        codnam = int(rec[2])
        codtxt = int(rec[3])

        name = ''
        for _ in range(codnam // 100):
            add = in_file.read_stringrec().strip()
            if isinstance(add, str):
                name += add
            else:
                name += add.decode()

        text = []
        for _ in range(codtxt // 100):
            text.append(in_file.read_stringrec().strip())


        textdata.add(text)
        recs.add(rec + (name,))

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

for rec_name in rec_names:
    readers[rec_name] = td_reader


# ---GNODE------------------------------------------------------------------- #


# descriptor
class gnode_rec(tb.IsDescription):
    nodex = tb.Int32Col(pos=1)
    nodeno = tb.Int32Col(pos=2)
    ndof = tb.Int32Col(pos=3)
    odof = tb.Int32Col(pos=4)

descriptors['GNODE'] = gnode_rec

readers['GNODE'] = standard_record_reader('GNODE', gnode_rec)


# ---GCOORD------------------------------------------------------------------ #

# descriptor
class gcoord_rec(tb.IsDescription):
    nodeno = tb.Int32Col(pos=0)
    xcoord = tb.Float32Col(pos=1)
    ycoord = tb.Float32Col(pos=2)
    zcoord = tb.Float32Col(pos=3)

descriptors['GCOORD'] = gcoord_rec

readers['GCOORD'] = standard_record_reader('GCOORD', gcoord_rec)


# ---GELMNT1----------------------------------------------------------------- #


class gelmnt1_rec(tb.IsDescription):
    elnox = tb.Int32Col(pos=0)
    elno = tb.Int32Col(pos=1)
    eltyp = tb.Int32Col(pos=2)
    eltyad = tb.Int32Col(pos=3)

descriptors['GELMNT1'] = gelmnt1_rec

atoms['GELMNT1'] = [('nodin', tb.Int32Atom(), ())]


def gelmnt1_reader(rec_name, rec, in_file, out_file):
    """..."""

    table, nodin_arr = get_table(rec_name, out_file)
    nodins = vldata(nodin_arr)
    recs = recdata(table, nodins)

    while rec_name == 'GELMNT1':

        # store element type number for later use
        eltyp = elem_data[int(rec[1])] = int(rec[2])

        nnodes = elem_types[eltyp][1]
        nodins.add(in_file.read_floats(nnodes, skipfirst=8))

        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['GELMNT1'] = gelmnt1_reader


# ---GELMNT2----------------------------------------------------------------- #


class gelmnt2_rec(tb.IsDescription):
    subno = tb.Int32Col(pos=0)
    slevel = tb.Int32Col(pos=1)
    stype = tb.Int32Col(pos=2)
    addno = tb.Int32Col(pos=3)
    t = tb.Float32Col(pos=4, shape=(12,))
    nnod = tb.Int32Col(pos=5)

descriptors['GELMNT2'] = gelmnt2_rec

atoms['GELMNT2'] = [('nod', tb.Int32Atom(), ())]


def gelmnt2_reader(rec_name, rec, in_file, out_file):
    """..."""

    table, nod_arr = get_table(rec_name, out_file)
    nods = vldata(nod_arr)
    recs = recdata(table, nods)

    while rec_name == 'GELMNT2':

        t = in_file.read_floats(12, skipfirst=8)
        nnod = int(in_file.read_floats(1, skipfirst=8)[0])

        nods.add(in_file.read_floats(nnod, skipfirst=8))

        recs.add(rec + (t,) + (nnod,))

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['GELMNT2'] = gelmnt2_reader


# ---GELREF1----------------------------------------------------------------- #


class gelref1_rec(tb.IsDescription):
    elno = tb.Int32Col(pos=0)
    matno = tb.Int32Col(pos=1)
    addno = tb.Int32Col(pos=2)
    intno = tb.Int32Col(pos=3)
    mintno = tb.Int32Col(pos=4)
    strano = tb.Int32Col(pos=5)
    streno = tb.Int32Col(pos=6)
    strepono = tb.Int32Col(pos=7)
    geono = tb.Int32Col(pos=8)
    fixno = tb.Int32Col(pos=9)
    eccno = tb.Int32Col(pos=10)
    transno = tb.Int32Col(pos=11)

descriptors['GELREF1'] = gelref1_rec

atoms['GELREF1'] = [('geono', tb.Int32Atom(), ()),
                    ('fixno', tb.Int32Atom(), ()),
                    ('eccno', tb.Int32Atom(), ()),
                    ('transno', tb.Int32Atom(), ())]


def gelref1_reader(rec_name, rec, in_file, out_file):
    """..."""

    table, geono_arr, fixno_arr, eccno_arr, transno_arr = get_table(rec_name,
                                                                    out_file)

    geonos = vldata(geono_arr)
    fixnos = vldata(fixno_arr)
    eccnos = vldata(eccno_arr)
    transnos = vldata(transno_arr)
    recs = recdata(table, geonos, fixnos, eccnos, transnos)

    while rec_name == 'GELREF1':
        rec2 = in_file.read_floatrec(skipfirst=8)
        rec3 = in_file.read_floatrec(skipfirst=8)

        eltyp = elem_data[int(rec[0])]
        nnodes = elem_types[eltyp][1]

        # geono
        if rec3[0] == -1.:
            geonos.add(in_file.read_floats(nnodes, skipfirst=8))
        else:
            geonos.add(None)

        # fixno
        if rec3[1] == -1.:
            fixnos.add(in_file.read_floats(nnodes, skipfirst=8))
        else:
            fixnos.add(None)

        # eccno
        if rec3[2] == -1.:
            eccnos.add(in_file.read_floats(nnodes, skipfirst=8))
        else:
            eccnos.add(None)

        # transno
        if rec3[3] == -1.:
            transnos.add(in_file.read_floats(nnodes, skipfirst=8))
        else:
            transnos.add(None)

        recs.add(rec + rec2 + rec3)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

readers['GELREF1'] = gelref1_reader


# ---GECCEN------------------------------------------------------------------ #

class geccen_rec(tb.IsDescription):
    eccno = tb.Int32Col(pos=0)
    ex = tb.Float32Col(pos=1)
    ey = tb.Float32Col(pos=2)
    ez = tb.Float32Col(pos=3)

descriptors['GECCEN'] = geccen_rec
readers['GECCEN'] = standard_record_reader('GECCEN', geccen_rec)


# ---GUNIVEC----------------------------------------------------------------- #

class gunivec_rec(tb.IsDescription):
    transno = tb.Int32Col(pos=0)
    unix = tb.Float32Col(pos=1)
    uniy = tb.Float32Col(pos=2)
    uniz = tb.Float32Col(pos=3)

descriptors['GUNIVEC'] = gunivec_rec
readers['GUNIVEC'] = standard_record_reader('GUNIVEC', gunivec_rec)


# ---GELTH------------------------------------------------------------------- #

class gelth_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    th = tb.Float32Col(pos=1)
    nint = tb.Int32Col(pos=2, dflt=0)

descriptors['GELTH'] = gelth_rec
readers['GELTH'] = standard_record_reader('GELTH', gelth_rec)


# ---GBEAMG------------------------------------------------------------------ #

class gbeamg_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    void = tb.Float32Col(pos=1)
    area = tb.Float32Col(pos=2)
    ix = tb.Float32Col(pos=3)
    iy = tb.Float32Col(pos=4)
    iz = tb.Float32Col(pos=5)
    iyz = tb.Float32Col(pos=6)
    wxmin = tb.Float32Col(pos=7)
    wymin = tb.Float32Col(pos=8)
    wzmin = tb.Float32Col(pos=9)
    shary = tb.Float32Col(pos=10)
    sharz = tb.Float32Col(pos=11)
    shceny = tb.Float32Col(pos=12)
    shcenz = tb.Float32Col(pos=13)
    sy = tb.Float32Col(pos=14)
    sz = tb.Float32Col(pos=15)

descriptors['GBEAMG'] = gbeamg_rec
readers['GBEAMG'] = standard_record_reader('GBEAMG', gbeamg_rec, skipfirst=8)


# ---GBARM------------------------------------------------------------------- #

class GBARM_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    bt = tb.Float32Col(pos=2)
    bb = tb.Float32Col(pos=3)
    sfy = tb.Float32Col(pos=4)
    sfz = tb.Float32Col(pos=5)
    nloby = tb.Int32Col(pos=6)
    nlobz = tb.Int32Col(pos=7)

descriptors['GBARM'] = GBARM_rec
readers['GBARM'] = reader_optionalfields('GBARM', GBARM_rec, 8, 2)


# ---GBOX-------------------------------------------------------------------- #

class GBOX_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    t_b = tb.Float32Col(pos=3)
    tt = tb.Float32Col(pos=4)
    by = tb.Float32Col(pos=5)
    sfy = tb.Int32Col(pos=6)
    sfz = tb.Int32Col(pos=7)
    nloby = tb.Int32Col(pos=8)
    nlobz = tb.Int32Col(pos=9)

descriptors['GBOX'] = GBOX_rec
readers['GBOX'] = reader_optionalfields('GBOX', GBOX_rec, 8, 2)


# ---GCHAN------------------------------------------------------------------- #

class GCHAN_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    by = tb.Float32Col(pos=3)
    tz = tb.Float32Col(pos=4)
    sfy = tb.Int32Col(pos=5)
    sfz = tb.Int32Col(pos=6)
    empty = tb.Int32Col(pos=7)
    k = tb.Int32Col(pos=8)
    nloby = tb.Int32Col(pos=9)
    nlobz = tb.Int32Col(pos=10)

descriptors['GCHAN'] = GCHAN_rec
readers['GCHAN'] = reader_optionalfields('GCHAN', GCHAN_rec, 8, 2)


# ---GCHANR------------------------------------------------------------------ #

class GCHANR_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    by = tb.Float32Col(pos=3)
    tz = tb.Float32Col(pos=4)
    sfy = tb.Int32Col(pos=5)
    sfz = tb.Int32Col(pos=6)
    empty = tb.Int32Col(pos=7)
    k = tb.Int32Col(pos=8)
    r = tb.Int32Col(pos=9)
    nloby = tb.Int32Col(pos=10)
    nlobz = tb.Int32Col(pos=11)

descriptors['GCHANR'] = GCHANR_rec
readers['GCHANR'] = reader_optionalfields('GCHANR', GCHANR_rec, 8, 2)


# ---GDOBO------------------------------------------------------------------- #

class GDOBO_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    by = tb.Float32Col(pos=3)
    tt = tb.Float32Col(pos=4)
    t_b = tb.Float32Col(pos=5)
    sfy = tb.Int32Col(pos=6)
    sfz = tb.Int32Col(pos=7)
    nloby = tb.Int32Col(pos=8)
    nlobz = tb.Int32Col(pos=9)

descriptors['GDOBO'] = GDOBO_rec
readers['GDOBO'] = reader_optionalfields('GDOBO', GDOBO_rec, 8, 2)


# ---GIORH------------------------------------------------------------------- #

class GIORH_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    bt = tb.Float32Col(pos=3)
    tt = tb.Float32Col(pos=4)
    bb = tb.Float32Col(pos=5)
    t_b = tb.Float32Col(pos=6)
    sfy = tb.Int32Col(pos=7)
    sfz = tb.Int32Col(pos=8)
    nlobyt = tb.Int32Col(pos=9)
    nlobyb = tb.Int32Col(pos=10)
    nlobz = tb.Int32Col(pos=11)

descriptors['GIORH'] = GIORH_rec
readers['GIORH'] = reader_optionalfields('GIORH', GIORH_rec, 8, 3)


# ---GIORHR------------------------------------------------------------------ #

class GIORHR_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    bt = tb.Float32Col(pos=3)
    tt = tb.Float32Col(pos=4)
    bb = tb.Float32Col(pos=5)
    t_b = tb.Float32Col(pos=6)
    sfy = tb.Int32Col(pos=7)
    sfz = tb.Int32Col(pos=8)
    rt = tb.Float32Col(pos=9)
    rb = tb.Float32Col(pos=10)
    nlobyt = tb.Int32Col(pos=11)
    nlobyb = tb.Int32Col(pos=12)
    nlobz = tb.Int32Col(pos=13)

descriptors['GIORHR'] = GIORHR_rec
readers['GIORHR'] = reader_optionalfields('GIORHR', GIORHR_rec, 8, 3)


# ---GLSEC------------------------------------------------------------------- #

class GLSEC_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    by = tb.Float32Col(pos=3)
    tz = tb.Float32Col(pos=4)
    sfy = tb.Int32Col(pos=5)
    sfz = tb.Int32Col(pos=6)
    k = tb.Float32Col(pos=7)
    nloby = tb.Int32Col(pos=8)
    nlobz = tb.Int32Col(pos=9)

descriptors['GLSEC'] = GLSEC_rec
readers['GLSEC'] = reader_optionalfields('GLSEC', GLSEC_rec, 8, 2)


# ---GLSECR------------------------------------------------------------------ #

class GLSECR_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    by = tb.Float32Col(pos=3)
    tz = tb.Float32Col(pos=4)
    sfy = tb.Int32Col(pos=5)
    sfz = tb.Int32Col(pos=6)
    k = tb.Float32Col(pos=7)
    r = tb.Float32Col(pos=8)
    nloby = tb.Int32Col(pos=9)
    nlobz = tb.Int32Col(pos=10)

descriptors['GLSECR'] = GLSECR_rec
readers['GLSECR'] = reader_optionalfields('GLSECR', GLSECR_rec, 8, 2)


# ---GPIPE------------------------------------------------------------------- #

class GPIPE_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    di = tb.Float32Col(pos=1)
    dy = tb.Float32Col(pos=2)
    t = tb.Float32Col(pos=3)
    sfy = tb.Int32Col(pos=4)
    sfz = tb.Int32Col(pos=5)
    ncir = tb.Int32Col(pos=6)
    nrad = tb.Int32Col(pos=7)

descriptors['GPIPE'] = GPIPE_rec
readers['GPIPE'] = reader_optionalfields('GPIPE', GPIPE_rec, 8, 2)


# ---GUSYI------------------------------------------------------------------- #

class GUSYI_rec(tb.IsDescription):
    geono = tb.Int32Col(pos=0)
    hz = tb.Float32Col(pos=1)
    ty = tb.Float32Col(pos=2)
    bt = tb.Float32Col(pos=3)
    b1 = tb.Float32Col(pos=4)
    tt = tb.Float32Col(pos=5)
    bb = tb.Float32Col(pos=6)
    b2 = tb.Float32Col(pos=7)
    t_b = tb.Float32Col(pos=8)
    sfy = tb.Int32Col(pos=9)
    sfz = tb.Int32Col(pos=10)
    nlobyt = tb.Int32Col(pos=11)
    nlobyb = tb.Int32Col(pos=12)
    nlobz = tb.Int32Col(pos=13)

descriptors['GUSYI'] = GUSYI_rec
readers['GUSYI'] = reader_optionalfields('GUSYI', GUSYI_rec, 8, 2)


# ---GSETMEMB---------------------------------------------------------------- #


class gsetmemb_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    isref = tb.Int32Col(pos=1)
    index = tb.Int32Col(pos=2)
    istype = tb.Int32Col(pos=3)
    isorig = tb.Int32Col(pos=4)

descriptors['GSETMEMB'] = gsetmemb_rec

atoms['GSETMEMB'] = [('irmemb', tb.Int32Atom(), ())]


def gsetmemb_reader(rec_name, rec, in_file, out_file):
    """..."""

    table, irmemb_arr = get_table(rec_name, out_file)
    irmembs = vldata(irmemb_arr)
    recs = recdata(table, irmembs)

    while rec_name == 'GSETMEMB':

        rec += in_file.read_floats(1)
        n_irmemb = int(rec[0]) - 5

        irmembs.add(in_file.read_floats(n_irmemb))
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

readers['GSETMEMB'] = gsetmemb_reader


# ---MISOSEL----------------------------------------------------------------- #


# descriptor
class misosel_rec(tb.IsDescription):
    matno = tb.Int32Col(pos=0)
    young = tb.Float32Col(pos=1)
    poiss = tb.Float32Col(pos=2)
    rho = tb.Float32Col(pos=3)
    damp = tb.Float32Col(pos=4)
    alpha = tb.Float32Col(pos=5)

descriptors['MISOSEL'] = misosel_rec

readers['MISOSEL'] = standard_record_reader('MISOSEL', misosel_rec,
                                            skipfirst=8)

# ---SCONCEPT---------------------------------------------------------------- #


class sconcept_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ircon = tb.Int32Col(pos=1)
    scontype = tb.Int32Col(pos=2)
    sconrole = tb.Int32Col(pos=3)
    irparent = tb.Int32Col(pos=4)
    npart = tb.Int32Col(pos=5)
    njoin = tb.Int32Col(pos=6)

descriptors['SCONCEPT'] = sconcept_rec

atoms['SCONCEPT'] = [('irpart', tb.Int32Atom(), ()),
                     ('irjoin', tb.Int32Atom(), ())]


def sconcept_reader(rec_name, rec, in_file, out_file):
    """..."""

    table, irpart_arr, irjoin_arr = get_table(rec_name, out_file)
    irparts = vldata(irpart_arr)
    irjoins = vldata(irjoin_arr)
    recs = recdata(table, irparts, irjoins)

    while rec_name == 'SCONCEPT':

        scontype = int(rec[2])
        #sconrole = int(rec[3])

        rec += in_file.read_floats(1)

        # the SCONCEPT card will have different data depending on
        # type (and role?)

        if scontype == 7:  # segmented beam
            rec2 = in_file.read_floats(2)
            npart = int(rec2[0])
            njoin = int(rec2[1])

            rec += rec2

            # read part references
            irparts.add(in_file.read_floats(npart))
            irjoins.add(in_file.read_floats(njoin))

        elif scontype in [2]:  # segment
            rec += (0., 0.)
            irparts.add(None)
            irjoins.add(None)

        else:
            raise NotImplementedError('SCONCEPT: concept type {} not '
                                      'implemented'.format(scontype))

        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

readers['SCONCEPT'] = sconcept_reader


# ---SCONMESH---------------------------------------------------------------- #


class sconmesh_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ircon = tb.Int32Col(pos=1)
    numrep = tb.Int32Col(pos=2)
    typrep = tb.Int32Col(pos=3)
    nferep = tb.Int32Col(pos=4)


descriptors['SCONMESH'] = sconmesh_rec

atoms['SCONMESH'] = [('irferep', tb.Int32Atom(), ())]


def sconmesh_reader(rec_name, rec, in_file, out_file):
    """..."""

    # assume numrep is allways 1. Throw NotImplementedError if above 1
    if int(rec[2]) > 1:
        raise NotImplementedError('SCONMESH: NUMREP > 1 not implemented')

    table, irferep_arr = get_table(rec_name, out_file)
    irfereps = vldata(irferep_arr)
    recs = recdata(table, irfereps)

    while rec_name == 'SCONMESH':

        rec += in_file.read_floats(1)

        nferep = int(rec[4])
        irfereps.add(in_file.read_floats(nferep))
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['SCONMESH'] = sconmesh_reader


# ---BNBCD------------------------------------------------------------------- #


class bnbcd_rec(tb.IsDescription):
    nodeno = tb.Int32Col(pos=0)
    ndof = tb.Int32Col(pos=1)

descriptors['BNBCD'] = bnbcd_rec

atoms['BNBCD'] = [('fix', tb.Int32Atom(), ())]


def bnbcd_reader(rec_name, rec, in_file, out_file):

    table, fix_arr = get_table(rec_name, out_file)
    fixs = vldata(fix_arr)
    recs = recdata(table, fixs)

    while rec_name == 'BNBCD':

        ndof = int(rec[1])

        if ndof < 3:
            fixs.add(rec[2:2+ndof])
        else:
            fixs.add(rec[2:] + in_file.read_floats(ndof-2, skipfirst=8))

        recs.add(rec[:2])

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['BNBCD'] = bnbcd_reader


# ---BLDEP------------------------------------------------------------------- #

class bldep_rec(tb.IsDescription):
    nodeno = tb.Int32Col(pos=0)
    cnod = tb.Int32Col(pos=1)
    ndof = tb.Int32Col(pos=2)
    ndep = tb.Int32Col(pos=3)

descriptors['BLDEP'] = bldep_rec

atoms['BLDEP'] = [('deps', tb.Float32Atom(), (3,))]

# could actually use a table for these kind of variable length data..
# instead of arrays. Allows for (int, int, float) records

def bldep_reader(rec_name, rec, in_file, out_file):

    table, deps_arr = get_table(rec_name, out_file)
    deps = vldata(deps_arr)
    recs = recdata(table, deps)

    while rec_name == 'BLDEP':

        ndep = int(rec[3])
        deps.add([in_file.read_floatrec(skipfirst=8)[:3] for _ in range(ndep)])
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['BLDEP'] = bldep_reader


# ---BEUSLO------------------------------------------------------------------ #


class beuslo_rec(tb.IsDescription):
    llc = tb.Int32Col(pos=0)
    lotyp = tb.Int32Col(pos=1)
    complx = tb.Int32Col(pos=2)
    layer = tb.Int32Col(pos=3)
    elno = tb.Int32Col(pos=4)
    ndof = tb.Int32Col(pos=5)
    intno = tb.Int32Col(pos=6)
    side = tb.Int32Col(pos=7)

descriptors['BEUSLO'] = beuslo_rec

atoms['BEUSLO'] = [('rload', tb.Float32Atom(), ()),
                   ('iload', tb.Float32Atom(), ())]


def beuslo_reader(rec_name, rec, in_file, out_file):

    table, rload_arr, iload_arr = get_table(rec_name, out_file)
    rloads = vldata(rload_arr)
    iloads = vldata(iload_arr)
    recs = recdata(table, rloads, iloads)

    while rec_name == 'BEUSLO':

        rec += in_file.read_floatrec(skipfirst=8)
        ndof = int(rec[5])

        rloads.add(in_file.read_floats(ndof, skipfirst=8))

        if rec[2]:  # complex
            iloads.add(in_file.read_floats(ndof, skipfirst=8))
        else:
            iloads.add(None)

        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

readers['BEUSLO'] = beuslo_reader


# ---BGRAV------------------------------------------------------------------- #


class BGRAV_rec(tb.IsDescription):
    llc = tb.Int32Col(pos=0)
    empty1 = tb.Float32Col(pos=1)
    empty2 = tb.Float32Col(pos=2)
    opt = tb.Int32Col(pos=3)
    gx = tb.Float32Col(pos=4)
    gy = tb.Float32Col(pos=5)
    gz = tb.Float32Col(pos=6)

descriptors['BGRAV'] = BGRAV_rec

readers['BGRAV'] = standard_record_reader('BGRAV', BGRAV_rec,
                                          skipfirst=8)


# ---BRIGAC------------------------------------------------------------------ #


class BRIGAC_rec(tb.IsDescription):
    llc = tb.Int32Col(pos=0)
    empty1 = tb.Float32Col(pos=1)
    complx = tb.Int32Col(pos=2)
    empty2 = tb.Float32Col(pos=3)
    xcoord = tb.Float32Col(pos=4)
    ycoord = tb.Float32Col(pos=5)
    zcoord = tb.Float32Col(pos=6)
    empty3 = tb.Float32Col(pos=7)
    raccl1 = tb.Float32Col(pos=7)
    raccl2 = tb.Float32Col(pos=8)
    raccl3 = tb.Float32Col(pos=9)
    raccl4 = tb.Float32Col(pos=10)
    raccl5 = tb.Float32Col(pos=11)
    raccl6 = tb.Float32Col(pos=12)
    iaccl1 = tb.Float32Col(pos=13)
    iaccl2 = tb.Float32Col(pos=14)
    iaccl3 = tb.Float32Col(pos=15)
    iaccl4 = tb.Float32Col(pos=16)
    iaccl5 = tb.Float32Col(pos=17)
    iaccl6 = tb.Float32Col(pos=18)

descriptors['BRIGAC'] = BRIGAC_rec


def brigac_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)

    while rec_name == 'BRIGAC':

        rec += in_file.read_floats(10, skipfirst=8)

        if rec[2]:  # complex
            rec += in_file.read_floats(6, skipfirst=8)
        else:
            rec += (0., 0., 0., 0., 0., 0.)

        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

readers['BRIGAC'] = brigac_reader


# ---HSUPTRAN, RSUPTRAN------------------------------------------------------ #


# descriptor
class XSUPTRAN_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    itref = tb.Int32Col(pos=1)
    t = tb.Float32Col(pos=2, shape=(4, 4))

descriptors['RSUPTRAN'] = XSUPTRAN_rec
descriptors['HSUPTRAN'] = XSUPTRAN_rec


def xsuptran_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)

    suptran_name = rec_name
    while rec_name == suptran_name:

        table.row['nfield'] = int(rec[0])
        table.row['itref'] = int(rec[1])

        t = rec[2:] + in_file.read_floats(14)

        table.row['t'] = (t[:4], t[4:8], t[8:12], t[12:])

        table.row.append()
        rec_name, rec = in_file.read_headerrec()

    table.flush()
    return rec_name, rec


readers['RSUPTRAN'] = xsuptran_reader
readers['HSUPTRAN'] = xsuptran_reader


# ---RDFORCES, RDNODRES, RDSTRAIN, RDSTRESS, RDELNFOR------------------------ #

rec_names = ('RDFORCES', 'RDNODRES', 'RDSTRAIN', 'RDSTRESS', 'RDELNFOR')

class compdef_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ir = tb.Int32Col(pos=1)
    lenrec = tb.Int32Col(pos=2)


def compdef_reader(rec_name, rec, in_file, out_file):

    table, icomp_arr = get_table(rec_name, out_file)
    icomps = vldata(icomp_arr)
    recs = recdata(table, icomps)

    compdef_rec_name = rec_name
    while rec_name == compdef_rec_name:

        lenrec = int(rec[2])
        ir = int(rec[1])

        # store lenrec for later use
        lenrec_data[('RV' + rec_name[2:], ir)] = lenrec

        icomps.add(rec[3:] + in_file.read_floats(lenrec-1))
        recs.add(rec[:3])

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

for rec_name in rec_names:
    descriptors[rec_name] = compdef_rec
    atoms[rec_name] = [('icomp', tb.Int32Atom(), ())]
    readers[rec_name] = compdef_reader


# ---RDPOINTS---------------------------------------------------------------- #


# descriptor
class rdpoints_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ispalt = tb.Int32Col(pos=1)
    iielno = tb.Int32Col(pos=2)
    icoref = tb.Int32Col(pos=3)
    ieltyp = tb.Int32Col(pos=4)
    nsp = tb.Int32Col(pos=5)
    ijkdim = tb.Int32Col(pos=6)
    nsptra = tb.Int32Col(pos=7)
    nlay = tb.Int32Col(pos=8)

descriptors['RDPOINTS'] = rdpoints_rec

atoms['RDPOINTS'] = [('respoint', tb.Float32Atom(), (4,)),
                     ('xform', tb.Float32Atom(), (3, 3))]


def rdpoints_reader(rec_name, rec, in_file, out_file):

    table, respoints_arr, xforms_arr = get_table(rec_name, out_file)
    respoints = vldata(respoints_arr)
    xforms = vldata(xforms_arr)
    recs = recdata(table, respoints, xforms)

    while rec_name == 'RDPOINTS':

        rec += in_file.read_floats(5)

        ispalt = int(rec[1])
        iielno = int(rec[2])
        icoref = int(rec[3])
        nsp = int(rec[5])
        ijkdim = int(rec[6])
        nsptra = int(rec[7])
        nlay = int(rec[8])

        nok = ijkdim // 10000
        noj = np.mod(ijkdim, 10000) // 100
        noi = np.mod(ijkdim, 100)

        # store nlay, nok, noj, noi for later use
        rdpoints_data[icoref] = (nlay, nok, noj, noi)

        # store nsp for later use
        nsp_data[(ispalt, iielno)] = nsp

        # respoints
        rpts = []
        for k in range(nok):
            for j in range(noj):
                for i in range(noi):
                    ipoint = in_file.read_floats(1)
                    if ipoint[0] > 0.:
                        rpts.append(ipoint + in_file.read_floats(3))
#                    else:
#                        respoints2.append((-1., 0., 0., 0.))


        respoints.add(rpts)

        # xforms
        xfs = []
        for _ in range(nsptra):
            t = in_file.read_floats(9)
            xfs.append((t[:3], t[3:6], t[6:]))

        xforms.add(xfs)

        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['RDPOINTS'] = rdpoints_reader


# ---RDIELCOR---------------------------------------------------------------- #


# descriptor
class rdielcor_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    icoref = tb.Int32Col(pos=1)
    igrid = tb.Int32Col(pos=2)

descriptors['RDIELCOR'] = rdielcor_rec

atoms['RDIELCOR'] = [('gamma', tb.Float32Atom(), ()),
                     ('beta', tb.Float32Atom(), ()),
                     ('alpha', tb.Float32Atom(), ()),
                     ('layers', tb.Float32Atom(), (2,))]


def rdielcor_reader(rec_name, rec, in_file, out_file):

    table, gamma_arr, beta_arr, alpha_arr, layers_arr = get_table(rec_name,
                                                                  out_file)

    gammas = vldata(gamma_arr)
    betas = vldata(beta_arr)
    alphas = vldata(alpha_arr)
    layers = vldata(layers_arr)
    recs = recdata(table, gammas, betas, alphas, layers)

    while rec_name == 'RDIELCOR':

        icoref = int(rec[1])

        # get layer and plane data
        try:  # if RDPOINTS data has been read
            nlay, nok, noj, noi = rdpoints_data[icoref]
        except KeyError:  # no rdpoints data, assume icoref==eltyp
            if icoref == 15:  # 3D 2 Node Beam
                nlay, nok, noj, noi = 0, 1, 1, 3
            elif icoref == 24:  # Flat Quadrilateral Thin Shell
                nlay, nok, noj, noi = 0, 2, 3, 3

        # gamma
        if nok == 1:
            gammas.add(rec[3:])
        elif nok > 1:
            gammas.add(rec[3:] + in_file.read_floats(nok-1))
        else:
            gammas.add(None)

        betas.add(in_file.read_floats(noj))
        alphas.add(in_file.read_floats(noi))

        l = []
        for _ in range(nlay):
            l.append(in_file.read_floats(2))

        layers.add(l)

        recs.add(rec[:3])

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['RDIELCOR'] = rdielcor_reader


# ---RDRESREF---------------------------------------------------------------- #


# descriptor
class rdresref_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ires = tb.Int32Col(pos=1)
    irno = tb.Int32Col(pos=2)
    ieres = tb.Int32Col(pos=3)
    icalty = tb.Int32Col(pos=4)
    icompl = tb.Int32Col(pos=5)
    numtyp = tb.Int32Col(pos=6)

descriptors['RDRESREF'] = rdresref_rec

atoms['RDRESREF'] = [('reftyps', tb.Float32Atom(), (3,))]


def rdresref_reader(rec_name, rec, in_file, out_file):

    table, reftyps_arr = get_table(rec_name, out_file)
    reftyps = vldata(reftyps_arr)
    recs = recdata(table, reftyps)

    while rec_name == 'RDRESREF':

        rec += in_file.read_floats(3)

        ires = int(rec[1])
        icompl = int(rec[5])
        numtyp = int(rec[6])

        # store complex flag for later use
        rdresref_data[ires] = icompl

        r = []
        for _ in range(numtyp):
            r.append(in_file.read_floats(3))

        reftyps.add(r)
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['RDRESREF'] = rdresref_reader


# ---RBLODCMB, RDRESCMB------------------------------------------------------ #


# descriptor
class rxxxxcmb_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ires = tb.Int32Col(pos=1)
    icompl = tb.Int32Col(pos=2)
    nblc = tb.Int32Col(pos=3)

descriptors['RBLODCMB'] = rxxxxcmb_rec
descriptors['RDRESCMB'] = rxxxxcmb_rec

atoms['RBLODCMB'] = [('blcs', tb.Float32Atom(), (3,))]
atoms['RDRESCMB'] = [('rcs', tb.Float32Atom(), (3,))]


def rxxxxcmb_reader(rec_name, rec, in_file, out_file):

    table, blcs_arr = get_table(rec_name, out_file)
    blcs = vldata(blcs_arr)
    recs = recdata(table, blcs)

    rxxxxcmb_rec_name = rec_name
    while rec_name == rxxxxcmb_rec_name:

        nblc = int(rec[3])

        blcs.add([in_file.read_floats(3) for _ in range(nblc)])
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['RBLODCMB'] = rxxxxcmb_reader
readers['RDRESCMB'] = rxxxxcmb_reader


# ---RVNODACC, RVNODDIS, RVNODVEL-------------------------------------------- #


class rvnodres_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ires = tb.Int32Col(pos=1)
    iinod = tb.Int32Col(pos=2)
    irdva = tb.Int32Col(pos=3)
    itrans = tb.Int32Col(pos=4)

descriptors['RVNODACC'] = rvnodres_rec
descriptors['RVNODDIS'] = rvnodres_rec
descriptors['RVNODVEL'] = rvnodres_rec

atoms['RVNODACC'] = [('acc', tb.Float32Atom(), ())]
atoms['RVNODDIS'] = [('dis', tb.Float32Atom(), ())]
atoms['RVNODVEL'] = [('vel', tb.Float32Atom(), ())]

# each row in a 'rvnodres' vlarray will have nres 32-bit floats, where
# nres = (nfield - 5) / (icompl + 1). If results are complex the dtype should
# be changed from float32 to complex64 before use


def rvnodres_reader(rec_name, rec, in_file, out_file):

    table, res_arr = get_table(rec_name, out_file)
    res = vldata(res_arr)
    recs = recdata(table, res)

    rvnodres_rec_name = rec_name
    while rec_name == rvnodres_rec_name:

        rec += in_file.read_floats(1)
        nfield = int(rec[0])
        ires = int(rec[1])
        icompl = rdresref_data[ires]
        nres = (nfield - 5) // (icompl + 1)

        res.add(in_file.read_floats(nres))
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

readers['RVNODACC'] = rvnodres_reader
readers['RVNODDIS'] = rvnodres_reader
readers['RVNODVEL'] = rvnodres_reader


# ---RVFORCES, RVSTRAIN, RVSTRESS-------------------------------------------- #


class rvelmres_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ires = tb.Int32Col(pos=1)
    iielno = tb.Int32Col(pos=2)
    ispalt = tb.Int32Col(pos=3)
    ir = tb.Int32Col(pos=4)

# ir = irforc for RVFORCES etc..

descriptors['RVFORCES'] = rvelmres_rec
descriptors['RVSTRAIN'] = rvelmres_rec
descriptors['RVSTRESS'] = rvelmres_rec

atoms['RVFORCES'] = [('force', tb.Float32Atom(), ())]
atoms['RVSTRAIN'] = [('strain', tb.Float32Atom(), ())]
atoms['RVSTRESS'] = [('stress', tb.Float32Atom(), ())]

# each row in a 'rvelmres' vlarray will have n 32-bit floats, where
# n = nsp*lenrec*(icompl + 1). If results are complex the dtype should
# be changed from float32 to complex64 before use.

# Reshaping the data to (nsp, lenrec)
# corresponds to a (respoint, rescomp) shape

def rvelmres_reader(rec_name, rec, in_file, out_file):

    table, res_arr = get_table(rec_name, out_file)
    res = vldata(res_arr)
    recs = recdata(table, res)

    rvnodres_rec_name = rec_name
    while rec_name == rvnodres_rec_name:

        rec += in_file.read_floats(1)

        ires = int(rec[1])
        iielno = int(rec[2])
        ispalt = int(rec[3])
        ir = int(rec[4])

        nsp = nsp_data[(ispalt, iielno)]
        lenrec = lenrec_data[(rec_name, ir)]
        icompl = rdresref_data[ires]

        res.add(in_file.read_floats(nsp*lenrec*(icompl + 1)))
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec

readers['RVFORCES'] = rvelmres_reader
readers['RVSTRAIN'] = rvelmres_reader
readers['RVSTRESS'] = rvelmres_reader


# ---WBODCON----------------------------------------------------------------- #


class wbodcon_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibcond = tb.Int32Col(pos=1)
    ibody = tb.Int32Col(pos=2)
    numcon = tb.Int32Col(pos=3)

descriptors['WBODCON'] = wbodcon_rec

atoms['WBODCON'] = [('refcond', tb.Float32Atom(), (3,))]


def wbodcon_reader(rec_name, rec, in_file, out_file):

    table, refcond_arr = get_table(rec_name, out_file)
    refconds = vldata(refcond_arr)
    recs = recdata(table, refconds)

    while rec_name == 'WBODCON':
        r = []
        for _ in range(int(rec[3])):
            r.append(in_file.read_floats(3))
        refconds.add(r)
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()

    recs.write()
    return rec_name, rec


readers['WBODCON'] = wbodcon_reader


# ---WDRESREF---------------------------------------------------------------- #


class wdresref_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    iwres = tb.Int32Col(pos=1)
    nresrf = tb.Int32Col(pos=2)
    numtyp = tb.Int32Col(pos=3)

descriptors['WDRESREF'] = wdresref_rec

atoms['WDRESREF'] = [('dir', tb.Float32Atom(), (2,)),
                     ('freq', tb.Float32Atom(), (2,)),
                     ('time', tb.Float32Atom(), (2,))]


def wdresref_reader(rec_name, rec, in_file, out_file):

    table, dir_arr, freq_arr, time_arr = get_table(rec_name, out_file)
    dirs = vldata(dir_arr)
    freqs = vldata(freq_arr)
    times = vldata(time_arr)
    recs = recdata(table, dirs, freqs, times)
    while rec_name == 'WDRESREF':
        for _ in range(int(rec[3])):
            irefty = int(in_file.read_floats(1)[0])
            dat = []
            for __ in range(int(rec[2])):
                dat.append(in_file.read_floats(2))
            if irefty == 1:
                dirs.add(dat)
            elif irefty == 2:
                freqs.add(dat)
            elif irefty == 3:
                times.add(dat)
        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec


readers['WDRESREF'] = wdresref_reader


# ---WGLOBDEF---------------------------------------------------------------- #


class wglobdef_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    g = tb.Float32Col(pos=1)
    ro = tb.Float32Col(pos=2)
    idepth = tb.Int32Col(pos=3)
    depth = tb.Float32Col(pos=4)

descriptors['WGLOBDEF'] = wglobdef_rec


def wglobdef_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)
    while rec_name == 'WGLOBDEF':
        if int(rec[3]):
            rec += in_file.read_floats(1)
        else:
            rec += (1.,)
        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec


readers['WGLOBDEF'] = wglobdef_reader


# ---WBODY------------------------------------------------------------------- #


class wbody_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibcond = tb.Int32Col(pos=1)
    clenb = tb.Float32Col(pos=2)
    bmass = tb.Float32Col(pos=3)
    wsarea = tb.Float32Col(pos=4)
    wparea = tb.Float32Col(pos=5)
    vol = tb.Float32Col(pos=6)
    xcg = tb.Float32Col(pos=7)
    ycg = tb.Float32Col(pos=8)
    zcg = tb.Float32Col(pos=9)
    xcb = tb.Float32Col(pos=10)
    ycb = tb.Float32Col(pos=11)
    zcb = tb.Float32Col(pos=12)
    xrb = tb.Float32Col(pos=13)
    yrb = tb.Float32Col(pos=14)
    zrb = tb.Float32Col(pos=15)
    t = tb.Float32Col(pos=16, shape=(4, 4))

descriptors['WBODY'] = wbody_rec

def wbody_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)
    while rec_name == 'WBODY':
        rec += in_file.read_floats(12)
        t = in_file.read_floats(16)
        rec += ((t[:4], t[4:8], t[8:12], t[12:]),)
        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['WBODY'] = wbody_reader


# ---W1MATRIX---------------------------------------------------------------- #


class w1matrix_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    imatrx = tb.Int32Col(pos=1)
    ibconi = tb.Int32Col(pos=2)
    ibconj = tb.Int32Col(pos=3)
    iwres = tb.Int32Col(pos=4)
    imtyp = tb.Int32Col(pos=5)
    hmat = tb.Float32Col(pos=6, shape=(6, 6))

descriptors['W1MATRIX'] = w1matrix_rec

def w1matrix_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)
    while rec_name == 'W1MATRIX':
        rec += in_file.read_floats(2)
        t = in_file.read_floats(36)
        rec += ((t[:6], t[6:12], t[12:18], t[18:24], t[24:30], t[30:]),)
        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['W1MATRIX'] = w1matrix_reader


# ---W1EXFORC, W1MOTION, W2EXFDIF, W2EXFSUM---------------------------------- #


rec_names = ('W1EXFORC', 'W1MOTION', 'W2EXFDIF', 'W2EXFSUM')

class wcommon1_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibcond = tb.Int32Col(pos=1)
    iwres = tb.Int32Col(pos=2)
    icompl = tb.Int32Col(pos=3)
    data = tb.ComplexCol(itemsize=8, pos=4, shape=(6,))


def wcommon1_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)

    wcommon1_rec_name = rec_name
    while rec_name == wcommon1_rec_name:
        if int(rec[3]):
            f = in_file.read_floats(12)
            rec += (tuple(complex(x,y) for x,y in zip(f[::2], f[1::2])),)
        else:
            rec += (in_file.read_floats(6),)

        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

for rec_name in rec_names:
    descriptors[rec_name] = wcommon1_rec
    readers[rec_name] = wcommon1_reader


# ---W1SFORCE---------------------------------------------------------------- #


class W1SFORCE_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibcond = tb.Int32Col(pos=1)
    isect = tb.Int32Col(pos=2)
    iwres = tb.Int32Col(pos=3)
    icompl = tb.Int32Col(pos=4)
    data = tb.ComplexCol(itemsize=8, pos=5, shape=(6,))

descriptors['W1SFORCE'] = W1SFORCE_rec


def W1SFORCE_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)

    while rec_name == 'W1SFORCE':
        rec += in_file.read_floats(1)
        if int(rec[4]):
            f = in_file.read_floats(12)
            rec += (tuple(complex(x,y) for x,y in zip(f[::2], f[1::2])),)
        else:
            rec += (in_file.read_floats(6),)

        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['W1SFORCE'] = W1SFORCE_reader


# ---W2HDRIFT---------------------------------------------------------------- #


class W2HDRIFT_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibcond = tb.Int32Col(pos=1)
    iwres = tb.Int32Col(pos=2)
    shfor = tb.Float32Col(pos=3, shape=(3,))

descriptors['W2HDRIFT'] = W2HDRIFT_rec


def W2HDRIFT_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)
    while rec_name == 'W2HDRIFT':
        recs.add(rec[:3] + (rec[3:] + in_file.read_floats(2),))
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['W2HDRIFT'] = W2HDRIFT_reader


# ---W2MDRIFT---------------------------------------------------------------- #


class W2MDRIFT_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibcond = tb.Int32Col(pos=1)
    iwres = tb.Int32Col(pos=2)
    smfor = tb.Float32Col(pos=3, shape=(6,))

descriptors['W2MDRIFT'] = W2MDRIFT_rec


def W2MDRIFT_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)
    while rec_name == 'W2MDRIFT':
        recs.add(rec[:3] + (rec[3:] + in_file.read_floats(5),))
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['W2MDRIFT'] = W2MDRIFT_reader


# ---W2FLUDIF, W2FLUSUM------------------------------------------------------ #

# not implemented..

# ---WFKPOINT---------------------------------------------------------------- #


class WFKPOINT_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=1)
    ifkpnt = tb.Int32Col(pos=2)
    fkpnt = tb.Float32Col(pos=3, shape=(3,))

descriptors['WFKPOINT'] = WFKPOINT_rec

def WFKPOINT_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)
    while rec_name == 'WFKPOINT':
        p = rec[2:] + in_file.read_floats(1)
        recs.add(rec[:2] + (p,))
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['WFKPOINT'] = WFKPOINT_reader


# ---WFLUIDKN---------------------------------------------------------------- #


class WFLUIDKN_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    iwres = tb.Int32Col(pos=1)
    ifkpnt = tb.Int32Col(pos=2)
    icompl = tb.Int32Col(pos=3)
    p = tb.ComplexCol(itemsize=8, pos=4)
    v = tb.ComplexCol(itemsize=8, pos=5, shape=(3,))

descriptors['WFLUIDKN'] = WFLUIDKN_rec


def WFLUIDKN_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)

    while rec_name == 'WFLUIDKN':
        if int(rec[3]):
            f = in_file.read_floats(2)
            rec += (complex(f[0], f[1]),)
            f = in_file.read_floats(6)
            rec += (tuple(complex(x,y) for x,y in zip(f[::2], f[1::2])),)
        else:
            rec += in_file.read_floats(1) + (in_file.read_floats(3),)

        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['WFLUIDKN'] = WFLUIDKN_reader


# ---WGRESPON---------------------------------------------------------------- #


# ---WSCATTER---------------------------------------------------------------- #


# ---WSECTION---------------------------------------------------------------- #


class WSECTION_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibody = tb.Int32Col(pos=1)
    isect = tb.Int32Col(pos=2)
    isecty = tb.Int32Col(pos=3)
    secusr = tb.Int32Col(pos=4)
    p1 = tb.Float32Col(pos=5, shape=(3,))
    p2 = tb.Float32Col(pos=6, shape=(3,))
    p3 = tb.Float32Col(pos=7, shape=(3,))

descriptors['WSECTION'] = WSECTION_rec

def WSECTION_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)
    while rec_name == 'WSECTION':
        rec += in_file.read_floats(10)
        recs.add(rec[:5] + (rec[5:8], rec[8:11], rec[11:14]))
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['WSECTION'] = WSECTION_reader

# ---WSURFACE---------------------------------------------------------------- #


class WSURFACE_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    iwres = tb.Int32Col(pos=1)
    indxwr = tb.Int32Col(pos=2)
    nnode = tb.Int32Col(pos=3)
    icompl = tb.Int32Col(pos=4)
    icompx = tb.Int32Col(pos=4)
    icompy = tb.Int32Col(pos=4)
    icompz = tb.Int32Col(pos=4)

descriptors['WSURFACE'] = WSURFACE_rec

atoms['WSURFACE'] = [('zelev', tb.ComplexAtom(itemsize=8), (4,))]


def WSURFACE_reader(rec_name, rec, in_file, out_file):

    table, zelev_arr = get_table(rec_name, out_file)
    zelevs = vldata(zelev_arr)
    recs = recdata(table, zelevs)

    while rec_name == 'WSURFACE':

        rec += in_file.read_floats(4)
        nfloats_tot = int(rec[0]) - 8
        nfloats_per = 1 + (int(rec[5] + rec[6] + rec[7])) * (1 + int(rec[4]))
        all_floats = in_file.read_floats(nfloats_tot)

        zelev_recs = []
        for start in range(0, nfloats_tot, nfloats_per):
            fl_it = iter(all_floats[start:start + nfloats_per])
            zelev_rec = (next(fl_it),)
            for comp in rec[5:]:
                if comp:
                    if rec[4]:
                        zelev_rec += (complex(next(fl_it), next(fl_it)),)
                    else:
                        zelev_rec += (next(fl_it),)
                else:
                    zelev_rec += (0.,)
            zelev_recs.append(zelev_rec)
        zelevs.add(zelev_recs)
        recs.add(rec)

        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['WSURFACE'] = WSURFACE_reader


# ---WINPUT------------------------------------------------------------------ #


class WINPUT_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibcond = tb.Int32Col(pos=1)  # ?? assume this field is ibcond
    t = tb.Float32Col(pos=2, shape=(4, 4))

descriptors['WINPUT'] = WINPUT_rec

def WINPUT_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)
    while rec_name == 'WINPUT':
        t = rec[2:] + in_file.read_floats(14)
        recs.add(rec[:2] + ((t[:4], t[4:8], t[8:12], t[12:]),))
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec

readers['WINPUT'] = WINPUT_reader

# ---W1PANPRE----------------------------------------------------------------- #

class W1PANPRE_rec(tb.IsDescription):
    nfield = tb.Int32Col(pos=0)
    ibcond = tb.Int32Col(pos=1)
    iwres = tb.Int32Col(pos=2)
    ippan = tb.Int32Col(pos=3)
    isymm1 = tb.Int32Col(pos=4)
    icompl = tb.Int32Col(pos=5)
    p = tb.ComplexCol(itemsize=8, pos=6)
    ipan = tb.Int32Col(pos=9)
    isymm2 = tb.Int32Col(pos=10)


descriptors['W1PANPRE'] = W1PANPRE_rec

def W1PANPRE_reader(rec_name, rec, in_file, out_file):

    table = get_table(rec_name, out_file)
    recs = recdata(table)

    while rec_name == 'W1PANPRE':
        f = in_file.read_floats(2)
        rec += f
        if int(f[1]):
            f = in_file.read_floats(2)
            rec += (complex(f[0], f[1]),)
        else:
            f = in_file.read_floats(1)
            rec += (complex(f[0], 0),)
        f = in_file.read_floats(2)
        rec += f

        recs.add(rec)
        rec_name, rec = in_file.read_headerrec()
    recs.write()
    return rec_name, rec


readers['W1PANPRE'] = W1PANPRE_reader
