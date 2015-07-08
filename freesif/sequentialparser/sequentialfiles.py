# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""This module contains file types used for reading/writing SESAM sequential
files.
"""

import struct
from math import ceil

# Set with allowable record names
allowed_record_names = set([
'DATE', 'IDENT', 'IEND', 'TDMATER', 'TDSECT', 'TDSETNAM', 'TDSUPNAM',
'TEXT', 'TSLAYER', 'ACFD', 'ADDATA', 'BEISTE', 'BELFIX', 'BELLAX',
'BELLO2', 'BELOAD1', 'BEDRAG1', 'BEMASS1', 'BEUSLO', 'BEUVLO', 'BEWAKIN',
'BEWALO1', 'BGRAV', 'BLDEP', 'BLDEP', 'BNACCLO', 'BNBCD', 'BNDISPL',
'BNDOF', 'BNINCO', 'BNLOAD', 'BNLOAX', 'BNMASS', 'BNTEMP', 'BNTRCOS',
'BNWALO', 'BRIGAC', 'BRIGDI', 'BRIGVE', 'BQDP', 'GBARM', 'GBEAMG',
'GBOX', 'GCHAN', 'GCHANR', 'GCOORD', 'GCROINT', 'GDOBO', 'GECC',
'GECCEN', 'GELINT', 'GELMNT1', 'GELREF1', 'GELSTRP', 'GELTH', 'GIORH',
'GIORHR', 'GLMASS', 'GLSEC', 'GLSECR', 'GNODE', 'GPIPE', 'GSEPSPEC',
'GSETMEMB', 'GSLAYER', 'GSLPLATE', 'GSLSTIFF', 'GTONP', 'GUNIVEC', 'GUSYI',
'MAXDMP', 'MAXSPR', 'MCNT', 'MGDAMP', 'MGLDAMP', 'MGLMASS', 'MGMASS',
'MGSPRNG', 'MISOAL', 'MISOEML', 'MISOHL', 'MISOHNL', 'MISOPL', 'MISOPL',
'MISOPL', 'MISOSEL', 'MISTEL', 'MORSMEL', 'MORSSEL', 'MORSSOL', 'MSHGLSP',
'MTEMP', 'MTENONL', 'MTRMEL', 'MTRSEL', 'MTRSOL', 'ADDATA', 'AMATRIX',
'AMDACCL', 'AMDDAMP', 'AMDDISP', 'AMDFREQ', 'AMDLOAD', 'AMDMASS', 'AMDSTIFF',
'AMDVELO', 'BLDEP', 'BNBCD', 'BNDISPL', 'BNDOF', 'BNINCO', 'BNLOAD',
'BNMASS', 'BNTRCOS', 'BQDP', 'BSELL', 'GCOORD', 'GELMNT1', 'GELMNT2',
'GELREF1', 'GNODE', 'HIERARCH', 'HSUPSTAT', 'HSUPTRAN', 'MAXDMP', 'MAXSPR',
'MGDAMP', 'MGSPRNG', 'RBLODCMB', 'RDELNFOR', 'RDFORCES', 'RDIELCOR',
'RDMLFACT', 'RDNODRES', 'RDPOINTS', 'RDRESCMB', 'RDRESREF', 'RDSERIES',
'RDSTRAIN', 'RDSTRESS', 'RDTRANS', 'RSUMLOAD', 'RSUMMASS', 'RSUPTRAN',
'RVABSCIS', 'RVELNFOR', 'RVFORCES', 'RVNODACC', 'RVNODDIS', 'RVNODVEL',
'RVORDINA', 'RVSTRAIN', 'RVSTRESS', 'TNODE', 'TELEMENT', 'TDRESREF',
'TDSERIES', 'TDSUPNAM', 'W1EXFORC', 'W1MATRIX', 'WIMOTION','WISFORCE',
'W2EXFDIF', 'W2EXFSUM', 'W2FLUDIF', 'W2FLUSUM', 'W2HDRIFT', 'W2MDRIFI',
'WBODCON', 'WBODY', 'WDRESREF', 'WFKPOINT', 'WFLUIDKN', 'WGLOBDEF',
'WGRESPON', 'WSCATTER', 'WSECTION', 'WSURFACE', 'TDBODNAM', 'TDRSNAM',
'TDSCATTER', 'TDSCONC', 'SCONCEPT', 'SCONMESH', 'SCONPLIS', 'SPROSELE',
'SPROMATR', 'SPROSEGM', 'SPROHYDR', 'SPROCODE', 'SPROORIE', 'SPROECCE',
'SPROPILE', 'SPROSOIL', 'IRECSIZE'])


class SequentialFile(file):
    """Base class for FormattedFile and UnformattedFile
    """
    # methods that need to be implemented in derived classes:

    # read_string_rec(self):
    # tofloatrec(self, string_rec):
    # read_floatrec(self, skipfirst=0):
    # read_headerrec(self):
    # read_floats(self, n, skipfirst=0):
    # at_headerrec(self):

    # common method implementations:
    def skip_restofrecord(self, rec_name, rec):
        """
        skips rest of the current record and returns the next header record
        """

#         two types of records:
#          1. No 'nfield' first value
#             Continuation records have 8 leading whitespaces
#             Continuation records are have max 4 floats
#             (input records)

#          2. The record has an 'nfield' as its first value
#             Will NOT have leading whitespaces on continuation records
#             Can have up to 32 floats on continuation records
#             (results records)

#         Note: This is only relevant for unformatted files. Formatted files
#               will always have 8 leading whitespaces and max 4 floats on
#                continuation records.

        # first float is negative:
        if rec[0] < 0.:
            if rec[0] == -5.:
                self.read_floatrec()
            return self.read_headerrec()

        # workaround to handle TDxxxxxx records (on unformatted files):
        if rec_name[:2] == 'TD':
            for i in range(int(rec[2]) / 100):
                self.read_stringrec()
            for i in range(int(rec[3]) / 100):
                self.read_stringrec()
            return self.read_headerrec()

        #read first record
        pos = self.tell()  # current position in in_file
        string_rec = self.read_stringrec()

        if string_rec[:8] == '        ':  # type 1 cont. record
            while string_rec[:8] == '        ':  # read cont. recs
                string_rec = self.read_stringrec()
            return string_rec[:8].strip(), self.tofloatrec(string_rec[8:])

        elif string_rec[:8].strip() in allowed_record_names: # no cont. records
            return string_rec[:8].strip(), self.tofloatrec(string_rec[8:])

        else:  # type 2 record
            nfloats = int(rec[0]) - 4  # nfield - 4
            self.seek(pos)  # rewind
            self.read_floats(nfloats)  # read remainder of record
            return self.read_headerrec()

    def get_recordcount(self):
        recs = {}

        # 1st record
        rec_name, rec = self.read_headerrec()
        recs[rec_name] = 1

        # loop until IEND record
        while rec_name:

            # next header record
            rec_name, rec = self.skip_restofrecord(rec_name, rec)

            if recs.has_key(rec_name):
                recs[rec_name] += 1
            else:
                recs[rec_name] = 1
        return recs



class UnformattedFile(SequentialFile):

    def __init__(self, fname, mode='r'):

        file.__init__(self, fname, mode+'b')

        if mode == 'r':
            self.seek(1)  # skip 1st byte (the 'K')

        elif mode == 'w':
            self.write('K')  # write mandatory(?) 1st byte

    def read_stringrec(self):

        self.hasrest = False

        rlen = struct.unpack('B', self.read(1))[0]
        if rlen == 129:
            rlen = 128
        rec = self.read(rlen)
        self.read(1)
        return rec

    def tofloatrec(self, stringrec):
        return struct.unpack('{}f'.format(len(stringrec)/4), stringrec)

    def read_floatrec(self, skipfirst=0):
        return self.tofloatrec(self.read_stringrec()[skipfirst:])

    def read_headerrec(self):
        stringrec = self.read_stringrec()
        return stringrec[:8].strip(), self.tofloatrec(stringrec[8:])

    def read_floats(self, n, skipfirst=0):
        """returns tuple of n floats"""

        if n == 0:
            return ()  # empty tuple

        if self.hasrest:
            floatrest = self.floatrest
        else:
            floatrest = ()

        floats_to_read = n - len(floatrest)

        if floats_to_read > 0:

            # continuation records are NOT necessary 4 floats
            # on unformatted files. Can be anything up tp 32 (128 Bytes)

            floats = ()

            while len(floats) < floats_to_read:
                floats += self.read_floatrec(skipfirst)

            floats = floatrest + floats

        else:
            floats = floatrest

        self.floatrest = floats[n:]

        if len(self.floatrest) > 0:
            self.hasrest = True

        return floats[:n]

    def at_headerrec(self):
        """check if current pos in file is at start of a header record.
        """

        pos = self.tell()
        rec_name = self.read(9)[1:].strip()
        self.seek(pos)

        if not self.hasrest and rec_name in allowed_record_names:
            return True

        return False


class FormattedFile(SequentialFile):

    def __init__(self, fname, mode='r'):

        file.__init__(self, fname, mode)

        # remaining float places on current record (used for write mode)
        self.float_places_left = 0

    def read_stringrec(self):
        self.hasrest = False

        return self.readline()[:-1]  # skip the newline character

    def tofloatrec(self, stringrec):
        return tuple([float(stringrec[i:i+16])
                     for i in range(0, len(stringrec), 16)])


    def read_floatrec(self, skipfirst=0):
        return self.tofloatrec(self.read_stringrec()[8:])

    def read_headerrec(self):
        """
        returns a tuple (rec_name, rec) where rec_name is the name
        of the record (string) and rec is a tuple of floats (the
        first 4 numeric values on the record).
        """

        stringrec = self.read_stringrec()
        return stringrec[:8].strip(), self.tofloatrec(stringrec[8:])

    def read_floats(self, n, skipfirst=0):
        """returns array of n floats."""

        if n == 0:
            return ()

        if self.hasrest:
            floatrest = self.floatrest
        else:
            floatrest = ()

        floats_to_read = n - len(floatrest)

        if floats_to_read > 0:
            nrecs = int(ceil(floats_to_read/4.))
            floats = ()

            for _ in xrange(nrecs):
                floats += self.read_floatrec(skipfirst)

            floats = floatrest + floats

        else:
            floats = floatrest

        self.floatrest = floats[n:]

        if len(self.floatrest) > 0:
            self.hasrest = True

        return floats[:n]

    def at_headerrec(self):
        """check if current pos in file is at start of a header record.
        """

        pos = self.tell()
        rec_name = self.read(8).strip()
        self.seek(pos)

        if not self.hasrest and rec_name in allowed_record_names:
            return True

        return False

    def _write_line(self, header, rec):
        if self.float_places_left:
            # fill with blanks:
            # self.write(' '*self.float_places_left*16 + '\n')
            self.write('\n')

        self.write('{:<8}'.format(header))
        for f in rec:
            self.write('{:>16.8e}'.format(f))

        self.float_places_left = 4 - len(rec)
        if not self.float_places_left:
            self.write('\n')

    def write_string(self, string, leading_blanks=8):
        if self.float_places_left:
            self.write('\n')

        self.write(' '*leading_blanks + string + '\n')

    def write_floats(self, floats):
        if len(floats) == 0:
            return

        fpl = self.float_places_left

        if fpl:  # places left on current line
            for f in floats[:fpl]:
                self.write('{:>16.8e}'.format(f))

            if fpl > len(floats):  # still float places left?
                self.float_places_left = fpl - len(floats)
                return
            elif fpl == len(floats):
                self.float_places_left = 0
                self.write('\n')
                return

        n_floats = len(floats[fpl:])
        n_lines = n_floats / 4
        n_rest = n_floats % 4
        if n_rest:
            n_lines += 1

        # write new lines
        for i in range(0, 4*n_lines, 4):
            self._write_line('        ', floats[fpl:][i:i+4])

        # finally, update self.float_places_left
        # self.float_places_left = 4 - n_rest if n_rest else 0

    def writeRecord(self, rec_name, rec):
        self._write_line(rec_name, rec[:4])
        self.write_floats(rec[4:])
