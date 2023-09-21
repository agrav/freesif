# -*- coding: utf-8 -*-
# Copyright (c) 2015 Audun Gravdal Johansen
"""This module contains file types used for reading/writing SESAM sequential
files.
"""

from __future__ import division

import struct
from math import ceil
import io

# Set with allowable record names
allowed_record_names_str = set([
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
'RDMLFACT', 'RDNODBOC', 'RDNODREA', 'RDNODRES', 'RDPOINTS', 'RDRESCMB',
'RDRESREF', 'RDSERIES', 'RDSTRAIN', 'RDSTRESS', 'RDTRANS', 'RSUMLOAD',
'RSUMMASS', 'RSUMREAC', 'RSUPTRAN', 'RVABSCIS', 'RVELNFOR', 'RVFORCES',
'RVNODACC', 'RVNODDIS', 'RVNODVEL', 'RVNODREA', 'RVORDINA', 'RVSTRAIN',
'RVSTRESS', 'TNODE', 'TELEMENT', 'TDRESREF', 'TDSERIES', 'TDSUPNAM',
'W1EXFORC', 'W1MATRIX', 'WIMOTION','WISFORCE', 'W2EXFDIF', 'W2EXFSUM',
'W2FLUDIF', 'W2FLUSUM', 'W2HDRIFT', 'W2MDRIFI', 'WBODCON', 'WBODY', 'WDRESREF',
'WFKPOINT', 'WFLUIDKN', 'WGLOBDEF', 'WGRESPON', 'WSCATTER', 'WSECTION',
'WSURFACE', 'TDBODNAM', 'TDRSNAM', 'TDSCATTER', 'TDSCONC', 'SCONCEPT',
'SCONMESH', 'SCONPLIS', 'SPROSELE', 'SPROMATR', 'SPROSEGM', 'SPROHYDR',
'SPROCODE', 'SPROORIE', 'SPROECCE', 'SPROPILE', 'SPROSOIL', 'IRECSIZE', 'W1PANPRE'])

allowed_record_names_b = set([
b'DATE', b'IDENT', b'IEND', b'TDMATER', b'TDSECT', b'TDSETNAM', b'TDSUPNAM',
b'TEXT', b'TSLAYER', b'ACFD', b'ADDATA', b'BEISTE', b'BELFIX', b'BELLAX',
b'BELLO2', b'BELOAD1', b'BEDRAG1', b'BEMASS1', b'BEUSLO', b'BEUVLO', b'BEWAKIN',
b'BEWALO1', b'BGRAV', b'BLDEP', b'BLDEP', b'BNACCLO', b'BNBCD', b'BNDISPL',
b'BNDOF', b'BNINCO', b'BNLOAD', b'BNLOAX', b'BNMASS', b'BNTEMP', b'BNTRCOS',
b'BNWALO', b'BRIGAC', b'BRIGDI', b'BRIGVE', b'BQDP', b'GBARM', b'GBEAMG',
b'GBOX', b'GCHAN', b'GCHANR', b'GCOORD', b'GCROINT', b'GDOBO', b'GECC',
b'GECCEN', b'GELINT', b'GELMNT1', b'GELREF1', b'GELSTRP', b'GELTH', b'GIORH',
b'GIORHR', b'GLMASS', b'GLSEC', b'GLSECR', b'GNODE', b'GPIPE', b'GSEPSPEC',
b'GSETMEMB', b'GSLAYER', b'GSLPLATE', b'GSLSTIFF', b'GTONP', b'GUNIVEC', b'GUSYI',
b'MAXDMP', b'MAXSPR', b'MCNT', b'MGDAMP', b'MGLDAMP', b'MGLMASS', b'MGMASS',
b'MGSPRNG', b'MISOAL', b'MISOEML', b'MISOHL', b'MISOHNL', b'MISOPL', b'MISOPL',
b'MISOPL', b'MISOSEL', b'MISTEL', b'MORSMEL', b'MORSSEL', b'MORSSOL', b'MSHGLSP',
b'MTEMP', b'MTENONL', b'MTRMEL', b'MTRSEL', b'MTRSOL', b'ADDATA', b'AMATRIX',
b'AMDACCL', b'AMDDAMP', b'AMDDISP', b'AMDFREQ', b'AMDLOAD', b'AMDMASS', b'AMDSTIFF',
b'AMDVELO', b'BLDEP', b'BNBCD', b'BNDISPL', b'BNDOF', b'BNINCO', b'BNLOAD',
b'BNMASS', b'BNTRCOS', b'BQDP', b'BSELL', b'GCOORD', b'GELMNT1', b'GELMNT2',
b'GELREF1', b'GNODE', b'HIERARCH', b'HSUPSTAT', b'HSUPTRAN', b'MAXDMP', b'MAXSPR',
b'MGDAMP', b'MGSPRNG', b'RBLODCMB', b'RDELNFOR', b'RDFORCES', b'RDIELCOR',
b'RDMLFACT', b'RDNODBOC', b'RDNODREA', b'RDNODRES', b'RDPOINTS', b'RDRESCMB',
b'RDRESREF', b'RDSERIES', b'RDSTRAIN', b'RDSTRESS', b'RDTRANS', b'RSUMLOAD',
b'RSUMMASS', b'RSUMREAC', b'RSUPTRAN', b'RVABSCIS', b'RVELNFOR', b'RVFORCES',
b'RVNODACC', b'RVNODDIS', b'RVNODVEL', b'RVNODREA', b'RVORDINA', b'RVSTRAIN',
b'RVSTRESS', b'TNODE', b'TELEMENT', b'TDRESREF', b'TDSERIES', b'TDSUPNAM',
b'W1EXFORC', b'W1MATRIX', b'WIMOTION','WISFORCE', b'W2EXFDIF', b'W2EXFSUM',
b'W2FLUDIF', b'W2FLUSUM', b'W2HDRIFT', b'W2MDRIFI', b'WBODCON', b'WBODY', b'WDRESREF',
b'WFKPOINT', b'WFLUIDKN', b'WGLOBDEF', b'WGRESPON', b'WSCATTER', b'WSECTION',
b'WSURFACE', b'TDBODNAM', b'TDRSNAM', b'TDSCATTER', b'TDSCONC', b'SCONCEPT',
b'SCONMESH', b'SCONPLIS', b'SPROSELE', b'SPROMATR', b'SPROSEGM', b'SPROHYDR',
b'SPROCODE', b'SPROORIE', b'SPROECCE', b'SPROPILE', b'SPROSOIL', b'IRECSIZE', b'W1PANPRE'])


class SequentialFile(object):
    """Base class for FormattedFile and UnformattedFile
    """
    # methods that need to be implemented in derived classes:

    # read_string_rec(self):
    # tofloatrec(self, string_rec):
    # read_floatrec(self, skipfirst=0):
    # read_headerrec(self):
    # read_floats(self, n, skipfirst=0):
    # at_headerrec(self):

    def __init__(self, name, mode='r'):
        self.f = io.open(name, mode)

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
            for i in range(int(rec[2]) // 100):
                self.read_stringrec()
            for i in range(int(rec[3]) // 100):
                self.read_stringrec()
            return self.read_headerrec()

        #read first record
        pos = self.f.tell()  # current position in in_file
        string_rec = self.read_stringrec()

        if string_rec[:8] == self._contrec_leadingspace:  # type 1 cont. record
#            while string_rec[:8] == self._contrec_leadingspace:  # read cont. recs
#                string_rec = self.read_stringrec()
#            return string_rec[:8].decode().strip(), self.tofloatrec(string_rec[8:])
            return self.skip_contrecs(string_rec)

        elif string_rec[:8].strip() in self.allowed_record_names: # no cont. records
            return self.strip_recname(string_rec[:8]), self.tofloatrec(string_rec[8:])

        else:  # type 2 record
            nfloats = int(rec[0]) - 4  # nfield - 4
            self.f.seek(pos)  # rewind
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

            if rec_name in recs:
                recs[rec_name] += 1
            else:
                recs[rec_name] = 1
        return recs

    def close(self):
        self.f.close()


class UnformattedFile(SequentialFile):

    def __init__(self, fname, mode='r'):
        super(UnformattedFile, self).__init__(fname, mode+'b')
        if mode == 'r':
            self.f.seek(1)  # skip 1st byte (the 'K')
        elif mode == 'w':
            self.f.write(b'K')  # write mandatory(?) 1st byte
        self.allowed_record_names = allowed_record_names_b
        self._contrec_leadingspace = b'        '

    def read_stringrec(self):

        self.hasrest = False

        rlen = struct.unpack('B', self.f.read(1))[0]
        if rlen == 129:
            rlen = 128
        rec = self.f.read(rlen)
        self.f.read(1)
        return rec

    def tofloatrec(self, stringrec):
        return struct.unpack('{}f'.format(len(stringrec)//4), stringrec)

    def read_floatrec(self, skipfirst=0):
        return self.tofloatrec(self.read_stringrec()[skipfirst:])

    def read_headerrec(self):
        stringrec = self.read_stringrec()
        return stringrec[:8].decode().strip(), self.tofloatrec(stringrec[8:])

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

        pos = self.f.tell()
#        rec_name = self.f.read(9)[1:].strip().decode()
        rec_name = self.f.read(9)[1:].strip()
        self.f.seek(pos)

        if not self.hasrest and rec_name in self.allowed_record_names:
            return True

        return False

    def strip_recname(self, rec_name):
        return rec_name.decode().strip()

    def skip_contrecs(self, string_rec):
        while string_rec[:8] == self._contrec_leadingspace:
            string_rec = self.read_stringrec()
        return self.strip_recname(string_rec[:8]), self.tofloatrec(string_rec[8:])



class FormattedFile(SequentialFile):

    def __init__(self, fname, mode='r'):

        super(FormattedFile, self).__init__(fname, mode)
        # remaining float places on current record (used for write mode)
        self.float_places_left = 0
        self.allowed_record_names = allowed_record_names_str
        self._contrec_leadingspace = '        '


    def read_stringrec(self):
        self.hasrest = False

        return self.f.readline()[:-1]  # skip the newline character

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

            for _ in range(nrecs):
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

        pos = self.f.tell()
        rec_name = self.f.read(8).strip()
        self.f.seek(pos)

        if not self.hasrest and rec_name in self.allowed_record_names:
            return True

        return False

    def strip_recname(self, rec_name):
        return rec_name.strip()

    def skip_contrecs(self, string_rec):
        while string_rec[:8] == self._contrec_leadingspace:
            string_rec = self.read_stringrec()
        return self.strip_recname(string_rec[:8]), self.tofloatrec(string_rec[8:])

    def _write_line(self, header, rec):
        if self.float_places_left:
            # fill with blanks:
            # self.write(' '*self.float_places_left*16 + '\n')
            self.f.write('\n')

        self.write('{:<8}'.format(header))
        for f in rec:
            self.f.write('{:>16.8e}'.format(f))

        self.float_places_left = 4 - len(rec)
        if not self.float_places_left:
            self.f.write('\n')

    def write_string(self, string, leading_blanks=8):
        if self.float_places_left:
            self.f.write('\n')

        self.f.write(' '*leading_blanks + string + '\n')

    def write_floats(self, floats):
        if len(floats) == 0:
            return

        fpl = self.float_places_left

        if fpl:  # places left on current line
            for f in floats[:fpl]:
                self.f.write('{:>16.8e}'.format(f))

            if fpl > len(floats):  # still float places left?
                self.float_places_left = fpl - len(floats)
                return
            elif fpl == len(floats):
                self.float_places_left = 0
                self.f.write('\n')
                return

        n_floats = len(floats[fpl:])
        n_lines = n_floats // 4
        n_rest = n_floats % 4
        if n_rest:
            n_lines += 1

        # write new lines
        for i in range(0, 4*n_lines, 4):
            self._write_line('        ', floats[fpl:][i:i+4])

        # finally, update self.float_places_left
        # self.float_places_left = 4 - n_rest if n_rest else 0

    def write_record(self, rec_name, rec):
        self._write_line(rec_name, rec[:4])
        self.write_floats(rec[4:])
