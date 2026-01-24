#
#    The high-throughput toolkit (httk)
#    Copyright (C) 2012-2025 The httk AUTHORS
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, sys, re, math
from decimal import Decimal
from fractions import Fraction
from collections import OrderedDict

class rewindable_iterator(object):

    def __init__(self, iterator):
        self._iter = iter(iterator)
        self._rewind = False
        self._cache = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._rewind:
            self._rewind = False
        else:
            self._cache = next(self._iter)
        return self._cache

    def rewind(self, rewindstr=None):
        if self._rewind:
            raise RuntimeError("Tried to backup more than one step.")
        elif self._cache is None:
            raise RuntimeError("Can't backup past the beginning.")
        self._rewind = True
        if rewindstr is not None:
            self._cache = rewindstr


def _read_cif_rewind_if_needed(f, row, done_fields):
    splitstr = row.lstrip().split(None, done_fields)
    if len(splitstr) > 1:
        rest = splitstr[-1]
        if rest.strip() != "":
            f.rewind(rest)
            return True
        return False
    else:
        return False


def _read_cif_loop(f, pragmatic=True, allow_cif2=False, use_types=False):
    noteol = False
    loop_data = OrderedDict()
    header = []
    for row in f:
        striprow = row.strip()
        lowrow = striprow.lower()
        if lowrow.startswith("_"):
            loop_data[lowrow[1:]] = []
            header += [lowrow[1:]]
            noteol = _read_cif_rewind_if_needed(f, row, 1)
        else:
            f.rewind()
            break

    while True and len(header)>0:
        for i in range(len(header)):
            try:
                row = next(f)
                while row.isspace():
                    row = next(f)
            except StopIteration:
                break
            striprow = row.strip()
            lowrow = striprow.lower()
            if not row or row.startswith("_") or lowrow.startswith("data_") or lowrow.startswith("loop_"):
                f.rewind()
                break
            f.rewind()
            val, noteol = _read_cif_data_value(f, noteol, pragmatic, allow_cif2, use_types, inloop=True)
            if val is None:
                # Could be a comment line, etc.
                continue
            loop_data[header[i]].append(val)
        else:
            continue
        break
    return loop_data


def _read_cif_data_value(f, noteol, pragmatic=True, allow_cif2=False, use_types=False, inloop=False, inlist=False):
    data_value = None
    for row in f:
        striprow = row.strip()
        if striprow.startswith("#") or striprow == "":
            noteol = False
            continue
        elif (not noteol) and row.startswith(';'):
            folded = False
            newline = False
            data_value = ""
            if row[1] == "\\" and row[2:].rstrip("\r\n") == "":
                folded = True
            elif row[1:].isspace():
                if not pragmatic:
                    data_value = row.lstrip().rstrip('\r\n')
                    newline = True
            else:
                data_value = row.lstrip()[1:].rstrip('\r\n')
                newline = True
            stripirow = ""
            for irow in f:
                stripirow = irow.strip()
                if irow.startswith(';'):
                    break
                if newline:
                    data_value += '\n'
                    newline = True
                if folded and irow.rstrip('\r\n').endswith("\\"):
                    data_value += irow.rstrip('\r\n').rstrip("\\")
                    newline = False
                else:
                    data_value += irow.rstrip('\r\n')
                    newline = True
            if len(stripirow) > 1:
                f.rewind(stripirow[1:])
                noteol = True
            else:
                noteol = False
            break
        elif striprow.startswith("'") or striprow.startswith('"'):
            # The cif quoting rules are ... weird. Quotes are "escaped" if they are not followed by whitespace.
            quote = striprow[0]
            starti = 1
            for chari in range(1, len(striprow)-1):
                if striprow[chari] == quote and str(striprow[chari+1]).isspace():
                    endi = chari
                    endq = chari+1
                    break
            else:
                if striprow[-1] != quote:
                    starti = 0
                    endi = len(striprow)
                    endq = len(striprow)
                else:
                    endi = len(striprow)-1
                    endq = len(striprow)
            data_value = striprow[starti:endi]
            if endq != len(striprow):
                f.rewind(striprow[endq:])
                noteol = True
            else:
                noteol = False
            break
        elif allow_cif2 and inlist and striprow.startswith("]"):
            # TODO: Is ] allowed without whitespace after? I need to check the spec
            splitstr = striprow.split("]", 1)
            if len(splitstr)>1 and len(splitstr[1])>0:
                f.rewind(splitstr[1])
                noteol = True
            data_value = None
            break
        elif allow_cif2 and striprow.startswith("["):
            if len(striprow)>1:
                f.rewind(striprow[1:])
                noteol = True
            data_value = []
            while True:
                innerval, noteol = _read_cif_data_value(f, noteol, pragmatic, allow_cif2, use_types, inloop=False, inlist=True)
                if innerval is None:
                    break
                data_value += [innerval]
            break
        elif allow_cif2 and inlist and ("]" in striprow):
            splitstr2 = striprow.split("]", 1)
            splitstr = splitstr2[0].split(None, 1)
            data_value = splitstr[0].strip()
            rightside = ""
            if len(splitstr) > 1:
                f.rewind(splitstr[1]+"]"+splitstr2[1])
            else:
                f.rewind("]"+splitstr2[1])
            noteol = True
            break
        else:
            if pragmatic and not inloop:
                # In pragmatic mode, if we are not in a loop and there is more than one data value
                # separated by whitespace, read all of it. This should always be ok to do, since
                # multiple data values in this situation would be an
                # error in the file otherwise, but if there is whitespace + underscore/data_/loop_ we parse that
                # as a new symbol, since otherwise we COULD misread valid files (with very weird formatting...).
                splitstr = re.split(r'\s+_|\s+data_|\s+loop_', striprow, maxsplit=1)
            else:
                splitstr = striprow.split(None, 1)
            # "Data on a line following a hash character `#' is considered to be a comment,
            # except if it is contained within a text string."
            data_value = splitstr[0].partition("#")[0].strip()
            rightside = ""
            if len(splitstr) > 1:
                rightside = splitstr[1].strip()
            if rightside != "":
                f.rewind(rightside)
                noteol = True
            else:
                noteol = False
            break
    if use_types:
        if _cif_is_int(data_value):
            data_value = cif_to_int(data_value)
        elif _cif_is_float(data_value):
            data_value = cif_to_float(data_value)

    return data_value, noteol


def _read_cif_data_block(f, pragmatic=True, allow_cif2=False, use_types=False):
    data_items = OrderedDict()
    loops = 0
    for row in f:
        striprow = row.strip()
        lowrow = striprow.lower()
        if striprow.startswith("#"):
            continue
        elif lowrow.startswith("data_"):
            f.rewind()
            return data_items
        elif lowrow.startswith("loop_"):
            _read_cif_rewind_if_needed(f, row, 1)
            loopdata = _read_cif_loop(f, pragmatic, allow_cif2, use_types)
            data_items['loop_'+str(loops)] = list(loopdata.keys())
            loops += 1
            data_items.update(loopdata)
        elif striprow.startswith(";"):
            # Multi-line string that we've failed to tie to a name, lets just skip it, maybe we should warn
            for irow in f:
                if irow.rstrip() == ";":
                    break
        elif striprow.startswith("_"):
            lowsplit = lowrow.split()
            data_name = lowsplit[0][1:]
            if len(lowsplit) > 1:
                noteol = True
                rightside = striprow.split(None, 1)[1].strip()
                f.rewind(rightside)
            else:
                noteol = False
            data_value, noteol = _read_cif_data_value(f, noteol, pragmatic, allow_cif2, use_types, inloop=False)
            data_items[data_name] = data_value
    return data_items


def read_cif(fs, pragmatic=True, allow_cif2=False, use_types=False):
    """
    Generic cif reader, given a filename / ioadapter it places all data in a python dictionary.

    It returns a tuple: (header, list)
    Where list are pairs of data blocks names and data blocks

    Each data block is a dictionary with tag_name:value

    For loops, value is another dictionary with format column_name:value

    The optional parameter pragmatic regulates handling of some counter-intuitive aspects of the cif specification, where
    the default pragmatic=True handles these features the way people usually use them, whereas pragmatic=False means
    to read the cif file precisely according to the spec. For example, in a multiline text field::

        ;
        some text
        ;

    Means the string '\\nsome text'. For this specific case pragmatic=True removes the leading newline.

    set use_types to True to convert things that look like floats and integers to those respective types
    """
    if isinstance(fs, str):
        fs = open(fs, "r", encoding="utf-8", errors="surrogateescape")
        f = rewindable_iterator(fs)
    else:
        f = rewindable_iterator(fs)
    try:
        header = ""
        datalist = []
        for row in f:
            if row.strip().startswith("#"):
                header += row
            else:
                f.rewind()
                break

        for row in f:
            lowrow = row.strip().lower()
            if lowrow.startswith("data_"):
                data_block_name = lowrow.partition('_')[2].split()[0].strip()
                _read_cif_rewind_if_needed(f, row, 1)
                data_block = _read_cif_data_block(f, pragmatic, allow_cif2, use_types)
                datalist += [(data_block_name, data_block)]
    finally:
        fs.close()
    return datalist, header

