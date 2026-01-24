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

_cif_ordinary_char = r'!%&()*+,-./0123456789:<=>?@ABCDEFGHIHJKLMNOPQRSTUVWXYZ\^`abcdefghijklmnopqrstuvwxyz{|}~'
_cif_non_blank_char = _cif_ordinary_char+'"'+"#$"+"'"+"_"+";[]"
_cif_text_lead_char = _cif_ordinary_char+'"'+"#$"+"'"+"_ \t[]"
_cif_any_print_char = _cif_ordinary_char+'"'+"#$"+"'"+"_ \t;[]"
_cif_helper_table = str.maketrans('', '')
_cif_non_blank_char_table = str.maketrans(_cif_non_blank_char, _cif_non_blank_char)
_cif_unicode_translation_table = {}
for i in range(sys.maxunicode+1):
    _cif_unicode_translation_table[i] = None
for key, value in _cif_non_blank_char_table.items():
    _cif_unicode_translation_table[key] = value
_cif_integer_regex = re.compile(r'^[+-]?[0-9]+$')
_cif_float_regex = re.compile(r'^[+-]?[0-9]+[eE][+-]?[0-9]+|([+-]?[0-9]*\.[0-9]+|[+-]?[0-9]\.)([eE][+-]?[0-9]+)?$')
_cif_simplestring_regex = re.compile(r'^[A-Za-z0-9()][A-Za-z0-9()+-]*$')

def _cif_validate_name(name_unfiltered, context=None):
    if context is not None:
        context = context+": "+name_unfiltered
    name = _cif_validate_non_blank_char(name_unfiltered, context)
    if len(name) > 75:
        sys.stderr.write("***Warning: write_cif: name length > 75, surplus characters removed in "+context+": "+name_unfiltered)
        name = name[:75]
    return name

def _cif_is_float(data_value):
    return (_cif_float_regex.match(data_value) is not None)

def _cif_is_simplestring(data_value):
    return (_cif_simplestring_regex.match(data_value) is not None)

def _cif_is_int(data_value):
    return (_cif_integer_regex.match(data_value) is not None)

def _cif_validate_non_blank_char(s, context=None):
    if sys.version_info[0] == 3:
        out = s.translate(_cif_unicode_translation_table)
    else:
        out = s.translate(_cif_helper_table, _cif_non_blank_char_table)
    if out != s:
        if context is not None:
            sys.stderr.write("***Warning: write_cif: non-permitted characters in "+context+" removed.")
        else:
            sys.stderr.write("***Warning: write_cif: non-permitted characters removed.")
    return out

def _cif_write_semicolontextfield(f, lines, noteol, max_line_length):
    if noteol:
        f.write("\n")
        noteol = False
    for i in range(len(lines)):
        lines[i] = lines[i].rstrip("\r\n")
        if lines[i][0] == ';':
            sys.stderr.write("***Warning: write_cif: had to insert space before semicolon at the start of a line of a multi-line string to fulfill arcane quoting rules.")
            lines[i] = ' '+lines[i]
        if len(lines[i]) > max_line_length:
            f.write(";\\"+"\n")
            break
    else:
        f.write(";")
    for line in lines:
        if len(line) > max_line_length:
            sublines = [line[i:i+max_line_length-2] for i in range(0, len(line), max_line_length-2)]
            # Handle a wonderful corner case: the line splitting for length creates lines that start with one, or more, semi-colons..., sigh...
            for i in range(1, len(sublines)):
                if sublines[i][0] == ";":
                    if len(sublines[i]) > 1 and sublines[i][1] != ";":
                        # If its just a single semi-colon, move it to the previous line, which we saved space for by splitting at max_line_length-2
                        sublines[i-1] += ";"
                        sublines[i] = sublines[i][1:]
                    else:
                        # Multiple semi-colons in a row, or a semi-colon + newline, this is a possibly unresolvable case (think long string of only semi-colons)
                        # fudge a solution by inserting a space
                        sys.stderr.write("***Warning: write_cif: had to insert space before semicolon in a long string to fulfill arcane quoting rules.")
                        sublines[i] = " "+sublines[i]
            for subline in sublines:
                f.write(subline+"\\"+"\n")
        else:
            f.write(line+"\n")

    f.write(";\n")
    return False


def _cif_write_data_value(f, orig_data_value, noteol, max_line_length, use_types, inloop):
    if orig_data_value is None:
        data_value = ""
    else:
        data_value = str(orig_data_value)
    has_whitespace = len(data_value.split()) > 1
    lines = data_value.splitlines()
    has_lines = len(lines) > 1
    has_single_quote = data_value.find("'") != -1
    has_double_quote = data_value.find('"') != -1
    too_long = len(data_value) + 2 > max_line_length
    if has_lines or (has_single_quote and has_double_quote) or too_long:
        noteol = _cif_write_semicolontextfield(f, lines, noteol, max_line_length)
        return noteol
    elif has_double_quote or (has_whitespace and not has_single_quote) or data_value == "":
        f.write("'"+data_value+"'")
        return True
    elif has_single_quote or (has_whitespace and not has_double_quote):
        f.write('"'+data_value+'"')
        return True
    elif not use_types:
        # Skip quotes if it looks like a number or is a simple string used in a loop
        if _cif_is_float(data_value):
            f.write(data_value)
            return True
        elif _cif_is_int(data_value):
            f.write(data_value)
            return True
        elif inloop and _cif_is_simplestring(data_value):
            f.write(data_value)
            return True
        else:
            f.write("'"+data_value+"'")
            return True
    else:
        # Always quote when a string, never quote otherwise
        if isinstance(orig_data_value, str):
            f.write("'"+data_value+"'")
        else:
            f.write(data_value)
        return True


def is_sequence(l):
    return isinstance(l, Iterable) and not isinstance(l, str)

def _cif_write_data_block(f, data_block, max_line_length, use_types):
    for key in data_block:
        val = data_block[key]
        if key.startswith("loop_"):
            f.write("loop_\n")
            outdata_columns = []
            for unfiltered_column in val:
                column = _cif_validate_non_blank_char(unfiltered_column, "column name: "+unfiltered_column)
                f.write("_"+column+"\n")
                outdata_columns += [data_block[unfiltered_column]]
            if len(outdata_columns) > 0:
                noteol = False
                for i in range(len(outdata_columns[0])):
                    column_count = 0
                    for j in range(len(outdata_columns)):
                        column_count += len(str(outdata_columns[j][i]))+2
                        if column_count > max_line_length and noteol:
                            f.write("\n")
                            column_count = 0
                            noteol = False
                        noteol = _cif_write_data_value(f, outdata_columns[j][i], noteol, max_line_length, use_types, inloop=True)
                        if noteol:
                            f.write(" ")
                            column_count += 1
                        else:
                            column_count = 0
                    if noteol:
                        noteol = False
                        f.write("\n")
        elif is_sequence(val):
            continue
        else:
            data_name = _cif_validate_name(key)
            # Do we have space _ + key + space + quote + the whole data value + quote?, if not, preemptively break line
            f.write("_"+data_name+" ")
            if len(data_name)+len(str(val))+4 > max_line_length:
                f.write("\n")
                noteol = False
            else:
                noteol = True
            noteol = _cif_write_data_value(f, val, noteol, max_line_length, use_types, inloop=False)
            if noteol:
                f.write("\n")
                noteol = False


def write_cif(ioa, data, header=None, max_line_length=80, use_types=False):
    """
    Generic cif writer, given a filename / ioadapter

    data = the cif data to write as an (ordered) dictionary of tag_name:value

    header = the header (comment) segment

    max_line_length = the maximum number of characters allowed on each line. This should not be set < 80
    (there is no point, and the length calculating algorithm breaks down at some small line length)

    use_types =

       if True: always quote values that are of string type. Numeric values are put in the file unquoted (as they should)
       if False (default): also strings that look like cif numbers are put in the file unquoted

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

    ioa = IoAdapterFileWriter.use(ioa)
    f = ioa.file

    if header is not None:
        lines = header.splitlines()
        for line in lines:
            if len(line) > max_line_length:
                header = "#\n" + header
                break
        for line in lines:
            if len(line) > max_line_length:
                sublines = [line[i:i+79] for i in range(0, len(line), 79)]
                for subline in sublines:
                    f.write(subline+"\\"+"\n")
            else:
                f.write(line+"\n")

    data_block_count = -1
    for data_block in data:
        data_block_count += 1
        data_block_name_unfiltered = data_block[0]
        if data_block_name_unfiltered is None:
            data_block_name = "data_"+str(data_block_count)
        else:
            data_block_name = _cif_validate_name(data_block_name_unfiltered, "data block name")
            if data_block_name == "":
                data_block_name = "data_"+str(data_block_count)

        f.write("data_"+data_block_name+"\n")
        _cif_write_data_block(f, data_block[1], max_line_length, use_types)
    ioa.close()
