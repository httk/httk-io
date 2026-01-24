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

def compose_ops_with_centerings(ops, centerings):
    """
    ops:         list of (R, t, time_flag) from _space_group_symop_magn_operation
    centerings:  list of (c, time_c) where c is 3-vector fractional translation,
                 time_c is 0 or 1 (0 = no time reversal)

    Returns a new list of (R, t', time') where t' = t + c and time' = (time_flag + time_c)%2.
    """
    composed = []
    for R, t, time_flag in ops:
        for Rc, c, time_c in centerings:
            if Rc != ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                raise Exception("Centering symop that includes rotation is invalid")
            t_new = (t[0] + c[0], t[1] + c[1], t[2] + c[2])
            time_new = (time_flag + time_c)%2 # time_flag * time_c for -1/+1 convention
            composed.append((R, t_new, time_new))
    return composed
