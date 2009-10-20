########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Dorival M. Pedroso                                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

from os.path import basename

# Read file with table
# ====================
# dat: dictionary with the following content:
#   dat = {'sx':[1,2,3],'ex':[0,1,2]}
def read_table(filename):
    file   = open(filename,'r')
    header = file.readline().split()
    dat    = {}
    for key in header: dat[key] = []
    for lin in file:
        res = lin.split()
        for i, key in enumerate(header):
            dat[key].append(float(res[i]))
    file.close()
    return dat

# Read many files with tables
# ===========================
# filenames: list with file names
# dat: dictionary with the following content:
#   dat = {'fkey1':{'sx':[1,2,3],'ex':[0,1,2]}}
def read_tables(filenames):
    dat = {}
    for fn in filenames:
        fkey      = basename(fn).replace('.dat','')
        file      = open(fn,'r')
        header    = file.readline().split()
        dat[fkey] = {}
        for key in header: dat[fkey][key] = []
        for lin in file:
            res = lin.split()
            for i, key in enumerate(header):
                dat[fkey][key].append(float(res[i]))
        file.close()
    return dat
