# Copyright (C) 2009 Dorival M Pedroso
# ------------------------------------
# Models

from os.path import basename
from numpy   import *

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
