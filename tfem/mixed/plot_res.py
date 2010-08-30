from numpy import array
from pylab import plot, show, grid, legend, xlabel, ylabel
from msys_readdata import *

plt = 1

if plt==1:
    fin    = read_table("mix01_fine.res")
    res    = read_table("mix01_nod_12_0.res")
    fin_ux = array(fin['ux'])
    fin_fx = array(fin['fx'])
    ux     = array(res['ux'])
    fx     = array(res['fx'])
    plot   (fin_ux,fin_fx,'k-', label='fx:fine')
    plot   (ux,fx,'r-', marker='o', label='fx')
    xlabel ('ux')
    ylabel ('fx')
    legend (loc='best')
    grid   ()
    show   ()

if plt==2:
    res    = read_table("mix02_nod_7_0.res")
    ux     = array(res['ux'])
    fx     = array(res['fx'])
    plot   (ux,fx,'ro', label='fx')
    xlabel ('ux')
    ylabel ('fx')
    legend (loc='best')
    grid   ()
    show   ()
