from numpy import array
from pylab import plot, show, grid, legend, xlabel, ylabel
from msys_readdata import *

plt = 1

if plt==1:
    fin     = read_table("mix01_fine.res")
    res     = read_table("mix01_nod_12_0.res")
    fin_ux  = array(fin['ux'])
    fin_fxe = array(fin['fx'])
    ux      = array(res['ux'])
    fxe     = array(res['fx'])
    fxi     = array(res['fx_int'])
    plot   (fin_ux,fin_fxe,'k-', label='fine')
    plot   (ux,fxe,'ro', label='fe')
    plot   (ux,fxi,'b-', label='fi')
    xlabel ('ux')
    ylabel ('fx')
    legend (loc='best')
    grid   ()
    show   ()

if plt==2:
    res     = read_table("mix02_nod_7_0.res")
    ux      = array(res['ux'])
    fxe     = array(res['fx'])
    fxi     = array(res['fx_int'])
    plot   (ux,fxe,'ro', label='fe')
    plot   (ux,fxi,'b-', label='fi')
    xlabel ('ux')
    ylabel ('fx')
    legend (loc='best')
    grid   ()
    show   ()
