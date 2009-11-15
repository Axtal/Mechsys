from numpy import linspace, pi, exp, sin, array
from pylab import plot, grid, show
from msys_readdata import *

if False:
    dat = read_table("wood_lewis_nod_34_0.res")

    def solution(x,t):
        H = 1.0
        for k in range(1,11):
            c  = pi*(2.0*k-1.0)
            H -= (4.0/c)*exp(-t*(c/8.0)**2.0)*sin(c*x/8.0)
        return H

    T = linspace(0.0,30.0,100)
    H = solution(1.0,T)

    plot(T,H,'b-',lw=2)
    plot(dat['Time'],dat['H'],'r-',lw=2)
    grid()
    show()

if False:
    dat = read_table("owen_hinton_03_nod_50_-200.res")
    t = array(dat['Time'])
    plot(1.e3*t,dat['ux'],'r-',lw=2)
    grid()
    show()

if True:
    #dat = read_table("zienk_shiomi_01_nod_20_-200.res")
    dat = read_table("zienk_shiomi_01_nod_2_0.res")
    #dat = read_table("zienk_shiomi_01_nod_10_0.res")
    #dat = read_table("zienk_shiomi_01_nod_18_0.res")
    T = array(dat['Time'])
    plot(T,dat['pw'],'r-',lw=2)
    grid()
    show()
