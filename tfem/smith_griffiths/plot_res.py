from numpy import *
from pylab import *
from msys_readdata import *

if False:
    dat = read_table("fig_11_01_nod_1_-200.res")

    tsw = 1.0
    def calc_U(t):
        if (t<tsw): return  0.441*sin(pi*t/tsw)-0.216*sin(2.0*pi*t/tsw);
        else:       return -0.432*sin(6.284*(t-tsw))

    T = linspace(0.0,1.8,100)
    U = zeros(len(T))
    for i in range(len(T)): U[i] = calc_U(T[i])

    plot(T,U,'b-')
    plot(dat['Time'],dat['uy'],'ro',lw=2)
    grid()
    show()

if True:
    res   = read_table("fig_11_04_nod_17_-100.res")
    p113d = read_table("sg_11_07_p113.dat")
    p113  = read_table("sg_11_07_p113.res")
    p114  = read_table("sg_11_07_p114.res")

    plot(res  ['Time'],res  ['uy'],'r-',lw=2)
    plot(p113d['Time'],p113d['uy'],'bo')
    plot(p113 ['Time'],p113 ['uy'],'b-')
    #plot(p114 ['Time'],p114 ['uy'],'b-',lw=1)
    #plot(p114 ['Time'],p114 ['uy'],'b+',lw=1)

    #legend(['MechSys','SG:p113d','SG:p113','SG:p114'],loc='upper left')
    legend(['MechSys','SG:p113d','SG:p113'],loc='upper left')
    grid()
    show()
