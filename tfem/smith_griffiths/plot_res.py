from numpy import *
from pylab import *
from msys_readdata import *

if True:
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

if False:
    dat = read_table("fig_11_04_nod_17_-100.res")

    plot(dat['Time'],dat['uy'],'r-',lw=2)
    plot(dat['Time'],dat['uy'],'ro',lw=2)
    grid()
    show()
