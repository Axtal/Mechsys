import optparse
from msys_fig import *

# input
op = optparse.OptionParser()
op.add_option('--tst',  '-t', dest='tst',  default='1', help='test number')
opts, args = op.parse_args()

if opts.tst=='0':
    dat = read_table("wood_lewis_nod_34.res")

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

elif opts.tst=='1':
    dat = read_table("owen_hinton_03_nod_50.res")
    t = array(dat['Time'])
    plot(1.e3*t,dat['ux'],'r-',lw=2)
    grid()
    show()

elif opts.tst=='2':
    nod = 2
    #dat = read_table("zienk_shiomi_01_nod_20.res")
    dat = read_table("zienk_shiomi_01_nod_%d.res"%nod)
    T = array(dat['Time'])
    plot(T,dat['pw'],'r-',lw=2)
    grid()
    show()

elif opts.tst=='3':
    #Time = 0.01
    Time = 1
    dat = read_table("zienk_shiomi_01_%06.3f.res"%Time)
    plot(dat['pw'],dat['y'],'r-',lw=2)
    grid()
    show()
