from numpy import *
from pylab import *
from data_handler import *

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
