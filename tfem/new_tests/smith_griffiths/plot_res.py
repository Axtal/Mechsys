from numpy import *
from pylab import *
from data_handler import *

dat = read_table("fig_11_04_nod_17.res")

plot(dat['Time'],dat['uy'],'r-',lw=2)
grid()
show()
