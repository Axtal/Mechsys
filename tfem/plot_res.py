from plotter import *
from data_handler import *

test = 1

if test==1:
    p = Plotter()
    p.plot ("labtest_ele_0.res", div_by_p=False)
    p.show()

if test==2:
    dat = read_table("labtest_nod_6.res")
    u   = array(dat['uz'])
    f   = array(dat['fz'])
    plot (-u,-f,'r-',lw=2)
    plot (-u,-f,'ro',lw=2)
    grid ()
    show ()
