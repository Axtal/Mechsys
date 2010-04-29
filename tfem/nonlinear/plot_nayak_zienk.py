from numpy import *
from pylab import *
from msys_readdata import *

plt = 1

if plt==0:
    dat = read_table("nayak_zienk_01.res")
    u    = array(dat['u'])
    fint = array(dat['fint'])
    fext = array(dat['fext'])

    plot(-u,-fext,'r-',lw=2)
    plot(-u,-fext,'ro')

    plot(-u,-fint,'b-',lw=2)
    plot(-u,-fint,'b*')

    grid()
    show()

elif plt==1:

    d1   = read_table("nayak_fig5a_ua.dat")
    u_a1 = array(d1['u_a'])
    p_y1 = array(d1['p_y'])


    r = read_table("nayak_zienk_01_nod_2_-200.res")
    a = 3.0
    u_a2 = 100.0*array(r['ux'])/a
    p_y2 = array(r['Time'])*1.4

    plot (u_a1, p_y1, 'r-', marker='s')
    plot (u_a2, p_y2, 'b-', lw=2)

    grid()
    show()
