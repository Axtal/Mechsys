from numpy import *
from pylab import *
from msys_readdata import *

plt = 0

if plt==0:
    dat = read_table("ep_test.res")
    u    = array(dat['u'])
    fint = array(dat['fint'])
    fext = array(dat['fext'])

    plot(-u,-fext,'r-',lw=4,label='Fext')
    plot(-u,-fext,'rs')

    plot(-u,-fint,'k-',lw=2,label='Fint')
    plot(-u,-fint,'k.')

    xlabel('U')
    ylabel('F')

    legend(loc='upper left')
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
