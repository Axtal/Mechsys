from numpy import *
from pylab import *
from msys_readdata import *

plt = 1

if plt==0:
    dat = read_table("nayak_zienk_01.res")
    u    = array(dat['u'])
    fint = array(dat['fint'])
    fext = array(dat['fext'])

    plot(u,fext,'r-',lw=2,label='Fext')
    plot(u,fext,'ro')

    plot(u,fint,'b-',lw=2,label='Fint')
    plot(u,fint,'b*')

    xlabel('U')
    ylabel('F')

    legend(loc='upper left')
    grid()
    show()

elif plt==1:

    d1   = read_table("nayak_fig5a_ua.dat")
    d2   = read_table("nayak_fig5a_ub.dat")
    u_a1 = array(d1['u_a'])
    p_y1 = array(d1['p_y'])
    u_b1 = array(d2['u_b'])
    p_y3 = array(d2['p_y'])

    r1 = read_table("nayak_zienk_01_nod_2_-200.res")
    r2 = read_table("nayak_zienk_01_nod_17_-200.res")
    a  = 3.0
    u_a2 = 100.0*array(r1['ux'])/a
    p_y2 = array(r1['Time'])*1.392
    u_b2 = 100.0*array(r2['ux'])/a
    p_y4 = array(r2['Time'])*1.392

    plot (u_b1, p_y3, 'r-', marker='s',lw=3, label='Nayak & Zienkiewicz (1972)')
    plot (u_a1, p_y1, 'r-', marker='s',lw=3, label='Nayak & Zienkiewicz (1972)')
    plot (u_a2, p_y2, 'k-', marker='.',      label='MechSys')
    plot (u_b2, p_y4, 'k-', marker='.',      label='MechSys')

    legend(loc='lower right')
    xlabel('ux')
    ylabel('p')
    grid()
    show()
