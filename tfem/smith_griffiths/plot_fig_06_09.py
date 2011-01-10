from msys_fig import *
from msys_plt import *

plt = 1

if plt==0:
    dat = read_table("fig_06_09.res")
    u    = array(dat['u'])
    fint = array(dat['fint'])
    fext = array(dat['fext'])

    plot(-fint,u,'b-',lw=4, clip_on=False)
    plot(-fint,u,'bs',      clip_on=False)

    plot(-fext,u,'r-',lw=2, clip_on=False)
    plot(-fext,u,'r.',      clip_on=False)

    grid()
    show()

if plt==1:
    cu    = 100.0
    DelP  = 520.0
    d0    = read_table("fig_06_09_FE.res")
    d1    = read_table("fig_06_09_ME.res")
    d2    = read_table("fig_06_09_NR.res")
    d3    = read_table("sg_fig611.dat")
    d4    = read_table("sg_fig611fine.dat")
    d5    = read_table("sg_fig611veryfine.dat")
    u0    = array(d0['u'])
    u1    = array(d1['u'])
    u2    = array(d2['u'])
    u3    = array(d3['u'])
    u4    = array(d4['u'])
    u5    = array(d5['u'])
    q_cu0 = DelP*array(d0['Time'])/cu
    q_cu1 = DelP*array(d1['Time'])/cu
    q_cu2 = DelP*array(d2['Time'])/cu
    q_cu3 = array(d3['q'])/cu
    q_cu4 = array(d4['q'])/cu
    q_cu5 = array(d5['q'])/cu

    plot(q_cu1,u1,'r-',lw=6, marker='s', ms=15, clip_on=False, label='MechSys (ME)')
    plot(q_cu2,u2,'k-',lw=3, marker='o', ms=12, clip_on=False, label='MechSys (NR)')
    plot(q_cu0,u0,'m-',lw=2, marker='+', ms=10, clip_on=False, label='MechSys (FE)')
    plot(q_cu3,u3,'b-',lw=1, marker='.', ms=8,  clip_on=False, label='Smith & Griffiths')
    #plot(q_cu4,u4,'y-',lw=2, marker='d', ms=8,  clip_on=False, label='Smith & Griffiths (fine)')
    plot(q_cu5,u5,'g-',lw=2            , ms=9,  clip_on=False, label='Smith & Griffiths (very fine)')
    xlabel (r'$q/c_u$');
    ylabel (r'$\delta$');

    axvline (5.14,color='green',label='Prandtl (q/cu=5.14)')

    legend(loc='lower left')
    grid()
    show()
