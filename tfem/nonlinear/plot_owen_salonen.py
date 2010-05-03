from numpy import array
from pylab import *
from msys_readdata import read_table
from msys_plotter  import Plotter
from msys_linfit   import LinFit

plt = 1

if plt==0:
    # data
    d1  = read_table("owen_salonen_fig5.dat")
    d2  = read_table("owen_salonen_fig5_bilinear.dat")
    Ea1 = array(d1['Ea'])
    Sa1 = array(d1['Sa'])
    Ea2 = array(d2['Ea'])
    Sa2 = array(d2['Sa'])

    i = 11
    l1 = LinFit (Ea2[:i+1],Sa2[:i+1], tls=True, cmx=False)
    l2 = LinFit (Ea2[i:],  Sa2[i:],   tls=True, cmx=True)
    X1 = linspace(0.,Ea2[i], 100)
    X2 = linspace(0.,Ea2[-1],100)

    #print "L1: c =", l1.c, "  m =", l1.m
    #print "L2: c =", l2.c, "  m =", l2.m

    sY = Sa2[i]
    E  = l1.m
    Et = l2.m
    Hp = 1.0/(1.0/Et-1.0/E)
    print "sY =", sY, "  E =", E, "  Et =", Et, "  Hp =", Hp

    plot(Ea1,Sa1,'r-',marker='o',label='Data')
    plot(Ea2,Sa2,'b-',lw=2,label='Owen and Salonen (1975)')
    plot(Ea2[i],Sa2[i],'y*')
    plot(X1,l1.y(X1),'g-')
    plot(X2,l2.y(X2),'g-')

    # simulation
    r1 = read_table("owen_salonen_uni_ele_0_-1.res")
    ea = array(r1['ez'])
    sa = array(r1['sz'])

    plot (ea,sa,'k-',marker='+',label='MechSys')

    legend(loc='lower right')
    grid()
    show()

elif plt==1:
    d1 = read_table("owen_salonen_fig6_pstrain.dat")
    p  = array(d1['p'])*2.0#2.0=>thickness = 0.5
    plot (d1['u'],p,'ro',label='Owen and Salonen (1975)')

    delp = 18000.0

    r1 = read_table("owen_salonen_3d_nod_0_-1.res")
    r2 = read_table("owen_salonen_3d_nod_2_-2.res")
    r3 = read_table("owen_salonen_3d_nod_81_-5.res")
    p1 = array(r1['Time'])*delp
    #plot (r1['uy'],p1,'b-',lw=3,marker='s',clip_on=False,label='MechSys (3D)')
    #plot (r2['uy'],p1,'r-',lw=2,marker='+',clip_on=False)
    #plot (r3['uz'],p1,'b-',lw=2,marker='.',clip_on=False)

    s1 = read_table("owen_salonen_2d_nod_0_0.res")
    s2 = read_table("owen_salonen_2d_nod_27_0.res")
    p2 = array(s1['Time'])*delp
    plot (s1['ux'],p2,'k-',marker='.',clip_on=False,label='MechSys (2D)')
    #plot (s2['uy'],p2,'r-',marker='d',clip_on=False)

    legend(loc='lower right')
    grid()
    show()
