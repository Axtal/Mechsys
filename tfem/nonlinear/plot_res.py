from numpy import array
from pylab import *
from msys_readdata import read_table
from msys_plotter import Plotter

if False:
    dat = read_table("owen_hinton_01.res")
    u    = array(dat['u'])
    fint = array(dat['fint'])
    fext = array(dat['fext'])

    plot(-u,-fext,'r-',lw=2)
    plot(-u,-fext,'ro')

    plot(-u,-fint,'b-',lw=2)
    plot(-u,-fint,'b*')

    grid()
    show()

if False:
    P   = 18
    res = read_table("owen_hinton_02_P%d.res"%P)
    dat = read_table("owen_hinton_02_P%d.dat"%P)

    plot(res['r'],res['st'],'b-',lw=2)
    plot(res['r'],res['st'],'bo')
    plot(dat['r'],dat['st'],'r*')

    grid()
    show()

if False:
    dat = read_table("owen_hinton_02_mesh.dat")
    for i in range(1,5):
        print (dat['x'][i]-dat['x'][i-1])/10.
    #plot(dat['x'],dat['y'],'bo')
    #grid()
    #show()

if True:
    res  = read_table("owen_hinton_02_n41.res")
    dat  = read_table("owen_hinton_02_pu.dat")
    da0  = read_table("jiang.dat")
    ud   = array(dat['u'])
    u    = array(res['ur'])
    fint = array(res['fr_int'])
    fext = array(res['fr_ext'])

    #subplot(1,2,1)
    #plot(u,fext,'r-',lw=2)
    #plot(u,fext,'ro')
    #plot(u,fint,'b-',lw=2)
    #plot(u,fint,'b*')
    #xlabel('u'); ylabel('force')
    #grid()

    #subplot(1,2,2)
    plot(u,res['P'],'r-',lw=2)
    plot(da0['u'],da0['p'],'gd')
    plot(ud/100,dat['p'],'bo')
    xlabel('u'); ylabel('P')
    grid()

    show()

if False:
    p        = Plotter()
    p.fc_cu  = 12.0
    p.fc_c   = 12.0
    p.fc_phi = 0.1
    p.fc_p   = 0.0
    p.fc_np  = 40
    p.plot ("owen_hinton_02_ele_4.res", draw_fl=True)
    show()
