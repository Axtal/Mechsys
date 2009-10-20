from numpy import array
from pylab import subplot, plot, grid, show
from msys_readdata import read_table
from msys_plotter import Plotter

if False:
    #dat = read_table("owen_hinton_nod_2.res")
    #uy = array(dat['uy'])
    #fy = array(dat['fy'])

    #plot(-uy,-fy,'r-',lw=2)
    #plot(-uy,-fy,'bo',lw=2)


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

if False:
    dat = read_table("owen_hinton_02_n41.res")
    u    = array(dat['ur'])
    fint = array(dat['fr_int'])
    fext = array(dat['fr_ext'])

    subplot(1,2,1)
    plot(u,fext,'r-',lw=2)
    plot(u,fext,'ro')
    plot(u,fint,'b-',lw=2)
    plot(u,fint,'b*')
    grid()

    subplot(1,2,2)
    plot(u,dat['P'],'r-',lw=2)
    grid()

    show()

if True:
    p        = Plotter()
    p.fc_c   = 12.0
    p.fc_phi = 0.5
    p.fc_p   = 0.0
    p.fc_np  = 40
    p.plot ("owen_hinton_02_ele_4.res", draw_fl=True)
    show()
