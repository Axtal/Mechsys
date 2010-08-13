from numpy import array
from pylab import *
from msys_readdata import read_table
from msys_plotter import Plotter

plt = 0

if plt==0:
    dat = read_table("owen_hinton_01.res")
    u    = array(dat['u'])
    fint = array(dat['fint'])
    fext = array(dat['fext'])

    #d2  = read_table("owen_hinton_01_disp.res")
    #u2   = array(d2 ['u'])
    #fi2  = array(d2 ['fint'])
    #fe2  = array(d2 ['fext'])

    plot(-u,-fint,'b-',lw=4, marker='s', clip_on=False)
    #plot(-u2,-fi2,'r-',lw=4, marker='*', clip_on=False)

    plot(-u,-fext,'c-',lw=2, clip_on=False)
    #plot(-u2,-fe2,'m-',lw=2, clip_on=False)

    grid()
    show()

elif plt==1:
    # 8, 12, 14, 18
    P   = 18
    res = read_table("owen_hinton_02_P%d.res"%P)
    dat = read_table("owen_hinton_02_P%d.dat"%P)

    plot(res['r'],res['st'],'b-',lw=2)
    plot(res['r'],res['st'],'bo')
    plot(dat['r'],dat['st'],'r*')

    grid()
    show()

elif plt==2:
    res  = read_table("owen_hinton_02_n41.res")
    dat  = read_table("owen_hinton_02_pu.dat")
    da0  = read_table("jiang.dat")
    ud   = array(dat['u'])
    u    = array(res['ur'])
    fint = array(res['fr_int'])
    fext = array(res['fr_ext'])

    subplot(1,2,1)
    plot(u,fext,'r-',lw=2,label='Fext')
    plot(u,fext,'ro')
    plot(u,fint,'b-',lw=2,label='Fint')
    plot(u,fint,'b*')
    xlabel('u'); ylabel('force')
    legend(loc='upper left')
    grid()

    subplot(1,2,2)
    plot(u,res['P'],'r-',lw=2)
    plot(da0['u'],da0['p'],'gd')
    plot(ud/100,dat['p'],'bo')
    xlabel('u'); ylabel('P')
    grid()

    show()

elif plt==3:
    p          = Plotter()
    p.fc_ty    = ['VM']
    p.fc_np    = 40
    p.fc_phi   = 0
    p.fc_cu    = 12.0
    p.rst_phi  = False
    p.rst_circ = False
    p.plot ("owen_hinton_02_ele_4_-1.res", draw_fl=True, draw_ros=True)
    show()

elif plt==4:
    dat = read_table("owen_hinton_02_mesh.dat")
    for i in range(1,5):
        print (dat['x'][i]-dat['x'][i-1])/10.
