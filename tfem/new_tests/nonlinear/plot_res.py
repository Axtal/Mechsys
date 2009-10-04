from numpy import *
from pylab import *
from data_handler import *

test = 2

if test==1:
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

if test==2:
    P   = 12
    res = read_table("owen_hinton_02_P%d.res"%P)
    dat = read_table("owen_hinton_02_P%d.dat"%P)

    plot(res['r'],res['st'],'b-',lw=2)
    plot(res['r'],res['st'],'bo')
    plot(dat['r'],dat['st'],'r*')

    grid()
    show()

if test==-2:
    dat = read_table("owen_hinton_02_mesh.dat")
    for i in range(1,5):
        print (dat['x'][i]-dat['x'][i-1])/10.
    #plot(dat['x'],dat['y'],'bo')
    #grid()
    #show()
