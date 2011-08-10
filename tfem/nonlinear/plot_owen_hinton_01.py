import optparse
from msys_fig import *

# input
op = optparse.OptionParser()
op.add_option('--tst',  '-t', dest='tst',  default='0', help='test number')
opts, args = op.parse_args()

if opts.tst=='0':
    res = read_table("owen_hinton_01.res")
    u    = array(res['u'])
    fint = array(res['fint'])
    fext = array(res['fext'])

    res = read_table("owen_hinton_p80_fig3.6.dat")
    udat = array(res['u'])
    Pdat = array(res['P'])

    #d2  = read_table("owen_hinton_01_disp.res")
    #u2   = array(d2 ['u'])
    #fi2  = array(d2 ['fint'])
    #fe2  = array(d2 ['fext'])

    plot(-u, -fint,'b-',lw=4, marker='s', clip_on=False)
    plot(udat,Pdat,'r-',lw=1, marker='*', clip_on=False)
    #plot(-u2,-fi2,'r-',lw=4, marker='*', clip_on=False)

    plot(-u,-fext,'c-',lw=2, clip_on=False)
    #plot(-u2,-fe2,'m-',lw=2, clip_on=False)

    xlabel('u')
    ylabel('P')
    grid()
    show()
