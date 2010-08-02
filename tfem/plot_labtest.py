from pylab import legend, show
from msys_plotter  import *
from msys_readdata import *
from msys_invariants import M_calc_phi

plt = 0

if plt==0:
    phi      = M_calc_phi(1.0,'cam')
    Moct     = phi_calc_M(phi,'oct')
    pf       = 150.0*sqrt(3.0)
    qf       = Moct*pf
    cu       = qf*sqrt(3.0)/(2.0*sqrt(2.0))
    p        = Plotter()
    p.show_k = True
    p.fc_p   = pf
    p.fc_phi = phi
    p.fc_cu  = cu
    p.fc_c   = 0.0
    p.plot ("labtest_ele_0_-1.res", draw_fl=True, draw_ros=True)

    # plot data
    dat = read_tables(['mdl_tst_01.dat','mdl_tst_02.dat'])
    for i, ed  in enumerate(dat['mdl_tst_01']['ed']):  dat['mdl_tst_01']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, ed  in enumerate(dat['mdl_tst_02']['ed']):  dat['mdl_tst_02']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, mev in enumerate(dat['mdl_tst_02']['mev']): dat['mdl_tst_02']['mev'][i] *= 100.0
    subplot(2,3,1); plot(dat['mdl_tst_01']['ed'],dat['mdl_tst_01']['q'],  'ko')
    subplot(2,3,4); plot(dat['mdl_tst_02']['ed'],dat['mdl_tst_02']['mev'],'ko')
    show()

elif plt==1:
    dat  = read_table("labtest_uf.res")
    u    = array(dat['u'])
    fint = array(dat['f_int'])
    fext = array(dat['f_ext'])
    plot(-u,-fext,'r-',lw=2)
    plot(-u,-fext,'ro')
    plot(-u,-fint,'b-',lw=2)
    plot(-u,-fint,'b*')
    legend(['f_ext','f_ext','f_int','f_int'], loc='lower right')
    grid()
    show()

elif plt==2:
    dat = read_table("labtest_nod_7_-7.res")
    u   = array(dat['uz'])
    f   = array(dat['fz'])
    plot (-u,-f,'r-',lw=2)
    plot (-u,-f,'ro',lw=2)
    grid ()
    show ()
