from plotter import *
from data_handler import *

test = 1

if test==1:
    p = Plotter()
    p.plot ("labtest_ele_0.res", div_by_p=False)

    # plot data
    dat = read_tables(['mdl_tst_01.dat','mdl_tst_02.dat'])
    for i, ed  in enumerate(dat['mdl_tst_01']['ed']):  dat['mdl_tst_01']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, ed  in enumerate(dat['mdl_tst_02']['ed']):  dat['mdl_tst_02']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, mev in enumerate(dat['mdl_tst_02']['mev']): dat['mdl_tst_02']['mev'][i] *= 100.0
    subplot(3,3,1); plot(dat['mdl_tst_01']['ed'],dat['mdl_tst_01']['q'],  'ko')
    subplot(3,3,4); plot(dat['mdl_tst_02']['ed'],dat['mdl_tst_02']['mev'],'ko')

    p.show()

if test==2:
    dat = read_table("labtest_nod_6.res")
    u   = array(dat['uz'])
    f   = array(dat['fz'])
    plot (-u,-f,'r-',lw=2)
    plot (-u,-f,'ro',lw=2)
    grid ()
    show ()
