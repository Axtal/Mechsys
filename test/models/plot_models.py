from msys_plotter    import *
from msys_invariants import *
from msys_readdata   import *

if False:
    p = Plotter()
    p.plot ("test_models.res", draw_ys=False,draw_fl=False,dpt_out=5,pqty='cam')
    p.plot ('FO1_CTR_01.kgf.pct.dat',fem_res=False)
    p.show()

else:
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi  = M_calc_phi(1,'cam')
    p.fc_poct = 150.0*sqrt(3.0)
    p.show_k  = True
    p.plot ("test_models_2.res",clr='red',   marker='+', draw_fl=True,draw_ros=True)
    p.plot ("test_models_3.res",clr='green', marker='+')
    p.plot ("test_models_4.res",clr='blue',  marker='+')

    # plot data
    dat = read_tables(['../../tfem/mdl_tst_01.dat','../../tfem/mdl_tst_02.dat'])
    for i, ed  in enumerate(dat['mdl_tst_01']['ed']):  dat['mdl_tst_01']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, ed  in enumerate(dat['mdl_tst_02']['ed']):  dat['mdl_tst_02']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, mev in enumerate(dat['mdl_tst_02']['mev']): dat['mdl_tst_02']['mev'][i] *= 100.0
    subplot(2,3,1); plot(dat['mdl_tst_01']['ed'],dat['mdl_tst_01']['q'],  'ko')
    subplot(2,3,4); plot(dat['mdl_tst_02']['ed'],dat['mdl_tst_02']['mev'],'ko')

    p.show()
