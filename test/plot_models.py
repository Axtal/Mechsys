from msys_plotter  import *
from msys_readdata import *

tst = 1

if tst==1 or tst==2:
    p = Plotter()
    p.show_k = True
    p.plot ("test_models.res")

    # plot data
    dat = read_tables(['../tfem/mdl_tst_01.dat','../tfem/mdl_tst_02.dat'])
    for i, ed  in enumerate(dat['mdl_tst_01']['ed']):  dat['mdl_tst_01']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, ed  in enumerate(dat['mdl_tst_02']['ed']):  dat['mdl_tst_02']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, mev in enumerate(dat['mdl_tst_02']['mev']): dat['mdl_tst_02']['mev'][i] *= 100.0
    subplot(2,3,1); plot(dat['mdl_tst_01']['ed'],dat['mdl_tst_01']['q'],  'ko')
    subplot(2,3,4); plot(dat['mdl_tst_02']['ed'],dat['mdl_tst_02']['mev'],'ko')

    p.show()

if tst==3:
    p = Plotter()
    p.plot ("test_models.res", draw_ys=False,draw_fl=False,dpt_out=5,pqty='cam')
    p.plot ('FO1_CTR_01.kgf.pct.dat',fem_res=False)
    p.show()
