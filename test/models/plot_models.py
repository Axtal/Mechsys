from msys_plotter    import *
from msys_invariants import *
from msys_readdata   import *
from numpy import array, log, sqrt
from pylab import plot, show, grid, legend, xlabel, ylabel

tst = 3

if tst==0:
    p = Plotter()
    p.plot ("test_models.res", draw_ys=False,draw_fl=False,dpt_out=5,pqty='cam')
    p.plot ('FO1_CTR_01.kgf.pct.dat',fem_res=False)
    p.show()

if tst==1:
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi   = M_calc_phi(1,'cam')
    p.fc_poct  = 150.0*sqrt(3.0)
    p.show_k   = True
    p.oct_sxyz = True
    p.fsz      = 14
    p.lwd=2; p.plot ("driver.res",clr='blue', markevery=10, label='CCM', draw_ros=True)
    legend()

    # plot data
    dat = read_tables(['../../tfem/mdl_tst_01.dat','../../tfem/mdl_tst_02.dat'])
    for i, ed  in enumerate(dat['mdl_tst_01']['ed']):  dat['mdl_tst_01']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, ed  in enumerate(dat['mdl_tst_02']['ed']):  dat['mdl_tst_02']['ed'] [i] *= (sqrt(3.0/2.0)*100.0)
    for i, mev in enumerate(dat['mdl_tst_02']['mev']): dat['mdl_tst_02']['mev'][i] *= 100.0
    subplot(2,3,1); plot(dat['mdl_tst_01']['ed'],dat['mdl_tst_01']['q'],  'ko')
    subplot(2,3,4); plot(dat['mdl_tst_02']['ed'],dat['mdl_tst_02']['mev'],'ko')

    p.show()

if tst==2:
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi  = M_calc_phi(1,'cam')
    p.fc_poct = 150.0*sqrt(3.0)
    p.show_k  = True
    #p.plot ("test_models_2.res",clr='red',   marker='+', draw_fl=True,draw_ros=True)
    p.plot ("test_models_1.res",clr='green', marker='+', label='Elastic')
    p.lwd=1; p.plot ("test_models_4.res",clr='blue', marker='s',    markevery=10, label='CCM')
    p.lwd=1; p.plot ("test_models_6.res",clr='red',  marker='None', markevery=10, label='Unconv01')
    legend()

    subplot(2,3,5)
    ccm = read_table("test_models_4.res")
    res = read_table("test_models_6.res")
    ccm_z0 = array(ccm['z0'])/sqrt(3.0)
    res_z0 = exp(array(res['z0']))/sqrt(3.0)
    ccm_ev = array(ccm['ex'])+array(ccm['ey'])+array(ccm['ez'])
    res_ev = array(res['ex'])+array(res['ey'])+array(res['ez'])
    plot (log(ccm_z0), ccm_ev,'blue')
    plot (log(res_z0), res_ev,'red')

    p.show()

if tst==3:
    p = Plotter()
    p.pq_ty   = 'oct'
    p.evd_ty  = 'oct'
    p.lnplus1 = True
    #p.justone = 4
    #p.fc_c    = 0.1
    p.fc_phi    = M_calc_phi(1,'oct')
    p.fc_poct   = 150.0*sqrt(3.0)
    p.show_k    = True
    p.only_four = True
    #p.set_eps = True
    p.lwd=2; p.plot ("test1.res", clr='blue',  markevery=10, label='U4', draw_fl=False,draw_ros=True)
    #subplot(2,3,3)
    legend()


    lam0  = 0.1
    lam1  = 1.0
    lam2  = 2.0
    x1    = 1.0
    x2    = 2.5
    ev1   = -2.0
    ev2   = 3.0
    psi0  = 2.0
    psi1  = 2.0
    Mso   = 3.0
    Mcs   = 1.0
    g1    = 5.0
    p0    = 1.0
    x0    = log(1.0+p0)
    xmax  = 6.
    edmax = 2.

    idx = 3 if p.only_four else 4
    num = 2 if p.only_four else 3
    subplot (2,num,idx)
    plot    ([0.,0.8],[0.,0.-psi0*0.8],label=r'$\psi_0$')
    plot    ([0.,0.8],[ev1,ev1+psi1*0.8],label=r'$\psi_1$')
    plot    ([0.,edmax],[ev2,ev2],label=r'$\varepsilon_{v2}$')
    legend  (loc='lower right')

    idx = 4 if p.only_four else 5
    subplot (2,num,idx)
    plot    ([x0,xmax],[0.,-lam0*(xmax-x0)],label=r'$\lambda_0$')
    plot    ([x1,xmax],[0.,-lam1*(xmax-x1)],label=r'$\lambda_1$')
    plot    ([x2,xmax],[0.,-lam2*(xmax-x2)],label=r'$\lambda_2$')

    p.show()
