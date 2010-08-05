from msys_plotter    import *
from msys_invariants import *
from msys_readdata   import *
from numpy import array, log, sqrt
from pylab import plot, show, grid, legend, xlabel, ylabel

tst = 10

if tst==0:
    p = Plotter()
    p.plot ("test_models.res", draw_ys=False,draw_fl=False,dpt_out=5,pqty='cam')
    p.plot ('FO1_CTR_01.kgf.pct.dat',fem_res=False)
    p.show()

if tst==1:
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi  = M_calc_phi(1,'cam')
    p.fc_poct = 150.0*sqrt(3.0)
    p.show_k  = True
    #p.plot ("test_models_2.res",clr='red',   marker='+', draw_fl=True,draw_ros=True)
    #p.plot ("test_models_3.res",clr='green', marker='+')
    p.lwd=1; p.plot ("test_models_4.res",clr='blue', marker='s',    markevery=10, label='CCM')
    p.lwd=1; p.plot ("test_models_6.res",clr='red',  marker='None', markevery=10, label='Unconv01')
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
    res = read_table("test_models_6.res")
    v0  = 2.0;
    p0  = 100.0*sqrt(3.0);
    l0  = 0.001;
    l1  = 0.01;
    l3  = 0.01;
    p   = -(array(res['sx'])+array(res['sy'])+array(res['sz']))/sqrt(3.0)
    z0  = array(res['z0'])
    xR3 = array(res['z1'])
    y   = log(array(res['v']))
    pf  = max([max(p),max(exp(z0)),max(exp(xR3))])
    plot   (log(p),y,'k+', label='p', markersize=20)
    plot   (z0,    y,      label='z0')
    plot   (xR3,   y,'ro', label='xR3')
    plot   ([log(p0),log(pf)],[log(v0),log(v0)-l0*log(pf/p0)],'m',marker='.',label='l0')
    plot   ([log(p0),log(pf)],[log(v0),log(v0)-l1*log(pf/p0)],'m',marker='^',label='l1')
    plot   ([log(p0),log(pf)],[log(v0),log(v0)-l3*log(pf/p0)],'m',marker='s',label='l3')
    xlabel (r'$\ln{p}$')
    ylabel (r'$\ln{v}$')
    legend ()
    grid   ()
    show   ()

if tst==4:
    ccm    = read_table("test_models_4.res")
    res    = read_table("test_models_6.res")
    ccm_p  = -(array(ccm['sx'])+array(ccm['sy'])+array(ccm['sz']))/sqrt(3.0)
    p      = -(array(res['sx'])+array(res['sy'])+array(res['sz']))/sqrt(3.0)
    ccm_ev =  (array(ccm['ex'])+array(ccm['ey'])+array(ccm['ez']))
    ev     =  (array(res['ex'])+array(res['ey'])+array(res['ez']))
    plot   (log(p),    ev,    'k+', label='Unconv01', markersize=20)
    plot   (log(ccm_p),ccm_ev,'ro', label='CCM')
    xlabel (r'$\ln{p}$')
    ylabel (r'$\varepsilon_{v}$')
    legend ()
    grid   ()
    show   ()

if tst==5:
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi  = M_calc_phi(1,'cam')
    p.fc_poct = 150.0*sqrt(3.0)
    p.show_k  = True
    #p.plot ("test_models_2.res",clr='red',   marker='+', draw_fl=True,draw_ros=True)
    p.lwd=1; p.plot ("test_models_4_ctt.res", clr='blue', marker='s',    markevery=10, label='CCM: CTT')
    p.lwd=1; p.plot ("test_models_4_pcte.res",clr='red',  marker='None', markevery=10, label='CCM: pcte')
    legend()
    p.show()

if tst==6:
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi  = M_calc_phi(1,'cam')
    p.fc_poct = 150.0*sqrt(3.0)
    p.show_k  = True
    p.lwd=1; p.plot ("test_models_4.res", clr='blue', marker='s',    markevery=10, label='CCM', draw_fl=True,draw_ros=True)
    p.lwd=1; p.plot ("test_models_7.res", clr='red',  marker='+',    markevery=10, label='Unconv02')
    legend()
    p.show()

if tst==7:
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi  = M_calc_phi(1,'cam')
    p.fc_poct = 150.0*sqrt(3.0)
    p.show_k  = True
    p.lwd=1; p.plot ("test_models_4.res", clr='blue', marker='s',    markevery=10, label='CCM', draw_fl=True,draw_ros=True)
    p.lwd=1; p.plot ("test_models_7.res", clr='red',  marker='+',    markevery=10, label='Unconv02')
    #p.lwd=1; p.plot ("test_models_8.res", clr='red',  marker='+',    markevery=10, label='Unconv03')
    legend()
    p.show()

if tst==8:
    l0   = 0.001;
    l1   = 0.005;
    l3   = 0.008;
    #res  = read_table("test_models_8.res")
    res  = read_table("test_models_7.res")
    z0   = array(res['z0'])
    z1   = array(res['z1']) # == xR3
    p    = -(array(res['sx'])+array(res['sy'])+array(res['sz']))/sqrt(3.0)
    v    = array(res['v'])
    x    = log(p)
    y    = log(v)
    x0   = x[0]
    y0   = y[0]
    xf   = max([max(z0),max(z1),max(x)])
    dx   = xf-x0
    xR10 = x0
    xR30 = z1[0]
    plot   ([x0,  xf],[y0,y0-l0*dx],'m',label=r'$\lambda_0$',marker='s')
    plot   ([xR10,xf],[y0,y0-l1*dx],'m',label=r'$\lambda_1$',marker='o')
    plot   ([xR30,xf],[y0,y0-l3*dx],'m',label=r'$\lambda_3$',marker='^')
    plot   (x, y,'b',label=r'p', marker='_',markersize=20)
    plot   (z0,y,'r',label=r'z0',marker='+')
    plot   (z1,y,'g',label=r'z1')
    xlabel (r'$x=\ln{p}$')
    ylabel (r'$y=\ln{v}$')
    legend ()
    grid   ()
    show   ()

if tst==9:
    #l0    = 0.0165
    #l1    = 0.208
    #l3    = 0.6
    #x0    = 4.99
    #v0    = 2.15
    #xR10  = 6.83
    #xR30  = 7.8

    l0    = 0.001
    l1    = 0.005
    l3    = 0.008
    betb  = 100.0
    betbb = 100.0
    v0    = 2.0
    xR10  = 100*sqrt(3.0)
    xR30  = xR10
    x0    = log(100.)

    res   = read_table("test_models_6.res")
    z0    = array(res['z0'])
    z1    = array(res['z1']) # == xR3
    p     = -(array(res['sx'])+array(res['sy'])+array(res['sz']))/sqrt(3.0)
    v     = array(res['v'])
    x     = log(p)
    y     = v
    xf    = 8.5#max([max(log(z0)),max(log(z1)),max(x)])
    dx    = xf-x0
    plot   ([x0,  xf],[v0,v0-l0*dx],'m',label=r'$\lambda_0$',marker='s')
    plot   ([xR10,xf],[v0,v0-l1*dx],'m',label=r'$\lambda_1$',marker='o')
    plot   ([xR30,xf],[v0,v0-l3*dx],'m',label=r'$\lambda_3$',marker='^')
    plot   (x,      y,'b',label=r'p')
    plot   (log(z0),y,'r',label=r'z0',marker='+')
    plot   (log(z1),y,'g',label=r'z1')
    xlabel (r'$x=\ln{p}$')
    ylabel (r'$y=v$')
    legend ()
    grid   ()
    show   ()

if tst==10:
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi  = M_calc_phi(1,'cam')
    p.fc_poct = 150.0*sqrt(3.0)
    p.show_k  = True
    #p.set_eps = True
    #p.lwd=1; p.plot ("test_models_4.res", clr='blue', marker='s',    markevery=10, label='Cam clay', draw_fl=True,draw_ros=True)
    p.lwd=1; p.plot ("test_models_7.res", clr='red',  marker='+',    markevery=10, label='Unconventional 2')
    #p.lwd=1; p.plot ("test_models_8.res", clr='red',  marker='+',    markevery=10, label='Unconventional 1')
    subplot(2,3,3)
    legend()
    p.show()
    #p.save_eps("dilatancy")
