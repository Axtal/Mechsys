# Copyright (C) 2009 Dorival M Pedroso
# ------------------------------------
# Models

from os.path    import basename
from numpy      import *
from pylab      import *
from invariants import *

class Plotter:
    # Plot results
    # obj: 'datafilename'
    #      TriaxTest
    #      Domain
    # ============
    def plot(self,
        obj,                     # object: data_filename, TriaxTest, Domain
        fem_res         = True,  # data_file is .res file from FEM analysis
        mdl             = None,  # constitutive model
        eid             = 0,     # element ID for output (if Domain)
        div_by_p        = True,  # divide q by p ?
        pqty            = 'cam', # use Cam clay invariants ?
        q_neg_extension = True,  # multiply q by -1 for extension (t<0, where t=sin(3th)
        closed_form     = False, # show closed-form solution (only with mdl)
        show_k          = False, # show k=dq/dp ?
        draw_ys         = True,  # draw yield surface
        draw_fl         = True,  # draw failure line
        dpt_out         = 1,     # draw YS after every delta point
        use_gen_p_q     = False, # use mdl.gen_p_q instead of contour when drawing YS
        npts_ys         = 20,    # number of points for drawing YS
        npts_fl         = 10,    # number of points for drawing failure line
        with_strains    = True,  # with strains
        isxyz           = (-1,0), # indices for sxyz plot, use negative numbers for principal components
        lst_plot_coor   = 'y'):  # last plot coordinate: x, y, z, xy

        if obj.__class__.__name__=='list': # list with tree lists: Sig, Eps, and Ivs
            Sig = obj[0]
            Ivs = obj[1]
            dot = ['r-','r-','r-','r-','r-','r-','r-','r-','r-']
            with_strains = False

        elif obj.__class__.__name__=='str': # data file/fem results
            # load data
            f   = open(obj,'r')
            h   = f.readline().split() # header
            Sig = []
            Eps = []
            sq2 = sqrt(2.0)
            if fem_res:
                for l in f:
                    r = l.split()
                    Sig.append(matrix([[float(r[h.index('sx' )])],
                                       [float(r[h.index('sy' )])],
                                       [float(r[h.index('sz' )])],
                                       [float(r[h.index('sxy')])*sq2]]))
                    Eps.append(matrix([[float(r[h.index('ex' )])],
                                       [float(r[h.index('ey' )])],
                                       [float(r[h.index('ez' )])],
                                       [float(r[h.index('exy')])*sq2]]))
                dot = ['r-','r-','r-','r-','r-','r-','r-','r-','r-']
            else:
                res = basename(obj).split('.')
                kgf = res[1]=='kgf'         # stresses in kgf/cm^2 ?
                pct = res[2]=='pct'         # strains in percentage ?
                mul = 98.0  if kgf else 1.0 # multiplier for stresses
                div = 100.0 if pct else 1.0 # divider for strains
                for l in f:
                    r = l.split()
                    Sig.append(mul*matrix([[-float(r[h.index('Sx' )])],
                                           [-float(r[h.index('Sy' )])],
                                           [-float(r[h.index('Sz' )])],
                                           [-float(r[h.index('Sxy')])*sq2]]))
                    Eps.append(    matrix([[-float(r[h.index('Ex' )])/div],
                                           [-float(r[h.index('Ey' )])/div],
                                           [-float(r[h.index('Ez' )])/div],
                                           [-float(r[h.index('Exy')])*sq2/div]]))
                dot = ['k+','k+','k+','k+','k+','k+','k+','k+','k+']
            draw_ys = False
            draw_fl = False

        elif obj.__class__.__name__=='TriaxTest':
            # point simulation results
            Sig = obj.Sig
            Eps = obj.Eps
            Ivs = obj.Ivs
            dot = ['r-','r-','r-','r-','r-','r-','r-','r-','r-']

        elif obj.__class__.__name__=='Domain':
            # FEM simulation results
            Sig = []
            Eps = []
            Ivs = []
            sq2 = sqrt(2.0)
            niv = 0
            for key in obj.Eout[eid].keys():
                for i in range(11):
                    if key=='z%d'%i: niv += 1
            z = zeros(niv)
            for k in range(len(obj.Eout[eid]['sx'])):
                Sig.append(matrix([[obj.Eout[eid]['sx' ][k]],
                                   [obj.Eout[eid]['sy' ][k]],
                                   [obj.Eout[eid]['sz' ][k]],
                                   [obj.Eout[eid]['sxy'][k]*sq2]]))
                Eps.append(matrix([[obj.Eout[eid]['ex' ][k]],
                                   [obj.Eout[eid]['ey' ][k]],
                                   [obj.Eout[eid]['ez' ][k]],
                                   [obj.Eout[eid]['exy'][k]*sq2]]))
                for i in range(niv): z[i] = obj.Eout[eid]['z%d'%i][k]
                Ivs.append(z.copy())
            dot = ['b-','b-','b-','b-','b-','b-','b-','b-','b-']

        else: raise Exception('Plotter:plot: object type==%s is invalid'%obj.__class__.__name__)

        # calculate additional variables
        np         = len(Sig) # number of points
        P,  Q,  T  = zeros(np), zeros(np), zeros(np)
        Sx, Sy, Sz = zeros(np), zeros(np), zeros(np)
        Ex, Ey, Ez = zeros(np), zeros(np), zeros(np)
        S1, S2, S3 = zeros(np), zeros(np), zeros(np) # principal stresses
        Sa, Sb, Sc = zeros(np), zeros(np), zeros(np) # octahedral coordinates
        Si, Sj     = zeros(np), zeros(np)
        pmax, qmax = sig_calc_p_q(Sig[0], Type=pqty) # max p and q
        imax       = 0
        if with_strains:
            Ev,  Ed    = zeros(np), zeros(np)
            E1, E2, E3 = zeros(np), zeros(np), zeros(np) # principal strains
            if closed_form:
                Evs, Eds = zeros(np), zeros(np) # closed-form solution
        for i in range(np):
            P [i], Q [i]        = sig_calc_p_q   (Sig[i],Type=pqty)
            T [i]               = sig_calc_t     (Sig[i])
            s123                = sig_calc_s123  (Sig[i], do_sort=True)
            Sa[i], Sb[i], Sc[i] = s123_calc_oct  (s123)
            S1[i], S2[i], S3[i] = s123  [0], s123  [1], s123  [2]
            Sx[i], Sy[i], Sz[i] = Sig[i][0], Sig[i][1], Sig[i][2]
            if with_strains:
                Ev[i], Ed[i]        = eps_calc_ev_ed (Eps[i])
                e123                = eps_calc_e123  (Eps[i])
                E1[i], E2[i], E3[i] = e123  [0], e123  [1], e123  [2]
                Ex[i], Ey[i], Ez[i] = Eps[i][0], Eps[i][1], Eps[i][2]
                if closed_form:
                    sig0 = Sig[0]
                    if i>0: Evs[i], Eds[i] = mdl.closed_form_ev_ed(sig0,Sig[i])
            if q_neg_extension and T[i]<0.0: Q[i] = -Q[i]
            if P[i]>=pmax: pmax, imax = P[i], i
            if Q[i]>=qmax: qmax, imax = Q[i], i
            if isxyz[0]<0 or isxyz[1]<0:
                Si[i] = s123[-isxyz[0]]
                Sj[i] = s123[-isxyz[1]]
                ikeys = ['1','2','3']
            else:
                Si[i] = Sig[i][isxyz[0]]
                Sj[i] = Sig[i][isxyz[1]]
                ikeys = ['x','y','z']

        # convert strains to percentage
        if with_strains:
            Ev *= 100.0
            Ed *= 100.0
            if closed_form:
                Evs *= 100.0
                Eds *= 100.0

        # constants
        rc('text', usetex=True)    # set LaTeX
        rc('font', family='serif') # set font
        fsz   = 14                 # font size
        lwd   = 2                  # linewidth
        cform = 'g*'               # closed-form format
        if with_strains:
            nhplt = 3 # number of horizontal plots
            nvplt = 3 # number of vertical plots
        else:
            nhplt = 1 # number of horizontal plots
            nvplt = 3 # number of vertical plots

        # q/p, Ed
        iplot = 1
        if with_strains:
            subplot(nhplt,nvplt,iplot);  iplot += 1
            Y    = Q/P if div_by_p else Q
            Ylbl = r'$q_{%s}/p_{%s}$'%(pqty,pqty) if div_by_p else r'$q_{%s}$'%(pqty)
            plot   (Ed, Y, dot[0],linewidth=lwd)
            xlabel (r'$\varepsilon_d$ [\%]',fontsize=fsz);  ylabel(Ylbl,fontsize=fsz);  grid()
            if closed_form: plot (Eds, Y, cform,linewidth=1)

        # q/p, Ev
        if with_strains:
            subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ev, Y, dot[1],linewidth=lwd)
            xlabel (r'$\varepsilon_v$ [\%]',fontsize=fsz);  ylabel(Ylbl,fontsize=fsz);  grid()
            if closed_form: plot (Evs, Y, cform,linewidth=1)

        # p, q
        subplot(nhplt,nvplt,iplot);  iplot += 1
        axhline(0.0,color='black'); axvline(0.0,color='black')
        plot   (P, Q, dot[2],linewidth=lwd)
        xlabel (r'$p_{%s}$'%(pqty),fontsize=fsz);  ylabel(r'$q_{%s}$'%(pqty),fontsize=fsz);  grid()
        axis   ('equal')
        if draw_ys:
            ptout = 0
            for i in range(np):
                if i>=ptout:
                    if use_gen_p_q:
                        p, q = mdl.gen_p_q (Ivs[i], T[i], pqty)
                        plot (p, q, 'g-')
                    else:
                        self.pq_ysurf_fline (mdl, pmax, qmax, T[i], Ivs[i], pqty, np=npts_ys)
                    plot ([P[i]], [Q[i]], 'yo')
                    ptout += dpt_out
        if show_k:
            k = (Q[-1]-Q[0])/(P[-1]-P[0])
            text ((P[0]+P[-1])/2.0,(Q[0]+Q[-1])/2.0,'k = %g'%k,fontsize=14,color='black',ha='left')
        if draw_fl:
            self.pq_ysurf_fline (mdl, pmax, qmax, 1.0, Ivs[0], pqty, fcrit=True, np=npts_fl)

        # Ed, Ev
        if with_strains:
            subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ed, Ev, dot[3],linewidth=lwd)
            xlabel (r'$\varepsilon_d$ [\%]',fontsize=fsz);  ylabel(r'$\varepsilon_v$ [\%]',fontsize=fsz); grid()
            if closed_form: plot (Eds, Evs, cform,linewidth=1)

        # lnp, Ev
        if with_strains:
            subplot(nhplt,nvplt,iplot);  iplot += 1
            X = log(P)
            plot   (X, Ev, dot[4],linewidth=lwd)
            xlabel (r'$\ln{(p_{%s})}$'%(pqty),fontsize=fsz);  ylabel(r'$\varepsilon_v$ [\%]',fontsize=fsz);  grid()
            if closed_form: plot (X, Evs, cform,linewidth=1)

        # Sa, Sb
        subplot(nhplt,nvplt,iplot);  iplot += 1
        plot   (Sa, Sb, dot[5],linewidth=lwd)
        axis   ('equal')
        xlabel (r'$\sigma_a$',fontsize=fsz);  ylabel(r'$\sigma_b$',fontsize=fsz);  grid()
        radius = max(sqrt(Sa**2.0+Sb**2.0))
        self.oct_rosette(radius)
        if draw_ys:
            ptout = 0
            for i in range(np):
                if i>=ptout:
                    sa, sb, sc = s123_calc_oct(array([S1[i],S2[i],S3[i]]))
                    self.oct_ysurf_fline (mdl, radius, P[i], Ivs[i], pqty, fcrit=False)
                    plot ([sa], [sb], 'yo')
                    ptout += dpt_out

        # Ek, Q/P
        if with_strains:
            subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (E1, Y, dot[6],linewidth=lwd,linestyle='-')
            plot   (E2, Y, dot[6],linewidth=lwd,linestyle='--')
            plot   (E3, Y, dot[6],linewidth=lwd,linestyle='-.')
            xlabel (r'$\varepsilon_1$[--], $\varepsilon_2$[- -], $\varepsilon_3$[- .]',fontsize=fsz)
            ylabel (Ylbl,fontsize=fsz);  grid()

        # Ek, Sk
        if with_strains:
            subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ex, -Sx, dot[7],linewidth=lwd,linestyle='-')
            plot   (Ey, -Sy, dot[7],linewidth=lwd,linestyle='--')
            plot   (Ez, -Sz, dot[7],linewidth=lwd,linestyle='-.')
            xlabel (r'$\varepsilon_x$[--], $\varepsilon_y$[- -], $\varepsilon_z$[- .]',fontsize=fsz)
            ylabel (r'$-\sigma_x$[--], $-\sigma_y$[- -], $-\sigma_z$[- .]',fontsize=fsz);  grid()

        # sqrt(2.0)*Si, Sj
        subplot(nhplt,nvplt,iplot);  iplot += 1
        if draw_ys and not use_gen_p_q:
            ptout = 0
            for i in range(np):
                if i>=ptout:
                    self.sxyz_ysurf_fline (mdl, Sig[imax], T[i], Ivs[i], np=npts_ys)
                    plot ([-sqrt(2.0)*Si[i]], [-Sj[i]], 'yo')
                    ptout += dpt_out
        if draw_fl:
            self.sxyz_ysurf_fline (mdl, Sig[imax], 1.0, Ivs[0], fcrit=True, np=npts_fl)
        #axhline(0.0,color='black'); axvline(0.0,color='black')
        plot   (-sqrt(2.0)*Si, -Sj, dot[8], linewidth=lwd)
        xlabel (r'$-\sqrt{2}\sigma_%s$'%(ikeys[abs(isxyz[0])]),fontsize=fsz);  ylabel(r'$-\sigma_%s$'%(ikeys[abs(isxyz[1])]),fontsize=fsz);  grid()
        axis   ('scaled')

        #new = Ivs[1:]
        #old = Ivs[:-1]
        #for i in range(len(new)):
            #print new[i]-old[i]

    # Plot node
    # =========
    def plot_node(self, obj):
        if obj.__class__.__name__=='str': # data file/fem results
            # load data
            f = open(obj,'r')
            h = f.readline().split() # header
            Ux, Uy = [], []
            Fx, Fy = [], []
            for l in f:
                r = l.split()
                Ux.append(-float(r[h.index('ux')]))
                Uy.append(-float(r[h.index('uy')]))
                Fx.append(-float(r[h.index('fx')]))
                Fy.append(-float(r[h.index('fy')]))
        else: raise Exception('Plotter:plot_node works only with obj=str')

        # constants
        rc('text', usetex=True)    # set LaTeX
        rc('font', family='serif') # set font
        fsz   = 14                 # font size
        lwd   = 2                  # linewidth
        nhplt = 1                  # number of horizontal plots
        nvplt = 2                  # number of vertical plots
        iplot = 1

        # ux, fx
        subplot(nhplt,nvplt,iplot);  iplot += 1
        plot   (Ux,Fx,'r-',linewidth=lwd)
        xlabel (r'$u_x$',fontsize=fsz)
        ylabel (r'$f_x$',fontsize=fsz);  grid()

        # uy, fy
        subplot(nhplt,nvplt,iplot);  iplot += 1
        plot   (Uy,Fy,'r-',linewidth=lwd)
        xlabel (r'$u_y$',fontsize=fsz)
        ylabel (r'$f_y$',fontsize=fsz);  grid()

    # Show figure
    # ===========
    def show(self): show()

    # Plot octahedral rosette
    # =======================
    def oct_rosette(self, radius, coef=0.2):
        l1 = (                    0.0  ,      radius            ) # line: 1 end points
        l2 = (-coef*radius*cos(pi/6.0) ,-coef*radius*sin(pi/6.0)) # line: 2 end points
        l3 = ( coef*radius*cos(pi/6.0) ,-coef*radius*sin(pi/6.0)) # line: 3 end points
        l4 = (     -radius*cos(pi/6.0) ,      radius*sin(pi/6.0)) # line: 4 = neg 1 end points
        lo = (     -radius*cos(pi/3.0) ,      radius*sin(pi/3.0)) # line: origin of cylindrical system
        # main lines
        plot([0.0,l1[0]],[0.0,l1[1]],'k-')
        plot([0.0,l2[0]],[0.0,l2[1]],'k-')
        plot([0.0,l3[0]],[0.0,l3[1]],'k-')
        # reference
        plot([0.0,l4[0]],[0.0, l4[1]],'k--')
        plot([0.0,lo[0]],[0.0, lo[1]],'k--')
        # text
        text(l1[0],l1[1],r'$-\sigma_1,\theta=+30^o$')
        text(l2[0],l2[1],r'$-\sigma_2$',ha='right')
        text(l3[0],l3[1],r'$-\sigma_3$')
        text(lo[0],lo[1],r'$\theta=0^o$',ha='center')
        text(l4[0],l4[1],r'$\theta=-30^o$',ha='right')

    # Plot yield surface or failure line in p-q plane
    # ===============================================
    def pq_ysurf_fline(self, mdl, pmax, qmax, t, ivs, pqty, fcrit=False, np=20):
        try:    cbar = mdl.cbar
        except: cbar = 0.0
        sig0   = (cbar/sqrt(3.0))*matrix([[1.0],[1.0],[1.0],[0.0]])
        p0, q0 = sig_calc_p_q (sig0, pqty)
        dp     = 1.6*(pmax-p0)/np
        dq     = 1.3*(qmax-q0)/np
        p      = zeros ((np,np))
        q      = zeros ((np,np))
        f      = zeros ((np,np))
        for i in range(np):
            for j in range(np):
                p[i,j] = p0 + i*dp
                q[i,j] = q0 + j*dq
                s123   = pqt_calc_s123 (p[i,j], q[i,j], t, pqty)
                sig    = matrix([[s123[0]],[s123[1]],[s123[2]],[0.0]])
                if fcrit: clr, f[i,j] = 'm', mdl.failure_crit (sig)
                else:     clr, f[i,j] = 'g', mdl.yield_func   (sig, ivs)
                if t<0.0: q[i,j] = -q[i,j]
        contour (p, q, f, [0.0], linewidths=1, colors=clr)
        if fcrit: text(p0, 0.0, r'$p_0=%g$'%p0)

    # Plot yield surface or failure line in sxyz plane
    # ================================================
    def sxyz_ysurf_fline(self, mdl, sigf, t, ivs, fcrit=False, np=20):
        try:    cbar = mdl.cbar
        except: cbar = 0.0
        sig0   = (cbar/sqrt(3.0))*matrix([[1.0],[1.0],[1.0],[0.0]])
        p0, q0 = sig_calc_p_q (sig0)
        pf, qf = sig_calc_p_q (sigf)
        xmin   = sqrt(2.0)*(p0/sqrt(3.0)+q0/sqrt(6.0))
        ymin   = 0.0
        xmax   = sqrt(2.0)*(pf/sqrt(3.0)+qf/sqrt(6.0))
        ymax   = (pf/sqrt(3.0)+2.0*qf/sqrt(6.0))
        dx     = 1.5*(xmax-xmin)/np
        dy     = 1.3*(ymax-ymin)/np
        x      = zeros ((np,np))
        y      = zeros ((np,np))
        f      = zeros ((np,np))
        for i in range(np):
            for j in range(np):
                x[i,j] = xmin + i*dx
                y[i,j] = ymin + j*dy
                sig    = matrix([[-y[i,j]], [-x[i,j]/sqrt(2.0)], [-x[i,j]/sqrt(2.0)], [0.0]])
                if fcrit: clr, f[i,j] = 'm', mdl.failure_crit (sig)
                else:     clr, f[i,j] = 'g', mdl.yield_func   (sig, ivs)
        plot    ([0.0,sqrt(2.0)*xmax],[0.0,xmax],'k-') # hydrostatic line
        contour (x, y, f, [0.0], linewidths=1, colors=clr)
        if fcrit: text(-sqrt(2.0)*sig0[0,0], 0.0, r'$\sigma_0=%g$'%(sig0[0,0]))

    # Plot yield surface or failure line in octahedral plane
    # ======================================================
    def oct_ysurf_fline(self, mdl, radius, p, ivs, pqty, fcrit=False, np=20, coef=0.2):
        samin =  -1.2*radius*cos(pi/6.0)
        samax =  coef*radius*cos(pi/6.0)
        sbmin = -coef*radius*sin(pi/6.0)
        sbmax =   1.2*radius
        sc,q  = convert_p_q(p,0.0,1.0,pqty,'oct')
        dsa   = 1.2*(samax-samin)/np
        dsb   = 1.2*(sbmax-sbmin)/np
        sa    = zeros ((np,np))
        sb    = zeros ((np,np))
        f     = zeros ((np,np))
        for i in range(np):
            for j in range(np):
                sa[i,j] = samin + i*dsa
                sb[i,j] = sbmin + j*dsb
                s123    = oct_calc_s123 (sa[i,j], sb[i,j], sc)
                sig     = matrix([[s123[0]],[s123[1]],[s123[2]],[0.0]])
                if fcrit: clr, f[i,j] = 'm', mdl.failure_crit (sig)
                else:     clr, f[i,j] = 'g', mdl.yield_func   (sig, ivs)
        contour (sa, sb, f, [0.0], linewidths=1, colors=clr)


# Test Plotter
# ============
if __name__=='__main__':
    from camclay       import *
    from elastoplastic import *

    # variables
    phi  = M_calc_phi(1.0,'cam')
    tst  = 'DP'
    pat  = 'psa'
    b    = 1.0#0.05
    sphi = sin(phi*pi/180.0)

    # path
    if pat=='psa':
        s1  = -200.0
        s3  = s1*(1.0-sphi)/(1.0+sphi)
        s2  = b*s1 + (1.0-b)*s3
        Sig = [matrix([[-100.0],[-100.0],[-100.0],[0.0]]),
               matrix([[s1],[s2],[s3],[0.0]])]
    elif pat=='com':
        Sig = [matrix([[-100.0],[-100.0],[-100.0],[0.0]]),
               matrix([[-250.0],[-100.0],[-100.0],[0.0]])]
    elif pat=='ext':
        Sig = [matrix([[-100.0],[-100.0],[-100.0],[0.0]]),
               matrix([[-100.0],[-200.0],[-200.0],[0.0]])]

    # model
    if tst=='CAM':
        mdl = CamClay({'lam':0.01, 'kap':0.001, 'nu':0.3, 'phi':phi})
        z0f = mdl.init_ivs({'v':2.0},Sig[-1])
        z0i = mdl.init_ivs({'v':2.0},Sig[0])
    elif tst=='MC' or tst=='DP':
        mdl = ElastoPlastic({'E':6000.0, 'nu':0.3, 'fcty':tst, 'c':0.0, 'phi':phi})
        z0f = mdl.init_ivs({},Sig[-1])
        z0i = mdl.init_ivs({},Sig[0])

    # plot
    Ivs = [z0i,z0f]
    plt = Plotter()
    plt.plot([Sig,Ivs],mdl,draw_ys=True,draw_fl=True,npts_ys=20,npts_fl=20,pqty='oct',isxyz=(-2,0))

    # add gradients
    s123                = sig_calc_s123 (Sig[-1], do_sort=True)
    sa, sb, sc          = s123_calc_oct (s123)
    p, q                = sig_calc_p_q  (Sig[-1])
    V, Y                = mdl.gradients (Sig[-1],Ivs[-1])
    dfdp, dfdq          = sig_calc_p_q  (V)
    dfds123             = sig_calc_s123 (V)
    dfdsa, dfdsb, dfdsc = s123_calc_oct (dfds123)

    print 'sigf =', Sig[-1].T, s123.T

    subplot (1,3,1)
    plt.pq_ysurf_fline(mdl, p, q,  1.0, Ivs[-1], 'oct', fcrit=True, np=20)
    plt.pq_ysurf_fline(mdl, p, q, -1.0, Ivs[-1], 'oct', fcrit=True, np=20)
    if b<=0.5: arrow (p,  q, 200.0*dfdp,  200.0*dfdq)
    else:      arrow (p, -q, 200.0*dfdp, -200.0*dfdq)

    subplot (1,3,2)
    arrow   (sa, sb, 200.0*dfdsa, 200.0*dfdsb)

    show()
