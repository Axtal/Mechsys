########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Dorival M. Pedroso                                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

from os.path import basename
from numpy   import sqrt, matrix, zeros, log, cos, pi, sin, tan
from pylab   import rc, subplot, plot, xlabel, ylabel, grid, axhline, axvline, axis, text, contour, show, rcParams, savefig, axes, legend
from msys_invariants import *
from msys_readdata   import *

class Plotter:
    # Constructor
    # ===========
    def __init__(self):
        # data
        self.show_k    = False                                    # show k=dq/dp ?
        self.div_by_p  = False                                    # divide q by p ?
        self.log_p     = True                                     # use log(p) instead of p ?
        self.q_neg_ext = False                                    # multiply q by -1 for extension (t<0, where t=sin(3th)
        self.pq_ty     = 'cam'                                    # invariants type
        self.fc_ty     = ['MC','MN','LD','VM']                    # failure criteria type (VM, DP, MC, MN)
        self.fc_clr    = ['k', 'k', 'k', 'k']                     # failure criteria colors
        self.fc_ls     = ['solid', 'dashed', 'dotted', 'dashdot'] # failure criteria linestyles
        self.fc_p      = 100.0                                    # p(oct) to be used when plotting FC in octahedral plane
        self.fc_t      = 1.0                                      # t=sin(3th) to be used when plotting FC in octahedral plane
        self.fc_cu     = 10.0                                     # cohesion for VM
        self.fc_phi    = 25.0                                     # friction angle for FC
        self.fc_c      = 0.0                                      # cohesion for FC
        self.fc_np     = 20                                       # number of points for drawing failure line
        self.isxyz     = (-1,0)                                   # indices for sxyz plot, use negative numbers for principal components
        self.devplot   = True                                     # plot s3-s1, s3-s2 instead of Ek, Sk
        self.pcte      = False                                    # pcte in Ev x p (logp) plot?
        self.justone   = -1                                       # all plots = -1
        self.maxed     = -1                                       # max Ed to stop the lines
        self.maxev     = -1                                       # max Ev to stop the lines
        self.maxidx    = -1                                       # max index of data to plot
        self.axesdat   = [0.11,0.15,0.98-0.11,0.95-0.15]          # geometry of figure
        self.ms        = 4                                        # marker size in points
        self.set_eps   = False                                    # set geometry for eps (savefig)
        self.only_six  = True                                     # only 2 x 3 instead of 3 x 3 plots?
        self.proport   = 0.75                                     # proportion (Aesthetic ratio): golden_mean = (sqrt(5)-1.0)/2.0 0.628
        self.lwd       = 2                                        # linewidth

    # Plot results
    # ============
    def plot(self, filename, clr='red', draw_fl=False, draw_ros=True, txtlst=True, txtmax=True,
                             label=None, marker='None', markevery=None):
        # load data
        dat = read_table(filename)
        Sig = []
        Eps = []
        sq2 = sqrt(2.0)
        if dat.has_key('Sx'): # old data file
            res = basename(filename).split('.')
            kgf = res[1]=='kgf'         # stresses in kgf/cm^2 ?
            pct = res[2]=='pct'         # strains in percentage ?
            mul = 98.0  if kgf else 1.0 # multiplier for stresses
            div = 100.0 if pct else 1.0 # divider for strains
            for i in range(len(dat['Sx'])):
                Sig.append(mul*matrix([[-float( dat['Sx' ][i] )],
                                       [-float( dat['Sy' ][i] )],
                                       [-float( dat['Sz' ][i] )],
                                       [-float( dat['Sxy'][i] )*sq2]]))
                Eps.append(    matrix([[-float( dat['Ex' ][i] )/div],
                                       [-float( dat['Ey' ][i] )/div],
                                       [-float( dat['Ez' ][i] )/div],
                                       [-float( dat['Exy'][i] )*sq2/div]]))
        else:
            for i in range(len(dat['sx'])):
                Sig.append(matrix([[float( dat['sx' ][i] )],
                                   [float( dat['sy' ][i] )],
                                   [float( dat['sz' ][i] )],
                                   [float( dat['sxy'][i] if dat.has_key('sxy') else 0.0 )*sq2]]))
                Eps.append(matrix([[float( dat['ex' ][i] )],
                                   [float( dat['ey' ][i] )],
                                   [float( dat['ez' ][i] )],
                                   [float( dat['exy'][i] if dat.has_key('sxy') else 0.0 )*sq2]]))

        # calculate additional variables
        np         = len(Sig) # number of points
        P,  Q,  T  = zeros(np), zeros(np), zeros(np)
        Sx, Sy, Sz = zeros(np), zeros(np), zeros(np)
        Ex, Ey, Ez = zeros(np), zeros(np), zeros(np)
        S1, S2, S3 = zeros(np), zeros(np), zeros(np) # principal stresses
        Sa, Sb, Sc = zeros(np), zeros(np), zeros(np) # octahedral coordinates
        Si, Sj     = zeros(np), zeros(np)
        Ev, Ed     = zeros(np), zeros(np)
        E1, E2, E3 = zeros(np), zeros(np), zeros(np) # principal strains
        for i in range(np):
            P [i], Q [i]        = sig_calc_p_q   (Sig[i], Type=self.pq_ty)
            T [i]               = sig_calc_t     (Sig[i])
            s123                = sig_calc_s123  (Sig[i], do_sort=True)
            Sa[i], Sb[i], Sc[i] = s123_calc_oct  (s123)
            S1[i], S2[i], S3[i] = s123  [0], s123  [1], s123  [2]
            Sx[i], Sy[i], Sz[i] = Sig[i][0], Sig[i][1], Sig[i][2]
            Ev[i], Ed[i]        = eps_calc_ev_ed (Eps[i])
            e123                = eps_calc_e123  (Eps[i])
            E1[i], E2[i], E3[i] = e123  [0], e123  [1], e123  [2]
            Ex[i], Ey[i], Ez[i] = Eps[i][0], Eps[i][1], Eps[i][2]
            if self.q_neg_ext and T[i]<0.0: Q[i] = -Q[i]
            if self.isxyz[0]<0 or self.isxyz[1]<0:
                Si[i] = s123[-self.isxyz[0]]
                Sj[i] = s123[-self.isxyz[1]]
                ikeys = ['1','2','3']
            else:
                Si[i] = Sig[i][self.isxyz[0]]
                Sj[i] = Sig[i][self.isxyz[1]]
                ikeys = ['x','y','z']
            Ev[i] *= 100.0 # convert strains to percentage
            Ed[i] *= 100.0

        # constants
        rc('text', usetex=True)    # set LaTeX
        rc('font', family='serif') # set font
        nhplt = 3                  # number of horizontal plots
        nvplt = 3                  # number of vertical plots
        iplot = 1                  # index of plot
        imaQ  = Q.argmax()         # index of q_max
        ilst  = len(P)-1           # index of last point 
        ima   = ilst if self.maxidx<0 else self.maxidx # max index of data to plot
        if self.maxed>0:
            imaEd = ima
            for k, ed in enumerate(Ed):
                if ed>=self.maxed:
                    imaEd = k
                    break
            ima = imaEd if imaEd<ima else ima
        if self.maxev>0:
            imaEv = ima
            for k, ev in enumerate(Ev):
                if ev>=self.maxev:
                    imaEv = k
                    break
            ima = imaEv if imaEv<ima else ima
        if self.only_six:
            nhplt = 2
            nvplt = 3

        # set for eps
        if self.set_eps:
            if self.justone>=0:
                self.set_fig_for_eps(multiplot=False)
                if self.justone==5: axes([0.05,0.05,.95,.95])
                else: axes(self.axesdat) # this needs to be after set_eps
            else:
                self.set_fig_for_eps(multiplot=True)
                axes([0.,0.,0.99,0.99]) # this needs to be after set_eps

        # 0) q/p, Ed ---------------------------------------------------------------------------
        if self.justone==0 or self.justone<0:
            Y    = Q/P if self.div_by_p else Q
            Ylbl = r'$q_{%s}/p_{%s}$'%(self.pq_ty,self.pq_ty) if self.div_by_p else r'$q_{%s}$'%(self.pq_ty)
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot (Ed[:ima], Y[:ima], color=clr, lw=self.lwd, label=label, marker=marker, markevery=markevery, ms=self.ms)
            if imaQ<=ima: plot (Ed[imaQ], Y[imaQ], '^', color=clr)
            if ilst<=ima: plot (Ed[-1],   Y[-1],   '^', color=clr)
            xlabel (r'$\varepsilon_d$ [\%]');  ylabel(Ylbl);  grid()
            if txtmax and imaQ<=ima: text (Ed[imaQ], Y[imaQ], '%.2f'%Y[imaQ], fontsize=8)
            if txtlst and ilst<=ima: text (Ed[-1],   Y[-1],   '%.2f'%Y[-1],   fontsize=8)

        # 1) q/p, Ev ---------------------------------------------------------------------------
        if self.justone==1 or self.justone<0:
            Y    = Q/P if self.div_by_p else Q
            Ylbl = r'$q_{%s}/p_{%s}$'%(self.pq_ty,self.pq_ty) if self.div_by_p else r'$q_{%s}$'%(self.pq_ty)
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ev[:ima], Y[:ima], lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            xlabel (r'$\varepsilon_v$ [\%]');  ylabel(Ylbl);  grid()

        # 2) p, q ---------------------------------------------------------------------------
        if self.justone==2 or self.justone<0:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            axhline (0.0,color='black'); axvline(0.0,color='black')
            if draw_fl: self.pq_fline(0.1*min(P),2.0*max(P),min(Q),max(Q))
            plot    (P[:ima], Q[:ima], lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            xlabel  (r'$p_{%s}$'%(self.pq_ty));  ylabel(r'$q_{%s}$'%(self.pq_ty));  grid()
            axis    ('equal')
            if self.show_k:
                k = (Q[-1]-Q[0])/(P[-1]-P[0])
                text ((P[0]+P[-1])/2.0,(Q[0]+Q[-1])/2.0,'k = %g'%k,color='black',ha='left', fontsize=10)

        # 3) Ed, Ev ---------------------------------------------------------------------------
        if self.justone==3 or self.justone<0:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ed[:ima], Ev[:ima], lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            xlabel (r'$\varepsilon_d$ [\%]');  ylabel(r'$\varepsilon_v$ [\%]'); grid()

        # 4) lnp, Ev ---------------------------------------------------------------------------
        if self.justone==4 or self.justone<0:
            if self.log_p:
                X    = log(P[:ima])
                xlbl = r'$\ln{(p_{%s})}$'%(self.pq_ty)
            else:
                X    = P[:ima]
                xlbl = r'$p_{%s}$'%(self.pq_ty)
            if self.pcte:
                for k, x in enumerate(X): X[k] = X[0]
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (X[:ima], Ev[:ima], lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            xlabel (xlbl);  ylabel(r'$\varepsilon_v$ [\%]');  grid()

        # 5) Sa, Sb ---------------------------------------------------------------------------
        if self.justone==5 or self.justone<0:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            if draw_ros: self.oct_rosette(min(Sa),max(Sa),min(Sb),max(Sb))
            if draw_fl:  self.oct_fline  (min(Sa),max(Sa),min(Sb),max(Sb))
            plot   (Sa[:ima], Sb[:ima], color=clr, lw=self.lwd, label=label, marker=marker, markevery=markevery, ms=self.ms)
            if imaQ<=ima: plot (Sa[imaQ], Sb[imaQ], '^', color=clr)
            if ilst<=ima: plot (Sa[ilst], Sb[ilst], '^', color=clr)
            axis ('equal')
            axis ('off')

        # 6) Ek, Q/P ---------------------------------------------------------------------------
        if self.justone==6 or self.justone<0 and not self.only_six:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (E1[:ima], Y[:ima], lw=self.lwd,linestyle='-',  color=clr)
            plot   (E2[:ima], Y[:ima], lw=self.lwd,linestyle='--', color=clr)
            plot   (E3[:ima], Y[:ima], lw=self.lwd,linestyle='-.', color=clr)
            xlabel (r'$\varepsilon_1$[--], $\varepsilon_2$[- -], $\varepsilon_3$[- .]')
            ylabel (Ylbl);  grid()

        if self.devplot:
            if self.justone==7 or self.justone<0 and not self.only_six:
                # 7) s3-s1, s3-s2
                if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
                if draw_fl: self.s123_fline(S3[imaQ])
                plot   (S3[:ima]-S1[:ima], S3[:ima]-S2[:ima], lw=self.lwd,linestyle='-', color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
                xlabel (r'$\sigma_3-\sigma_1$');  ylabel(r'$\sigma_3-\sigma_2$');  grid()
                axis   ('equal')
        else:
            # 7) Ek, Sk ---------------------------------------------------------------------------
            if self.justone==7 or self.justone<0 and not self.only_six:
                if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
                plot   (Ex[:ima], -Sx[:ima], lw=self.lwd,linestyle='-', color=clr)
                plot   (Ey[:ima], -Sy[:ima], lw=self.lwd,linestyle='--', color=clr)
                plot   (Ez[:ima], -Sz[:ima], lw=self.lwd,linestyle='-.', color=clr)
                xlabel (r'$\varepsilon_x$[--], $\varepsilon_y$[- -], $\varepsilon_z$[- .]')
                ylabel (r'$-\sigma_x$[--], $-\sigma_y$[- -], $-\sigma_z$[- .]');  grid()

        # 8) sqrt(2.0)*Si, Sj ---------------------------------------------------------------------------
        if self.justone==8 or self.justone<0 and not self.only_six:
            if self.justone<0: self.ax = subplot (nhplt,nvplt,iplot);  iplot += 1
            if draw_fl: self.sxyz_fline()
            plot   (-sqrt(2.0)*Si[:ima], -Sj[:ima],  lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            xlabel (r'$-\sqrt{2}\sigma_%s$'%(ikeys[abs(self.isxyz[0])]));  ylabel(r'$-\sigma_%s$'%(ikeys[abs(self.isxyz[1])]));  grid()
            axis   ('equal')

    # Plot octahedral rosette
    # =======================
    def oct_rosette(self, xmin, xmax, ymin, ymax):
        r  = sqrt((xmax-xmin)**2. + (ymax-ymin)**2.)
        cf = 0.2
        cr = 1.1
        l1 = (             0.0  , cr*r            ) # line: 1 end points
        l2 = (-cf*r*cos(pi/6.0) ,-cf*r*sin(pi/6.0)) # line: 2 end points
        l3 = ( cf*r*cos(pi/6.0) ,-cf*r*sin(pi/6.0)) # line: 3 end points
        l4 = (-cr*r*cos(pi/6.0) , cr*r*sin(pi/6.0)) # line: 4 = neg 1 end points
        lo = (-cr*r*cos(pi/3.0) , cr*r*sin(pi/3.0)) # line: origin of cylindrical system
        # main lines
        plot([0.0,l1[0]],[0.0,l1[1]],'k-', color='gray')
        plot([0.0,l2[0]],[0.0,l2[1]],'k-', color='gray')
        plot([0.0,l3[0]],[0.0,l3[1]],'k-', color='gray')
        # reference
        plot([0.0,l4[0]],[0.0, l4[1]],'k-', color='gray')
        plot([0.0,lo[0]],[0.0, lo[1]],'k-', color='gray')
        # text
        text(l1[0],l1[1],r'$-\sigma_1,\theta=+30^\circ$', ha='center', fontsize=8)
        text(l2[0],l2[1],r'$-\sigma_2$',                  ha='right',  fontsize=8)
        text(l3[0],l3[1],r'$-\sigma_3$',                  ha='left',   fontsize=8)
        text(lo[0],lo[1],r'$\theta=0^\circ$',             ha='center', fontsize=8)
        text(l4[0],l4[1],r'$\theta=-30^\circ$',           ha='left',   fontsize=8)

    # Plot failure line in p-q plane
    # ==============================
    def pq_fline(self, pmin, pmax, qmin, qmax):
        Dp    = pmax-pmin
        Dq    = qmax-qmin
        pmin -= 0.5*Dp
        pmax += 0.5*Dp
        qmax += 0.1*Dq
        dp    = (pmax-pmin)/self.fc_np
        dq    = (qmax-qmin)/self.fc_np
        p     = zeros ((self.fc_np,self.fc_np))
        q     = zeros ((self.fc_np,self.fc_np))
        f     = zeros ((self.fc_np,self.fc_np))
        for k in range(len(self.fc_ty)):
            for i in range(self.fc_np):
                for j in range(self.fc_np):
                    p[i,j] = pmin + i*dp
                    q[i,j] = qmin + j*dq
                    s123   = pqt_calc_s123 (p[i,j], q[i,j], self.fc_t, self.pq_ty)
                    sig    = matrix([[s123[0]],[s123[1]],[s123[2]],[0.0]])
                    f[i,j] = self.failure_crit (sig, self.fc_ty[k])
                    if self.fc_t<0.0: q[i,j] = -q[i,j]
            contour (p, q, f, [0.0], linewidths=1, colors=self.fc_clr[k], linestyles=self.fc_ls[k], label=self.fc_ty[k])

    # Plot failure line in octahedral plane
    # =====================================
    def oct_fline(self, samin, samax, sbmin, sbmax):
        r      = sqrt((samax-samin)**2. + (sbmax-sbmin)**2.)
        Dsa    = samax-samin
        Dsb    = sbmax-sbmin
        samin -= 0.1*Dsa + r*cos(pi/6.0)
        samax += 0.1*Dsa + 0.4*r*cos(pi/6.0)
        sbmin -= 0.1*Dsb
        sbmax += 0.1*Dsb
        dsa    = (samax-samin)/self.fc_np
        dsb    = (sbmax-sbmin)/self.fc_np
        sa     = zeros ((self.fc_np,self.fc_np))
        sb     = zeros ((self.fc_np,self.fc_np))
        sc     = self.fc_p
        f      = zeros ((self.fc_np,self.fc_np))
        for k in range(len(self.fc_ty)):
            for i in range(self.fc_np):
                for j in range(self.fc_np):
                    sa[i,j] = samin + i*dsa
                    sb[i,j] = sbmin + j*dsb
                    s123    = oct_calc_s123 (sa[i,j], sb[i,j], sc)
                    sig     = matrix([[s123[0]],[s123[1]],[s123[2]],[0.0]])
                    f[i,j]  = self.failure_crit (sig, self.fc_ty[k])
            contour (sa, sb, f, [0.0], linewidths=1, colors=self.fc_clr[k], linestyles=self.fc_ls[k], label=self.fc_ty[k])

    # Plot failure line in s1-s3, s2-s3 plane
    # =======================================
    def s123_fline(self, s3):
        xmin, xmax = self.ax.get_xbound()
        ymin, ymax = self.ax.get_ybound()
        smax = s3+max([2.2*abs(xmin),2.2*abs(xmax)])
        xmin, xmax = -smax, smax
        ymin, ymax = -smax, smax
        dx  = (xmax-xmin)/self.fc_np
        dy  = (ymax-ymin)/self.fc_np
        x   = zeros ((self.fc_np,self.fc_np))
        y   = zeros ((self.fc_np,self.fc_np))
        f   = zeros ((self.fc_np,self.fc_np))
        for k in range(len(self.fc_ty)):
            for i in range(self.fc_np):
                for j in range(self.fc_np):
                    x[i,j] = xmin + i*dx
                    y[i,j] = ymin + j*dy
                    s1     = s3 - x[i,j]
                    s2     = s3 - y[i,j]
                    sig    = matrix([[s1], [s2], [s3], [0.0]])
                    f[i,j] = self.failure_crit (sig, self.fc_ty[k])
            contour (x, y, f, [0.0], linewidths=1, colors=self.fc_clr[k], linestyles=self.fc_ls[k], label=self.fc_ty[k])

    # Plot failure line in sxyz plane
    # ===============================
    def sxyz_fline(self):
        xmin, xmax = self.ax.get_xbound()
        ymin, ymax = self.ax.get_ybound()
        Dx    = xmax-xmin
        Dy    = ymax-ymin
        DD    = max([Dx,Dy])
        xmin -= 2.0*DD
        xmax += 2.0*DD
        ymax += 0.5*DD
        ymin -= 0.5*DD
        dx    = (xmax-xmin)/self.fc_np
        dy    = (ymax-ymin)/self.fc_np
        x     = zeros ((self.fc_np,self.fc_np))
        y     = zeros ((self.fc_np,self.fc_np))
        f     = zeros ((self.fc_np,self.fc_np))
        sq2   = sqrt(2.0)
        #print xmin, xmax, ymin, ymax
        plot ([xmin,xmax], [xmin/sq2,xmax/sq2], 'k-') # hydrostatic line
        for k in range(len(self.fc_ty)):
            for i in range(self.fc_np):
                for j in range(self.fc_np):
                    x[i,j] = xmin + i*dx
                    y[i,j] = ymin + j*dy
                    sig    = matrix([[-x[i,j]/sq2], [-x[i,j]/sq2], [-y[i,j]], [0.0]])
                    f[i,j] = self.failure_crit (sig, self.fc_ty[k])
            contour (x, y, f, [0.0], linewidths=1, colors=self.fc_clr[k], linestyles=self.fc_ls[k], label=self.fc_ty[k])

    # Failure criterion
    # =================
    def failure_crit(self, sig, fc_ty):
        if fc_ty=='VM':
            p, q = sig_calc_p_q(sig)
            f    = q - 2.0*(sqrt(2.0)/sqrt(3.0))*self.fc_cu
        elif fc_ty=='DP':
            sphi = sin(self.fc_phi*pi/180.0)
            cbar = sqrt(3.0)*self.fc_c/tan(self.fc_phi*pi/180.0)
            kdp  = 2.0*sqrt(2.0)*sphi/(3.0-sphi)
            p, q = sig_calc_p_q(sig)
            f    = q - (p + cbar)*kdp
        elif fc_ty=='MC':
            p, q = sig_calc_p_q (sig)
            t    = sig_calc_t   (sig)
            th   = arcsin(t)/3.0
            sphi = sin(self.fc_phi*pi/180.0)
            cbar = sqrt(3.0)*self.fc_c/tan(self.fc_phi*pi/180.0)
            g    = sqrt(2.0)*sphi/(sqrt(3.0)*cos(th)-sphi*sin(th))
            f    = q - (p + cbar)*g
        elif fc_ty=='MN':
            sphi     = sin(self.fc_phi*pi/180.0)
            cbar     = sqrt(3.0)*self.fc_c/tan(self.fc_phi*pi/180.0)
            kmn      = (9.0-sphi**2.0)/(1.0-sphi**2.0)
            sig0     = cbar*matrix([[1.0],[1.0],[1.0],[0.0]])
            sig_     = sig+sig0
            l        = sig_calc_s123(sig_)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig_)
            f        = I1*I2 - kmn*I3
        elif fc_ty=='LD':
            sphi     = sin(self.fc_phi*pi/180.0)
            cbar     = sqrt(3.0)*self.fc_c/tan(self.fc_phi*pi/180.0)
            kld      = ((3.0-sphi)**3.0)/((1.0+sphi)*((1.0-sphi)**2))
            sig0     = cbar*matrix([[1.0],[1.0],[1.0],[0.0]])
            sig_     = sig+sig0
            l        = sig_calc_s123(sig_)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig_)
            f        = I1**3.0 - kld*I3
        else: raise Exception('failure_crit: fc_ty==%s is invalid' % fc_ty)
        return f

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
        lwd   = 2                  # linewidth
        nhplt = 1                  # number of horizontal plots
        nvplt = 2                  # number of vertical plots
        iplot = 1

        # ux, fx
        self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
        plot   (Ux,Fx,'r-',linewidth=lwd)
        xlabel (r'$u_x$')
        ylabel (r'$f_x$');  grid()

        # uy, fy
        self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
        plot   (Uy,Fy,'r-',linewidth=lwd)
        xlabel (r'$u_y$')
        ylabel (r'$f_y$');  grid()

    # Show figure
    # ===========
    def show(self): show()

    # Set figure for eps
    # ==================
    def set_fig_for_eps(self, small=False, multiplot=False):
        fig_width_pt  = 246.0 if not multiplot else 492.0 # Get this from LaTeX using \showthe\columnwidth
        inches_per_pt = 1.0/72.27                         # Convert pt to inch
        fig_width     = fig_width_pt*inches_per_pt  # width in inches
        fig_height    = fig_width*self.proport      # height in inches
        fig_size      = [fig_width,fig_height]
        if small:
            params = {'backend':         'ps',
                      'axes.labelsize':  8,
                      'text.fontsize':   8,
                      'legend.fontsize': 8,
                      'xtick.labelsize': 6,
                      'ytick.labelsize': 6,
                      'text.usetex':     True,
                      'figure.figsize': fig_size}
        else:
            params = {'backend':         'ps',
                      'axes.labelsize':  10,
                      'text.fontsize':   10,
                      'legend.fontsize': 10,
                      'xtick.labelsize': 8,
                      'ytick.labelsize': 8,
                      'text.usetex':     True,
                      'figure.figsize': fig_size}
        rcParams.update(params)

    # FC legend
    # =========
    def fc_leg(self, loc='center right'):
        lines = []
        for k in range(len(self.fc_ty)):
            lines.append (plot([0],[0],linestyle=self.fc_ls[k],color=self.fc_clr[k]))
        legend (lines, self.fc_ty, loc=loc)
