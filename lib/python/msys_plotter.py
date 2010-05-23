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

import os, subprocess
from os.path import basename
from numpy   import exp, sqrt, matrix, zeros, log, cos, pi, sin, tan, linspace, polyfit
from pylab   import rc, subplot, plot, xlabel, ylabel, grid, axhline, axvline, axis, text, contour, show
from pylab   import rcParams, savefig, axes, legend, gca, title, figure, clf, annotate
from pylab   import matplotlib as MPL
from msys_invariants import *
from msys_readdata   import *
from msys_linfit     import *

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
        self.fc_ty     = ['MC','MN']                              # failure criteria type (VM, DP, MC, MN)
        self.fc_clr    = ['k', 'k', 'k', 'k']                     # failure criteria colors
        self.fc_ls     = ['solid', 'dashed', 'dotted', 'dashdot'] # failure criteria linestyles
        self.fc_poct   = 100.0*sqrt(3.0)                          # p(oct) to be used when plotting FC in octahedral plane
        self.fc_t      = 1.0                                      # t=sin(3th) to be used when plotting FC in octahedral plane
        self.fc_cu     = 10.0                                     # cohesion for VM
        self.fc_phi    = -1                                       # friction angle for FC
        self.fc_c      = 0.0                                      # cohesion for FC
        self.fc_np     = 40                                       # number of points for drawing failure line
        self.oct_norm  = False                                    # normalize plot in octahedral plane by p ?
        self.rst_phi   = True                                     # show fc_phi in Rosette
        self.rst_circ  = True                                     # draw circle in Rosette
        self.isxyz     = (-1,0)                                   # indices for sxyz plot, use negative numbers for principal components
        self.devplot   = True                                     # plot s3-s1, s3-s2 instead of Ek, Sk
        self.pcte      = -1                                       # if pcte>0 => pcte in Ev x p (logp) plot?
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
        self.fc_prms   = {'A':None,'B':None,'c':None,'bet':None}  # nonlinear FC parameters

        # matplotlib's structures
        self.PH = MPL.path.Path
        self.PP = MPL.patches
        self.PC = MPL.patches.PathPatch
        # colors
        self.lgray = (180/255.,180/255.,180/255.)

    # Plot results
    # ============
    def plot(self, filename, clr='red', draw_fl=False, draw_ros=False, txtlst=False, txtmax=False,
                             label=None, marker='None', markevery=None, calc_phi=False):
        # load data
        Sig, Eps = self.load_data (filename)

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
            s123                = sig_calc_s123  (Sig[i], do_sort=False)
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
        QdivP = Q/P

        # constants
        rc('text', usetex=True)    # set LaTeX
        rc('font', family='serif') # set font
        nhplt = 3                  # number of horizontal plots
        nvplt = 3                  # number of vertical plots
        iplot = 1                  # index of plot
        if self.only_six: nhplt, nvplt = 2, 3

        # calc friction angle
        imaQP = QdivP.argmax() # index of QP used to calculate phi
        if calc_phi or self.fc_phi<0:
            self.fc_phi  = M_calc_phi (QdivP[imaQP], self.pq_ty)
            self.fc_poct = P[imaQP]*sqrt(3.0) if self.pq_ty=='cam' else P[imaQP]
            self.fc_cu   = qf_calc_cu (Q[imaQP], self.pq_ty)

        # set for eps
        if self.set_eps:
            if self.justone>=0:
                self.set_fig_for_eps(multiplot=False)
                #if self.justone==5: axes([0.01,0.01,.99,.99])
                #else: axes(self.axesdat) # this needs to be after set_eps
            else:
                self.set_fig_for_eps(multiplot=True)
                #axes([0.,0.,0.99,0.99]) # this needs to be after set_eps

        # q p ratio and label
        Y = QdivP if self.div_by_p else Q
        Ylbl = r'$q_{%s}/p_{%s}$'%(self.pq_ty,self.pq_ty) if self.div_by_p else r'$q_{%s}$'%(self.pq_ty)

        # 0) q/p, Ed ---------------------------------------------------------------------------
        if self.justone==0 or self.justone<0:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot (Ed, Y, color=clr, lw=self.lwd, label=label, marker=marker, markevery=markevery, ms=self.ms)
            plot (Ed[imaQP], Y[imaQP], '^', color=clr)
            #plot (Ed[-1],    Y[-1],    '^', color=clr)
            #xlabel (r'$\varepsilon_d$ [\%]');  ylabel(Ylbl);  grid()
            if txtmax: text (Ed[imaQP], Y[imaQP], '%.2f'%Y[imaQP], fontsize=8)
            if txtlst: text (Ed[-1],    Y[-1],    '%.2f'%Y[-1],    fontsize=8)

        # 1) q/p, Ev ---------------------------------------------------------------------------
        if self.justone==1 or self.justone<0:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ev, Y, lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            #xlabel (r'$\varepsilon_v$ [\%]');  ylabel(Ylbl);  grid()

        # 2) p, q ---------------------------------------------------------------------------
        if self.justone==2 or self.justone<0:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            axhline (0.0,color='black'); axvline(0.0,color='black')
            if draw_fl: self.pq_fline(0.1*min(P),2.0*max(P),min(Q),max(Q))
            plot    (P, Q, lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            #xlabel  (r'$p_{%s}$'%(self.pq_ty));  ylabel(r'$q_{%s}$'%(self.pq_ty));  grid()
            axis    ('equal')
            if self.show_k:
                k = (Q[-1]-Q[0])/(P[-1]-P[0])
                text ((P[0]+P[-1])/2.0,(Q[0]+Q[-1])/2.0,'k = %g'%k,color='black',ha='left', fontsize=10)

        # 3) Ed, Ev ---------------------------------------------------------------------------
        if self.justone==3 or self.justone<0:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ed, Ev, lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            #xlabel (r'$\varepsilon_d$ [\%]');  ylabel(r'$\varepsilon_v$ [\%]'); grid()

        # 4) lnp, Ev ---------------------------------------------------------------------------
        if self.justone==4 or self.justone<0:
            if self.log_p:
                X    = log(P)
                xlbl = r'$\ln{(p_{%s})}$'%(self.pq_ty)
            else:
                X    = P
                xlbl = r'$p_{%s}$'%(self.pq_ty)
            if self.pcte>0:
                for k, x in enumerate(X): X[k] = self.pcte
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (X, Ev, lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            #xlabel (xlbl);  ylabel(r'$\varepsilon_v$ [\%]');  grid()

        # 5) Sa, Sb ---------------------------------------------------------------------------
        if self.justone==5 or self.justone<0:
            pcoef = self.fc_poct if self.oct_norm else 1.0
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            if draw_ros: self.oct_rosette(min(Sa)/pcoef,max(Sa)/pcoef,min(Sb)/pcoef,max(Sb)/pcoef)
            if draw_fl:  self.oct_fline  (min(Sa)/pcoef,max(Sa)/pcoef,min(Sb)/pcoef,max(Sb)/pcoef)
            plot (Sa/pcoef, Sb/pcoef, color=clr, lw=self.lwd, label=label, marker=marker, markevery=markevery, ms=self.ms)
            plot (Sa[imaQP]/pcoef, Sb[imaQP]/pcoef, '^', color=clr)
            #plot (Sa[-1   ]/pcoef,  Sb[-1  ]/pcoef, '^', color=clr)
            axis ('equal')
            axis ('off')

        # 6) Ek, Q/P ---------------------------------------------------------------------------
        if self.justone==6 or self.justone<0 and not self.only_six:
            if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (E1, Y, lw=self.lwd,linestyle='-',  color=clr)
            plot   (E2, Y, lw=self.lwd,linestyle='--', color=clr)
            plot   (E3, Y, lw=self.lwd,linestyle='-.', color=clr)
            #xlabel (r'$\varepsilon_1$[--], $\varepsilon_2$[- -], $\varepsilon_3$[- .]')
            ylabel (Ylbl);  grid()

        if self.devplot:
            if self.justone==7 or self.justone<0 and not self.only_six:
                # 7) s3-s1, s3-s2
                if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
                if draw_fl: self.s123_fline(S3[imaQP])
                plot   (S3-S1, S3-S2, lw=self.lwd,linestyle='-', color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
                #xlabel (r'$\sigma_3-\sigma_1$');  ylabel(r'$\sigma_3-\sigma_2$');  grid()
                axis   ('equal')
        else:
            # 7) Ek, Sk ---------------------------------------------------------------------------
            if self.justone==7 or self.justone<0 and not self.only_six:
                if self.justone<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
                plot   (Ex, -Sx, lw=self.lwd,linestyle='-', color=clr)
                plot   (Ey, -Sy, lw=self.lwd,linestyle='--', color=clr)
                plot   (Ez, -Sz, lw=self.lwd,linestyle='-.', color=clr)
                #xlabel (r'$\varepsilon_x$[--], $\varepsilon_y$[- -], $\varepsilon_z$[- .]')
                ylabel (r'$-\sigma_x$[--], $-\sigma_y$[- -], $-\sigma_z$[- .]');  grid()

        # 8) sqrt(2.0)*Si, Sj ---------------------------------------------------------------------------
        if self.justone==8 or self.justone<0 and not self.only_six:
            if self.justone<0: self.ax = subplot (nhplt,nvplt,iplot);  iplot += 1
            if draw_fl: self.sxyz_fline()
            plot   (-sqrt(2.0)*Si, -Sj,  lw=self.lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=self.ms)
            #xlabel (r'$-\sqrt{2}\sigma_%s$'%(ikeys[abs(self.isxyz[0])]));  ylabel(r'$-\sigma_%s$'%(ikeys[abs(self.isxyz[1])]));  grid()
            axis   ('equal')

        # 9) Mohr-circles -----------------------------------------------------------------------------
        if self.justone==9:
            #each = 5
            #for k in range(0,np,each):
            max_s = max([max(-S1), max(-S2), max(-S3)])
            for k in [np-2,np-1]:
                s1 = -S1[k]
                s2 = -S2[k]
                s3 = -S3[k]
                C0 =     (s1+s2)/2.
                C1 =     (s2+s3)/2.
                C2 =     (s3+s1)/2.
                R0 = abs((s1-s2)/2.)
                R1 = abs((s2-s3)/2.)
                R2 = abs((s3-s1)/2.)
                self.draw_arc(gca(), C0,0.,R0, 0., pi, clr, res=30)
                self.draw_arc(gca(), C1,0.,R1, 0., pi, clr, res=30)
                self.draw_arc(gca(), C2,0.,R2, 0., pi, clr, res=30)
            p1, = plot ([0], [0], 'r-')
            p2, = plot ([0,max_s],[0,max_s*tan(self.fc_phi*pi/180.)],'-',color='orange')
            #xlabel (r'$-\sigma_i$');  ylabel(r'$\tau$');  grid()
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
        plot([0.0,l1[0]],[0.0,l1[1]],'k-', color=self.lgray)
        plot([0.0,l2[0]],[0.0,l2[1]],'k-', color=self.lgray)
        plot([0.0,l3[0]],[0.0,l3[1]],'k-', color=self.lgray)
        # reference
        plot([0.0,l4[0]],[0.0, l4[1]],'k-', color=self.lgray)
        plot([0.0,lo[0]],[0.0, lo[1]],'k-', color=self.lgray)
        # text
        txt = '/p_{oct}' if self.oct_norm else ''
        text(l1[0],l1[1],r'$-\sigma_1%s,\theta=+30^\circ$'%txt, ha='center', fontsize=8)
        text(l2[0],l2[1],r'$-\sigma_2%s$'%txt,                  ha='right',  fontsize=8)
        text(l3[0],l3[1],r'$-\sigma_3%s$'%txt,                  ha='left',   fontsize=8)
        text(lo[0],lo[1],r'$\theta=0^\circ$',                   ha='center', fontsize=8)
        text(l4[0],l4[1],r'$\theta=-30^\circ$',                 ha='left',   fontsize=8)
        if self.rst_phi:
            text(0.0,l2[1]-0.05*r,r'$\phi_{comp}=%2.1f^\circ$'%self.fc_phi, ha='center', va='top', fontsize=10)
        if self.rst_circ:
            M = phi_calc_M (self.fc_phi, 'oct')
            R = M if self.oct_norm else self.fc_poct*M
            self.draw_arc(gca(), 0.,0.,R, 0.85*pi/2., 1.05*(pi/2.+pi/3.), self.lgray)

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
        samin -= 0.01*Dsa +     r*cos(pi/6.0)
        samax += 0.01*Dsa + 0.3*r*cos(pi/6.0)
        sbmin -= 0.01*Dsb
        sbmax += 0.01*Dsb + 0.1*r
        dsa    = (samax-samin)/self.fc_np
        dsb    = (sbmax-sbmin)/self.fc_np
        sa     = zeros ((self.fc_np,self.fc_np))
        sb     = zeros ((self.fc_np,self.fc_np))
        sc     = 1.0 if self.oct_norm else self.fc_poct
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
            #sphi     = sin(self.fc_phi*pi/180.0)
            #cbar     = sqrt(3.0)*self.fc_c/tan(self.fc_phi*pi/180.0)
            #kmn      = (9.0-sphi**2.0)/(1.0-sphi**2.0)
            #sig0     = cbar*matrix([[1.0],[1.0],[1.0],[0.0]])
            #sig_     = sig+sig0
            tgphi2   = tan(self.fc_phi*pi/180.0)**2.0
            kmn      = 9.0 + 8.0*tgphi2
            sig_     = sig
            l        = sig_calc_s123(sig_)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig_)
            f        = I1*I2 - kmn*I3
        elif fc_ty=='MNnl':
            p, q = sig_calc_p_q (sig, 'cam')
            M    = self.refcurve(p)/p if p>0.0 else self.fc_prms['A']
            sphi = 3.0*M/(M+6.0)
            kmn  = (9.0-sphi**2.0)/(1.0-sphi**2.0)
            l    = sig_calc_s123(sig)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig)
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
    def fc_leg(self, loc='best'):
        lines = []
        for k in range(len(self.fc_ty)):
            lines.append (plot([0],[0],linestyle=self.fc_ls[k],color=self.fc_clr[k]))
        legend (lines, self.fc_ty, loc=loc, prop={'size':9})

    # Draw arc
    # ========
    def draw_arc(self, ax, xc,yc,R, alp_min,alp_max, eclr, fclr='None', lwd=1, res=200):
        A   = linspace(alp_min,alp_max,res)
        dat = [(self.PH.MOVETO, (xc+R*cos(alp_min),yc+R*sin(alp_min)))]
        for a in A[1:]:
            x = xc + R*cos(a)
            y = yc + R*sin(a)
            dat.append((self.PH.LINETO, (x,y)))
        cmd,vert = zip(*dat)
        ph = self.PH (vert, cmd)
        pc = self.PC (ph, facecolor=fclr, edgecolor=eclr, linewidth=lwd)
        ax.add_patch (pc)

    # Ref curve model
    # ===============
    def refcurve(self, x):
        A   = self.fc_prms['A']
        B   = self.fc_prms['B']
        c   = self.fc_prms['c']
        bet = self.fc_prms['bet']
        c1  = bet*(A-B)
        c2  = exp(-c*bet)
        c3  = 1.0-c2
        return A*x - log(c3+c2*exp(c1*x))/bet

    # Find phi
    # ========
    def find_phi(self, files, with_refcurve=False, find_refcte=False, mcirc=False, txt='comp'):
        # load data and calculate additional variables
        phi_ave = 0.0
        p_at_qpmax, q_at_qpmax                = [], []
        s1_at_qpmax, s2_at_qpmax, s3_at_qpmax = [], [], []
        for f in files:
            Sig, Eps   = self.load_data (f)
            np         = len(Sig) # number of points
            P,  Q      = zeros(np), zeros(np)
            S1, S2, S3 = zeros(np), zeros(np), zeros(np) # principal stresses
            for i in range(np):
                P[i], Q[i]          = sig_calc_p_q  (Sig[i], Type='cam')
                s123                = sig_calc_s123 (Sig[i], do_sort=True)
                S1[i], S2[i], S3[i] = s123[0], s123[1], s123[2]
            QdivP    = Q/P
            imaQP    = QdivP.argmax()
            phi_ave += M_calc_phi (QdivP[imaQP], 'cam')
            q_at_qpmax .append (Q [imaQP])
            p_at_qpmax .append (P [imaQP])
            s1_at_qpmax.append (S1[imaQP])
            s2_at_qpmax.append (S2[imaQP])
            s3_at_qpmax.append (S3[imaQP])

        if with_refcurve:
            if find_refcte:
                X, Y = array(p_at_qpmax[:2]), array(q_at_qpmax[:2])
                f0   = LinFit(X,Y, tls=False, cmx=False)
                self.fc_prms['A'] = f0.m
                X, Y = array(p_at_qpmax[-2:]), array(q_at_qpmax[-2:])
                f1   = LinFit(X,Y, tls=False, cmx=True)
                self.fc_prms['B']   = f1.m
                self.fc_prms['c']   = f1.c
                self.fc_prms['bet'] = 1.0
            phi_ini = M_calc_phi (self.fc_prms['A'],'cam')
            phi_fin = M_calc_phi (self.fc_prms['B'],'cam')

        # phi fit
        X, Y     = array(p_at_qpmax), array(q_at_qpmax)
        f        = LinFit (X, Y, tls=True, cmx=False)
        M        = f.m
        phi_fit  = M_calc_phi (M,'cam')
        phi_ave /= len(p_at_qpmax)

        # plot
        if mcirc: x = linspace(0., -min([min(s1_at_qpmax),min(s2_at_qpmax),min(s3_at_qpmax)]), 100)
        else:     x = linspace(0., max(X), 100)
        self.proport = 1.0
        self.set_fig_for_eps()
        if with_refcurve: axes([0.12,0.12,0.85,0.75])
        else:             axes([0.12,0.12,0.85,0.85])
        if mcirc: xlabel(r'$-\sigma_{i}$'); ylabel(r'$\tau$')
        else:     xlabel(r'$p_{cam}$'); ylabel(r'$q_{cam}$')
        grid ()
        if mcirc: axis('equal')
        if with_refcurve and not mcirc:
            p2, = plot (x, self.refcurve(x),'b-', marker='.', markevery=10, linewidth=1)
            p3, = plot ([0],[0],'b-',linewidth=2)
            lb  = legend([p2], [r'$A=%2.2f, B=%2.2f, c=%2.2f, \beta=%2.2f$'%(self.fc_prms['A'],self.fc_prms['B'],self.fc_prms['c'],self.fc_prms['bet'])],
                         bbox_to_anchor=(0,1.03,1,0.1), loc=3, mode='expand', borderaxespad=0.)
            lc  = legend([p3], [r'$\phi_{ini}=%2.2f^\circ, \phi_{fin}=%2.2f^\circ$'%(phi_ini,phi_fin)], loc='lower right')
        alp = tan(phi_fit*pi/180.) if mcirc else M
        p0, = plot (x, alp*x, 'g-', marker='None', markevery=10, linewidth=1)
        if mcirc:
            for k, s1 in enumerate(s1_at_qpmax):
                s1 = -s1
                #s2 = -s2_at_qpmax[k]
                s3 = -s3_at_qpmax[k]
                #C0 =     (s1+s2)/2.
                #C1 =     (s2+s3)/2.
                C2 =     (s3+s1)/2.
                #R0 = abs((s1-s2)/2.)
                #R1 = abs((s2-s3)/2.)
                R2 = abs((s3-s1)/2.)
                #self.draw_arc(gca(), C0,0.,R0, 0., pi, 'red', res=30)
                #self.draw_arc(gca(), C1,0.,R1, 0., pi, 'red', res=30)
                self.draw_arc(gca(), C2,0.,R2, 0., pi, 'red', res=30)
            p1, = plot ([0], [0], 'r-')
        else: p1, = plot (X, Y, 'r^', clip_on=False)
        la  = legend([p0,p1], [r'fit: $\phi=%2.2f^\circ$'%(phi_fit), 'DEM data'], loc='upper left')
        if with_refcurve:
            gca().add_artist(lb)
            gca().add_artist(lc)

        return phi_fit, phi_ave


    # Load Data
    # =========
    def load_data(self, filename):
        dat = read_table(filename)
        Sig = []
        Eps = []
        sq2 = sqrt(2.0)
        if dat.has_key('Sx'): # old data file
            raise Exception('load_data: old data file: with Sx: not ready yet')
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
        elif dat.has_key('Sa'): # old data file
            for i in range(len(dat['Sa'])):
                Sig.append(matrix([[-float( dat['Sr' ][i] )],
                                   [-float( dat['St' ][i] )],
                                   [-float( dat['Sa' ][i] )],
                                   [-float( dat['Sat'][i] if dat.has_key('Sat') else 0.0)*sq2]]))
                Eps.append(matrix([[-float( dat['Er' ][i] )/100.],
                                   [-float( dat['Et' ][i] )/100.],
                                   [-float( dat['Ea' ][i] )/100.],
                                   [-float( dat['Eat'][i] if dat.has_key('Eat') else 0.0)*sq2/100.]]))
                Ev, Ed = eps_calc_ev_ed (Eps[len(Eps)-1])
                Ev *= 100.0 # convert strains to percentage
                Ed *= 100.0
                if self.maxed>0 and Ed>self.maxed: break
                if self.maxev>0 and Ev>self.maxev: break
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
                Ev, Ed = eps_calc_ev_ed (Eps[len(Eps)-1])
                Ev *= 100.0 # convert strains to percentage
                Ed *= 100.0
                if self.maxed>0 and Ed>self.maxed: break
                if self.maxev>0 and Ev>self.maxev: break
        return Sig, Eps


    # Save and crop EPS
    # =================
    def save_eps(self, filekey):
        # save
        savefig (filekey+'.eps', bbox_inches='tight')

        # crop
        subprocess.check_call (['ps2eps', '-q', '-l', '-f', '%s.eps'%filekey])

        # rename
        os.rename ('%s.eps.eps'%filekey, '%s.eps'%filekey)

        print "<[1;34m%s.eps[0m> created"%filekey

    # Plot node
    # =========
    def plot_node(self, filename):
        dat = read_table(filename)

        # constants
        rc('text', usetex=True)               # set LaTeX
        rc('font', family='serif')            # set font
        lwd   = 2                             # linewidth
        nhplt = 1                             # number of horizontal plots
        nvplt = 3 if dat.has_key('fz') else 2 # number of vertical plots
        iplot = 1

        # ux, fx
        self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
        plot   (dat['ux'],dat['fx'],'r-',linewidth=lwd)
        xlabel (r'$u_x$')
        ylabel (r'$f_x$');  grid()

        # uy, fy
        self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
        plot   (dat['uy'],dat['fy'],'r-',linewidth=lwd)
        xlabel (r'$u_y$')
        ylabel (r'$f_y$');  grid()

        # uz, fz
        if dat.has_key('fz'):
            self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (dat['uz'],dat['fz'],'r-',linewidth=lwd)
            xlabel (r'$u_z$')
            ylabel (r'$f_z$');  grid()
