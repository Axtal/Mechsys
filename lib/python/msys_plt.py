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

from msys_invs import *
from msys_fig  import *

class Plotter:
    # Constructor
    # ===========
    def __init__(self):
        # data
        self.show_k     = False     # show k=dq/dp ?
        self.div_by_p   = False     # divide q by p ?
        self.log_p      = True      # use log(p) instead of p ?
        self.lnplus1    = False     # plot ln(p+1) instead of ln(p) ?
        self.q_neg_ext  = True      # multiply q by -1 for extension (t<0, where t=sin(3th)
        self.ed_neg_ext = False     # multiply ed by -1 for extension (t<0, where t=sin(3th)
        self.pq_ty      = 'cam'     # invariants type
        self.evd_ty     = 'cam'     # strain invariants type
        self.mark_max   = False     # mark max (failure) point ?
        self.mark_lst   = False     # mark residual (failure) point ?
        self.oct_norm   = False     # normalize plot in octahedral plane by p ?
        self.oct_sxyz   = True      # use Sx,Sy,Sz in oct plane instead of S1,S2,S3
        self.isxyz      = (1,0)     # indices for sxyz plot, use negative numbers for principal components
        self.devplot    = True      # plot s3-s1, s3-s2 instead of Ek, Sk
        self.pcte       = -1        # if pcte>0 => pcte in Ev x p (logp) plot?
        self.maxed      = -1        # max Ed to stop the lines
        self.maxev      = -1        # max Ev to stop the lines
        self.one        = -1        # all plots = -1
        self.six        = False     # only 2 x 3 instead of 3 x 3 plots?
        self.four       = False     # only 2 x 4 instead of others


    # Plot results
    # ============
    def plot(self, filename, clr='red', lwd=2, label=None, zorder=None,
             marker='None', markevery=None, ms=rcParams['lines.markersize']):

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
            S1[i], S2[i], S3[i] = s123  [0], s123  [1], s123  [2]
            Sx[i], Sy[i], Sz[i] = Sig[i][0], Sig[i][1], Sig[i][2]
            sxyz                = Sx[i], Sy[i], Sz[i]
            Sa[i], Sb[i], Sc[i] = sxyz_calc_oct  (sxyz) if self.oct_sxyz else s123_calc_oct (s123)
            Ev[i], Ed[i]        = eps_calc_ev_ed (Eps[i], Type=self.evd_ty)
            e123                = eps_calc_e123  (Eps[i])
            E1[i], E2[i], E3[i] = e123  [0], e123  [1], e123  [2]
            Ex[i], Ey[i], Ez[i] = Eps[i][0], Eps[i][1], Eps[i][2]
            if self.q_neg_ext  and T[i]<0.0: Q [i] = -Q [i]
            if self.ed_neg_ext and T[i]<0.0: Ed[i] = -Ed[i]
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
        nhplt = 3  # number of horizontal plots
        nvplt = 3  # number of vertical plots
        iplot = 1  # index of plot
        if self.six:  nhplt, nvplt = 2, 3
        if self.four: nhplt, nvplt = 2, 2

        # calc friction angle
        imaQP   = QdivP.argmax() # index of QP used to calculate phi
        fc_phi  = M_calc_phi (QdivP[imaQP], self.pq_ty)
        fc_poct = P[imaQP]*sqrt(3.0) if self.pq_ty=='cam' else P[imaQP]
        fc_cu   = qf_calc_cu (Q[imaQP], self.pq_ty)

        # q p ratio and label
        Y = QdivP if self.div_by_p else Q
        Ylbl = r'$q_{%s}/p_{%s}$'%(self.pq_ty,self.pq_ty) if self.div_by_p else r'$q_{%s}$'%(self.pq_ty)

        # 0) q/p, Ed ---------------------------------------------------------------------------
        if self.one==0 or self.one<0:
            if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot (Ed, Y, color=clr, lw=lwd, label=label, marker=marker, markevery=markevery, ms=ms)
            if self.mark_max: plot (Ed[imaQP], Y[imaQP], '^', color=clr)
            if self.mark_lst: plot (Ed[-1],    Y[-1],    '^', color=clr)
            xlabel (r'$\varepsilon_d$ [%]');  ylabel(Ylbl);  grid(True)

        # 1) q/p, Ev ---------------------------------------------------------------------------
        if self.one==1 or self.one<0:
            if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ev, Y, lw=lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=ms)
            xlabel (r'$\varepsilon_v$ [%]');  ylabel(Ylbl);  grid(True)

        # 2) p, q ---------------------------------------------------------------------------
        if self.one==2 or self.one<0 and not self.four:
            if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            axhline (0.0,color='black'); axvline(0.0,color='black')
            plot    (P, Q, lw=lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=ms)
            xlabel  (r'$p_{%s}$'%(self.pq_ty));  ylabel(r'$q_{%s}$'%(self.pq_ty));  grid(True)
            axis    ('equal')
            if self.show_k:
                k = (Q[-1]-Q[0])/(P[-1]-P[0])
                text ((P[0]+P[-1])/2.0,(Q[0]+Q[-1])/2.0,'k = %g'%k,color='black',ha='left', fontsize=10)

        # 3) Ed, Ev ---------------------------------------------------------------------------
        if self.one==3 or self.one<0:
            if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (Ed, Ev, lw=lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=ms)
            if self.mark_lst: plot (Ed[-1], Ev[-1], '^', color=clr)
            xlabel (r'$\varepsilon_d$ [%]');  ylabel(r'$\varepsilon_v$ [%]'); grid(True)

        # 4) lnp, Ev ---------------------------------------------------------------------------
        if self.one==4 or self.one<0:
            if self.log_p:
                if self.lnplus1:
                    X    = log(P+1.0)
                    xlbl = r'$\ln{(1+p_{%s})}$'%(self.pq_ty)
                else:
                    X    = log(P)
                    xlbl = r'$\ln{(p_{%s})}$'%(self.pq_ty)
            else:
                X    = P
                xlbl = r'$p_{%s}$'%(self.pq_ty)
            if self.pcte>0:
                for k, x in enumerate(X): X[k] = self.pcte
            if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (X, Ev, lw=lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=ms)
            xlabel (xlbl);  ylabel(r'$\varepsilon_v$ [%]');  grid(True)

        # 5) Sa, Sb ---------------------------------------------------------------------------
        if self.one==5 or self.one<0 and not self.four:
            pcoef = self.fc_poct if self.oct_norm else 1.0
            if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot (Sa/pcoef, Sb/pcoef, color=clr, lw=lwd, label=label, marker=marker, markevery=markevery, ms=ms)
            if self.mark_max: plot (Sa[imaQP]/pcoef, Sb[imaQP]/pcoef, '^', color=clr)
            if self.mark_lst: plot (Sa[-1   ]/pcoef,  Sb[-1  ]/pcoef, '^', color=clr)
            axis ('equal')

        # 6) Ek, Q/P ---------------------------------------------------------------------------
        if self.one==6 or self.one<0 and not self.six and not self.four:
            if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
            plot   (E1, Y, lw=lwd,linestyle='-',  color=clr)
            plot   (E2, Y, lw=lwd,linestyle='--', color=clr)
            plot   (E3, Y, lw=lwd,linestyle='-.', color=clr)
            xlabel (r'$\varepsilon_1$[--], $\varepsilon_2$[- -], $\varepsilon_3$[- .]')
            ylabel (Ylbl);  grid(True)

        if self.devplot:
            if self.one==7 or self.one<0 and not self.six and not self.four:
                # 7) s3-s1, s3-s2
                if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
                plot   (S3-S1, S3-S2, lw=lwd,linestyle='-', color=clr, label=label, marker=marker, markevery=markevery, ms=ms)
                xlabel (r'$\sigma_3-\sigma_1$');  ylabel(r'$\sigma_3-\sigma_2$');  grid(True)
                axis   ('equal')
        else:
            # 7) Ek, Sk ---------------------------------------------------------------------------
            if self.one==7 or self.one<0 and not self.six and not self.four:
                if self.one<0: self.ax = subplot(nhplt,nvplt,iplot);  iplot += 1
                plot   (Ex, -Sx, lw=lwd,linestyle='-', color=clr)
                plot   (Ey, -Sy, lw=lwd,linestyle='--', color=clr)
                plot   (Ez, -Sz, lw=lwd,linestyle='-.', color=clr)
                xlabel (r'$\varepsilon_x$[--], $\varepsilon_y$[- -], $\varepsilon_z$[- .]')
                ylabel (r'$-\sigma_x$[--], $-\sigma_y$[- -], $-\sigma_z$[- .]');  grid(True)

        # 8) sqrt(2.0)*Si, Sj ---------------------------------------------------------------------------
        if self.one==8 or self.one<0 and not self.six and not self.four:
            if self.one<0: self.ax = subplot (nhplt,nvplt,iplot);  iplot += 1
            plot   (-sqrt(2.0)*Si, -Sj,  lw=lwd, color=clr, label=label, marker=marker, markevery=markevery, ms=ms)
            xlabel (r'$-\sqrt{2}\sigma_%s$'%(ikeys[abs(self.isxyz[0])]));  ylabel(r'$-\sigma_%s$'%(ikeys[abs(self.isxyz[1])]));  grid(True)
            axis   ('equal')

        # 9) Mohr-circles -----------------------------------------------------------------------------
        if self.one==9:
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
                Arc (C0,0.,R0, 0.,pi, ec=clr, res=30)
                Arc (C1,0.,R1, 0.,pi, ec=clr, res=30)
                Arc (C2,0.,R2, 0.,pi, ec=clr, res=30)
            xlabel (r'$-\sigma_i$');  ylabel(r'$\tau$');  grid(True)
            axis   ('equal')

        # return failure criterion constants
        return fc_phi, fc_poct, fc_cu


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