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

from msys_invariants import *
from msys_fig        import *
from numpy           import ogrid

class FCrits:
    # Constructor
    # ===========
    def __init__(self):
        # data
        self.c       = 0.0                   # cohesion
        self.phi     = 30.0                  # friction angle for FC
        self.b       = None                  # anisotropic fc: b coefficient
        self.alp     = None                  # anisotropic fc: alpha coefficient
        self.a       = None                  # anisotropic fc: bedding planes normal
        self.sstar   = None                  # anisotropic fc: sig to pass fc through
        self.R       = None                  # anisotropic fc: R coefficient
        self.obliq1  = False                 # anisotropic fc: oblique projection type 1
        self.obliq2  = False                 # anisotropic fc: oblique projection type 2
        self.Imat    = matrix(diag(ones(3))) # identity matrix 3x3
        self.r       = None                  # radius in Pi plane
        self.sc      = 1.0                   # mean pressure => distance of cross-section to the origin of Haigh-Westergaard space
        self.samin   = None                  # min sa in Pi plane
        self.samax   = None                  # max sa in Pi plane
        self.sbmin   = None                  # min sb in Pi plane
        self.sbmax   = None                  # max sb in Pi plane
        self.scmin   = None                  # min sc in Pi plane (3D only)
        self.scmax   = None                  # max sc in Pi plane (3D only)
        self.sxyz    = True                  # consider sx(right), sy(left), sz(up) instead of s1(up), s2(right), s3(left)
        self.rst_txt = []                    # rosette text
        self.leg_plt = []                    # legend plots (items appended in plot)
        self.leg_txt = []                    # legend texts (items appended in plot)


    # Intersection and friction angle
    # ===============================
    def inter (self, sig_star=2, do_plot=True):
        sphi = sin(self.phi*pi/180.0)
        A    = (1.0+sphi)/(1.0-sphi)
        if   sig_star==0: l = array([-A , -1., -1.])
        elif sig_star==1: l = array([-1., -A , -1.])
        elif sig_star==2: l = array([-1., -1., -A ])
        dev_sig_tr = array([2.*l[0]-l[1]-l[2],
                            2.*l[1]-l[2]-l[0],
                            2.*l[2]-l[0]-l[1]])/3.0
        sig   = array([-1.,-1.,-1.])/sqrt(3.)
        a     = 0.5
        sig_a = sig + a*dev_sig_tr
        if do_plot:
            sa, sb, sc = sxyz_calc_oct (sig_a)
            print sa, sb, sc
            plot ([sa], [sb], 'go')




    # Set constants of anisotropic criterion
    # ======================================
    def aniso (self, phi_deg=30.0, b=None, alpha=0.1, a=[0.,0.,1.], sig_star=2):

        # set constants
        a = matrix(a).T
        self.phi   = phi_deg
        self.b     = b
        self.alp   = alpha
        self.a     = a / norm(a)
        self.sstar = sig_star

        # assemble l, vector with principal values
        sphi = sin(self.phi*pi/180.0)
        A    = (1.0+sphi)/(1.0-sphi)
        if   sig_star==0: l = array([-A , -1., -1.])
        elif sig_star==1: l = array([-1., -A , -1.])
        elif sig_star==2: l = array([-1., -1., -A ])
        else: raise Exception('set_anisocrit: sig_star must be 0, 1, or 2. %d is invalid'%sig_star)

        # invariants
        nvec, namp, sig, tau = self.get_nvec_namp_sig_tau (l)
        self.R               = tau/sig


    # Anisotropic criterion invariants
    # ================================
    def get_nvec_namp_sig_tau(self, l, Q=None):

        # calculate nvec
        if self.b==None:
            I2   = l[0]*l[1] + l[1]*l[2] + l[2]*l[0]
            I3   = l[0]*l[1]*l[2]
            nvec = -matrix(sqrt(I3/(I2*l))).T
        else:
            nnew = -matrix([[abs(l[0])**(-self.b)], [abs(l[1])**(-self.b)], [abs(l[2])**(-self.b)]])
            nvec = nnew / norm(nnew)

        # calculate namp
        if Q==None: namp = nvec + self.alp *       self.a
        else:       namp = nvec + self.alp * Q.T * self.a
        namp = namp / norm(namp)

        # calculate sig and tau
        if   self.obliq1: Pamp = (nvec * namp.T) / (nvec.T * namp)[0]
        elif self.obliq2: Pamp = (namp * nvec.T) / (namp.T * nvec)[0]
        else:             Pamp =  namp * namp.T
        Qamp = self.Imat - Pamp
        tamp = diag(l) * namp
        pamp = Pamp*tamp
        qamp = Qamp*tamp
        sp   = norm(pamp)
        sq   = norm(qamp)

        # return invariants
        return nvec, namp, sp, sq


    # Failure criteria (or yield surface)
    # ===================================
    def func (self, sig, typ):

        # yield surface
        ysurf = False
        if typ[0]=='y':
            ysurf = True
            typ   = typ[1:]
            pc    = self.sc/2.0

        # sin phi
        sphi = sin(self.phi*pi/180.0)

        # von Mises
        if typ=='VM':
            p, q = sig_calc_p_q(sig)
            if self.fc_psa: k = sqrt(2.0)*self.fc_cu
            else:           k = 2.0*(sqrt(2.0)/sqrt(3.0))*self.fc_cu
            if ysurf: raise Exception('func: ysurf is not available with VM')
            else: f = q - k

        # Drucker/Prager
        elif typ=='DP':
            cbar = sqrt(3.0)*self.c/tan(self.phi*pi/180.0)
            kdp  = 2.0*sqrt(2.0)*sphi/(3.0-sphi)
            p, q = sig_calc_p_q(sig)
            if ysurf: f = ((p-pc)/pc)**2.0 + (q/(kdp*pc))**2.0 - 1.0
            else:     f = q - (p + cbar)*kdp

        # Mohr/Coulomb
        elif typ=='MC':
            p, q = sig_calc_p_q (sig)
            t    = sig_calc_t   (sig)
            th   = arcsin(t)/3.0
            cbar = sqrt(3.0)*self.c/tan(self.phi*pi/180.0)
            g    = sqrt(2.0)*sphi/(sqrt(3.0)*cos(th)-sphi*sin(th))
            if ysurf: raise Exception('func: ysurf is not available with MC')
            else: f = q - (p + cbar)*g

        # Nakai/Matsuoka
        elif typ=='NM':
            #cbar     = sqrt(3.0)*self.fc_c/tan(self.fc_phi*pi/180.0)
            #kmn      = (9.0-sphi**2.0)/(1.0-sphi**2.0)
            #sig0     = cbar*matrix([[1.0],[1.0],[1.0],[0.0]])
            #sig_     = sig+sig0
            tgphi2   = tan(self.phi*pi/180.0)**2.0
            kmn      = 9.0 + 8.0*tgphi2
            sig_     = sig
            l        = sig_calc_s123(sig_)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig_)
            if ysurf: raise Exception('func: ysurf is not available with NM')
            else: f = I1*I2 - kmn*I3

        # Nakai/Matsuoka nonlinear
        elif typ=='NMnl':
            p, q = sig_calc_p_q (sig, 'cam')
            M    = self.refcurve(p)/p if p>0.0 else self.fc_prms['A']
            sphi = 3.0*M/(M+6.0)
            kmn  = (9.0-sphi**2.0)/(1.0-sphi**2.0)
            l    = sig_calc_s123(sig)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig)
            if ysurf: raise Exception('func: ysurf is not available with NMnl')
            else: f = I1*I2 - kmn*I3

        # Lade/Duncan
        elif typ=='LD':
            cbar     = sqrt(3.0)*self.fc_c/tan(self.fc_phi*pi/180.0)
            kld      = ((3.0-sphi)**3.0)/((1.0+sphi)*((1.0-sphi)**2))
            sig0     = cbar*matrix([[1.0],[1.0],[1.0],[0.0]])
            sig_     = sig+sig0
            l        = sig_calc_s123(sig_)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig_)
            if ysurf: raise Exception('func: ysurf is not available with LD')
            else: f = I1**3.0 - kld*I3

        # Argyris/Sheng
        elif typ=='AS':
            p, q, t = sig_calc_pqt (sig, 'cam')
            Mcs     = 6.0*sphi/(3.0-sphi)
            om      = ((3.0-sphi)/(3.0+sphi))**4.0
            M       = Mcs*(2.0*om/(1.0+om-(1.0-om)*t))**0.25;
            if ysurf: raise Exception('func: ysurf is not available with AS')
            else: f = q/p - M

        # Anisotropic
        elif typ=='AMP' or typ=='AMPb' or typ=='AMPba':
            l, Q = sig_calc_rot (sig)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return 1.0e+8
            nvec, namp, sp, sq = self.get_nvec_namp_sig_tau (l, Q)
            if ysurf:
                z0 = 2.0*pc/sqrt(3.0)
                f  = log(sp/z0) + (1.0/self.bet)*(sq/(self.R*sp))**self.bet
                #f  = ((sp-pc)/pc)**2.0 + (sq/(self.R*pc))**2.0 - 1.0
            else: f = sq - self.R * sp

        # error
        else: raise Exception('func: typ==%s is invalid' % typ)

        # return function value
        return f


    # Failure criteria names
    # ======================
    def names (self, typ):
        if   typ=='VM':    return 'von Mises'
        elif typ=='DP':    return 'Drucker/Prager'
        elif typ=='MC':    return 'Mohr/Coulomb'
        elif typ=='NM':    return 'Nakai/Matsuoka'
        elif typ=='NMnl':  return 'Nakai/Matsuoka (non-linear)'
        elif typ=='LD':    return 'Lade/Duncan'
        elif typ=='AS':    return 'Argyris/Sheng'
        elif typ=='AMP':   return 'Anisotropic'
        elif typ=='AMPb':  return r'$b=%g$'%self.b
        elif typ=='AMPba': return r'$b=%g,\,\alpha=%g$'%(self.b,self.alp)
        else: raise Exception('failure_crit_names: typ==%s is invalid' % typ)


    # Rosette
    # =======
    # th   : show theta angles
    # ref  : reference lines
    # pos  : positive values
    # fsz  : font size
    def rst (self, th=False, ref=False, pos=False, fsz=10):

        # radius
        r = 1.*phi_calc_M(self.phi,'oct') if self.r==None else self.r

        # constants
        cf = 0.2
        cr = 1.1
        l1 = (             0.0  , cr*r            ) # line: 1 end points
        l2 = (-cr*r*cos(pi/6.0) ,-cr*r*sin(pi/6.0)) # line: 2 end points
        l3 = ( cr*r*cos(pi/6.0) ,-cr*r*sin(pi/6.0)) # line: 3 end points
        l4 = (-cr*r*cos(pi/6.0) , cr*r*sin(pi/6.0)) # line: 4 = neg 1 end points
        lo = (-cr*r*cos(pi/3.0) , cr*r*sin(pi/3.0)) # line: origin of cylindrical system

        # main lines
        plot ([0.0,l1[0]],[0.0,l1[1]],'k-', color='grey', zorder=0)
        plot ([0.0,l2[0]],[0.0,l2[1]],'k-', color='grey', zorder=0)
        plot ([0.0,l3[0]],[0.0,l3[1]],'k-', color='grey', zorder=0)

        # reference
        plot ([0.0, l4[0]],[0.0, l4[1]],'--', color='grey', zorder=-1)
        plot ([0.0,-l4[0]],[0.0, l4[1]],'--', color='grey', zorder=-1)
        plot ([0.0,   0.0],[0.0,-l1[1]],'--', color='grey', zorder=-1)
        if ref:
            plot ([0.0, lo[0]],[0.0, lo[1]],'--', color='grey', zorder=-1)
            plot ([0.0,-lo[0]],[0.0, lo[1]],'--', color='grey', zorder=-1)
            plot ([-cr*r,cr*r],[0.0,0.0],   '--', color='grey', zorder=-1)

        # text
        if self.sxyz: k1,k2,k3 = 'z','y','x'
        else:         k1,k2,k3 = '1','3','2'
        if pos:
            if th: t1 = text(l1[0],l1[1],r'$\sigma_%s,\theta=+30^\circ$'%k1, ha='center', fontsize=fsz)
            else:  t1 = text(l1[0],l1[1],r'$\sigma_%s$'%k1,                  ha='center', fontsize=fsz)
            t2 = text(l2[0],l2[1],r'$\sigma_%s$'%k2,  ha='right',  fontsize=fsz)
            t3 = text(l3[0],l3[1],r'$\sigma_%s$'%k3,  ha='left',   fontsize=fsz)
        else:
            if th: t1 = text(l1[0],l1[1],r'$-\sigma_%s,\theta=+30^\circ$'%k1, ha='center', fontsize=fsz)
            else:  t1 = text(l1[0],l1[1],r'$-\sigma_%s$'%k1,                  ha='center', fontsize=fsz)
            t2 = text(l2[0],l2[1],r'$-\sigma_%s$'%k2,  ha='right',  fontsize=fsz)
            t3 = text(l3[0],l3[1],r'$-\sigma_%s$'%k3,  ha='left',   fontsize=fsz)
        self.rst_txt.append (t1)
        self.rst_txt.append (t2)
        self.rst_txt.append (t3)
        if th:
            t4 = text(lo[0],lo[1],r'$\theta=0^\circ$',   ha='center', fontsize=fsz)
            t5 = text(l4[0],l4[1],r'$\theta=-30^\circ$', ha='center', fontsize=fsz)
            self.rst_txt.append (t4)
            self.rst_txt.append (t5)


    # Plot failure criteria
    # =====================
    def plot (self, typ='NM', clr='red', lst='-', lwd=1, label=None, np=40, leg_set=True, show_phi=True, fsz=10):

        # radius
        r = 1.*phi_calc_M(self.phi,'oct') if self.r==None else self.r

        # contour
        f     = zeros ((np,np))
        sa    = zeros ((np,np))
        sb    = zeros ((np,np))
        samin = -1.1*r if self.samin==None else self.samin
        samax =  1.1*r if self.samax==None else self.samax
        sbmin = -1.1*r if self.sbmin==None else self.sbmin
        sbmax =  1.1*r if self.sbmax==None else self.sbmax
        dsa   = (samax-samin)/np
        dsb   = (sbmax-sbmin)/np
        for i in range(np):
            for j in range(np):
                sa[i,j] = samin + i*dsa
                sb[i,j] = sbmin + j*dsb
                if self.sxyz: s = oct_calc_sxyz (sa[i,j], sb[i,j], self.sc)
                else:         s = oct_calc_s123 (sa[i,j], sb[i,j], self.sc)
                sig    = matrix([[s[0]],[s[1]],[s[2]],[0.0]])
                f[i,j] = self.func (sig, typ)
        contour (sa,sb,f, [0.0], colors=clr, linestyles=lst, linewidths=lwd)

        # set axis
        gca().set_xticks([])
        gca().set_yticks([])
        gca().set_frame_on (False)
        axis ('equal')

        # legend
        if leg_set:
            self.leg_plt.append (plot ([0],[0], color=clr, linestyle=lst, linewidth=lwd))
            self.leg_txt.append (self.names (typ))

        # show phi
        if show_phi:
            fmt  = '%g'
            sphi = sin(self.phi*pi/180.0)
            q    = self.sc * 2.0*sqrt(2.0)*sphi/(3.0+sphi)
            s    = '$\phi=' + fmt + '^\circ$'
            text (0.0, -1.05*q, s%self.phi, fontsize=fsz, va='top')


    # Legend
    # ======
    def leg (self, fsz=8, ncol=2):
        legend (self.leg_plt, self.leg_txt, bbox_to_anchor=(0,0,1,1), loc=3, ncol=ncol, mode='expand',
                borderaxespad=0., handlelength=3, prop={'size':fsz})


    # Plot 3D failure criteria
    # ========================
    def plot3d (self, typs=['MC','NM'], clr='red', np=40):

        # radius
        r = 1.*phi_calc_M(self.phi,'oct') if self.r==None else self.r

        # contour
        F = []
        for typ in typs: F.append (zeros((np,np,np)))
        sa    = zeros ((np,np,np))
        sb    = zeros ((np,np,np))
        sc    = zeros ((np,np,np))
        samin = -1.1*r       if self.samin==None else self.samin
        samax =  1.1*r       if self.samax==None else self.samax
        sbmin = -1.1*r       if self.sbmin==None else self.sbmin
        sbmax =  1.1*r       if self.sbmax==None else self.sbmax
        scmin =  0.0         if self.scmin==None else self.scmin
        scmax =  1.1*self.sc if self.scmax==None else self.scmax
        dsa   = (samax-samin)/np
        dsb   = (sbmax-sbmin)/np
        dsc   = (scmax-scmin)/np
        for i in range(np):
            for j in range(np):
                for k in range(np):
                    sa[i,j,k] = samin + i*dsa
                    sb[i,j,k] = sbmin + j*dsb
                    sc[i,j,k] = scmin + k*dsc
                    if self.sxyz: s = oct_calc_sxyz (sa[i,j,k], sb[i,j,k], sc[i,j,k])
                    else:         s = oct_calc_s123 (sa[i,j,k], sb[i,j,k], sc[i,j,k])
                    sig      = matrix([[s[0]],[s[1]],[s[2]],[0.0]])
                    for m, typ in enumerate(typs): F[m][i,j,k] = self.func (sig, typ)

        from enthought.mayavi.mlab import contour3d
        from enthought.mayavi.mlab import show as mlab_show

        for m, typ in enumerate(typs): contour3d (sa, sb, sc, F[m], contours=[0.0])
        mlab_show ()
