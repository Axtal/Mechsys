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

# CamClay: closed-form solution: ev, ed
# =====================================
# Proportional loading: k = dq/dp (path)
def CamClay_ev_ed(prms, inis, sig0, sig):

    # parameters and initial values
    lam = prms['lam']
    kap = prms['kap']
    nu  = prms['nu']
    v0  = inis['v0']
    M   = phi_calc_M (prms['phi'], 'oct')
    chi = (kap-lam)/v0

    # auxiliary variables
    alp    = (3.0*(1.0-2.0*nu))/(2.0*(nu+1.0))
    p0, q0 = sig_calc_p_q (sig0, 'oct')
    p,  q  = sig_calc_p_q (sig,  'oct')
    dp, dq = p-p0, q-q0
    k      = dq/dp
    r      = q/(M*p)
    r0     = q0/(M*p0)

    # elastic strains
    eve = -(kap*log(p/p0))/v0
    ede = sqrt(3.0)*k*kap*log(p/p0)/(2.0*alp*v0)
    if dp<0 or dq<0: return eve, ede

    # elastoplastic strains
    evp = chi*(log((r**2.0+1.0)/(r0**2.0+1.0))+log(p/p0))
    edp = 2.0*chi*( k*M*log(p/p0)/(k**2.0-M**2.0) - k*log((r+1.0)/(r0+1.0))/(2.0*k+2.0*M) + k*log((r-1.0)/(r0-1.0))/(2.0*k-2.0*M) + arctan(r) - arctan(r0) )/(sqrt(3.0)*M)
    ev  = eve + evp
    ed  = ede + edp
    return ev, ed
