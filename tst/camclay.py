# mymodel.py

import math

def calc_p_q(sig):
    p = (sig[0]+sig[1]+sig[2])/3.0
    q = math.sqrt(((sig[0]-sig[1])**2.0 + (sig[1]-sig[2])**2.0 + (sig[2]-sig[0])**2.0 + 3.0*(sig[3]**2.0 + sig[4]**2.0 + sig[5]**2.0) )/2.0)
    return p, q

def calc_M(prms):
    phi    = prms['a2']
    sinphi = math.sin(phi*math.pi/180.0)
    M      = 6.0*sinphi/(3.0-sinphi)
    return M

def calc_De(E,nu):
    SQ2 = math.sqrt(2.0)
    c   = E/((1.0+nu)*(1.0-2.0*nu))
    c1  = c*(1.0-nu)
    c2  = c*(1.0-2.0*nu)/2.0
    c3  = c*nu;
    De  = (c1     , c3     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
           c3     , c1     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
           c3     , c3     , c1     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
           0.0*SQ2, 0.0*SQ2, 0.0*SQ2, c2 *2.0, 0.0*2.0, 0.0*2.0,
           0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, c2 *2.0, 0.0*2.0,
           0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, 0.0*2.0, c2 *2.0) # In Mandel's basis
    return De

def init_ivs(prms, ini,sig,eps):
    p, q = calc_p_q (sig)
    M    = calc_M   (prms)
    z0   = p+q*q/(M*M*p)
    z1   = ini['v']
    return [z0,z1]

def stiff(prms, deps, sig,eps,ivs):
    # parameters
    lam = prms['a0'] # Lambda
    kap = prms['a1'] # Kappa
    phi = prms['a2'] # Phi
    nu  = prms['a3'] # Poisson

    # invariants
    p, q = calc_p_q(sig)

    # Young
    E = -3.0*p*(1.0-2.0*nu)*ivs[1]/kap

    B0 = (0,0,0,0,0,0)
    B1 = (0,0,0,0,0,0)

    return calc_De(E,nu), [B0,B1]
