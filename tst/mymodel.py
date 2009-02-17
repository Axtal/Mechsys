# mymodel.py

import math

def stiff(prms, deps, sig,eps,ivs):

    SQ2 = math.sqrt(2.0)
    E   = prms['a0'] # Young
    nu  = prms['a1'] # Poisson

    c  = E/((1.0+nu)*(1.0-2.0*nu))
    c1 = c*(1.0-nu)
    c2 = c*(1.0-2.0*nu)/2.0
    c3 = c*nu;

    De = (c1     , c3     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
          c3     , c1     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
          c3     , c3     , c1     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
          0.0*SQ2, 0.0*SQ2, 0.0*SQ2, c2 *2.0, 0.0*2.0, 0.0*2.0,
          0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, c2 *2.0, 0.0*2.0,
          0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, 0.0*2.0, c2 *2.0) # In Mandel's basis

    return De, []
