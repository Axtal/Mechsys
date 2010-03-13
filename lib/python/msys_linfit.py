from numpy import array, sqrt

# Linear fitting
# ==============
# ls:  least-square
# tls: total least-square
# cmx: y = c + m*x
class LinFit:
    def __init__(self, X, Y, tls=True, cmx=True):
        # data
        self.cmx = cmx
        self.c   = 0.
        sX, sY, sXY, sXX = sum(X), sum(Y), sum(X*Y), sum(X*X)
        n = len(X)

        # do a least-square (LS) first
        if self.cmx:
            den    = sX**2. - n*sXX
            self.c = (sX*sXY - sXX*sY)/den
            self.m = (sX*sY  -  n*sXY)/den
        else:
            self.m = sXY/sXX

        # check if y increases with x
        increasing = True if self.m>=0.0 else False

        # total least-square (TLS)
        if tls:
            sYY = sum(Y*Y)
            if self.cmx:
                xb = sX/n
                yb = sY/n
                A  = n*xb*yb - sXY
                B  = sYY - sXX - n*yb**2. + n*xb**2.
                C  = sXY - n*xb*yb
            else:
                A = -sXY
                B = sYY - sXX
                C = sXY
            d = B**2. - 4.*A*C
            if d<0.0: raise Exception('LinFit::Constructor: (TLS) delta == %g is negative'%d)
            m1 = (-B + sqrt(d))/(2.*A)
            m2 = (-B - sqrt(d))/(2.*A)
            if increasing: self.m = m1 if m1>0.0 else m2
            else:          self.m = m1 if m1<0.0 else m2
            self.c = yb - self.m*xb if self.cmx else 0.0

    # Calculate y(x)
    # ==============
    def y(self, x): return self.c + self.m*x

    # Sum residual
    # ============
    def S(self, X, Y, tls=False):
        res = 0.0
        cf  = (1.+self.m**2.) if tls else 1.
        for i, x in enumerate(X):
            res += ((Y[i] - self.c - self.m*X[i])**2.)/cf
        return res


# Test
# ====
if __name__=='__main__':
    from pylab import *
    X  = array ([1., 2., 3., 4.])
    Y  = array ([6., 5., 7., 10.])
    f0 = LinFit (X,Y, tls=True,  cmx=True )
    f1 = LinFit (X,Y, tls=True,  cmx=False)
    f2 = LinFit (X,Y, tls=False, cmx=True )
    f3 = LinFit (X,Y, tls=False, cmx=False)
    x  = linspace(0,max(X),100)
    plot (X,Y,'ro',clip_on=False)
    plot (x,f0.y(x),'c-', marker='.', markevery=10, label='TLS: c=%2.2f, m=%2.2f, S=%2.2f, St=%2.2f'%(f0.c, f0.m, f0.S(X,Y), f0.S(X,Y,True)))
    plot (x,f1.y(x),'r-', marker='.', markevery=10, label='TLS: c=%2.2f, m=%2.2f, S=%2.2f, St=%2.2f'%(f1.c, f1.m, f1.S(X,Y), f1.S(X,Y,True)))
    plot (x,f2.y(x),'c-', marker='+', markevery=10, label='LS:  c=%2.2f, m=%2.2f, S=%2.2f, St=%2.2f'%(f2.c, f2.m, f2.S(X,Y), f2.S(X,Y,True)))
    plot (x,f3.y(x),'r-', marker='+', markevery=10, label='LS:  c=%2.2f, m=%2.2f, S=%2.2f, St=%2.2f'%(f3.c, f3.m, f3.S(X,Y), f3.S(X,Y,True)))
    axis ('equal')
    grid ()
    legend(loc='best')
    show ()
