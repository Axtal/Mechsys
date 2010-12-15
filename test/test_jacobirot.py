from mechsys  import *
from msys_fig import *
from scipy    import diag
import numpy.linalg as npyla

Q,L = [], []
A   = [[1., 2., 3.],
       [2., 3., 2.],
       [3., 2., 2.]]
it  = JacobiRot(A,Q,L)
A   = matrix(A)
Q   = matrix(Q)
L   = diag(L)
l,q = npyla.eig(A)
l   = diag(l)
print 'it = ', it
print 'A  =\n', A
print 'Q  =\n', Q
print 'q  =\n', q
print 'L  =\n', L
print 'l  =\n', l
print 'A = Q* L *Q.T = \n', Q*L*Q.T
print 'A = q* l *q.T = \n', q*l*q.T
print 'L = Q.T* A *Q = \n', Q.T*A*Q
print 'L = q.T* A *q = \n', q.T*A*q
