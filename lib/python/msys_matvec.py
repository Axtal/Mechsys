# CIVL4250 Numerical Methods in Geomechanics
# Copyright (C) 2009 Dorival M Pedroso
# ------------------------------------------
# Matrix and Vector classes

from numpy import matrix, linalg, sqrt, zeros

class Vector(matrix):
    # Constructor
    # ===========
    def __init__(self, vals, dtype=float, copy=False):
        super(Vector,self).__init__()
        if self.shape[1]>self.shape[0]:
            self.resize((self.shape[1],self.shape[0]))

    # Access item
    # ===========
    def __getitem__(self, key):
        return matrix.__getitem__(self, (key,0))

    # Euclidian norm
    # ==============
    def norm(self):
        return sqrt((self.T*self)[0])

    # Nice Print
    # ==========
    def write(self, nf='%10g', Tol=1.0e-13):
        m = self.shape[0] # number of rows
        lin = ''          # empty string
        for i in range(m):
            if abs(self[i])<Tol: lin += nf % 0
            else:                lin += nf % self[i]
            lin += '\n'
        print lin

class Matrix(matrix):
    # Determinant
    # ===========
    def det(self): return linalg.det(self)

    # Nice Print
    # ==========
    def write(self, nf='%10g', Tol=1.0e-13):
        m = self.shape[0] # number of rows
        n = self.shape[1] # number of columns
        lin = ''          # empty string
        for i in range(m):
            for j in range(n):
                if abs(self[i,j])<Tol: lin += nf % 0
                else:                  lin += nf % self[i,j]
            lin += '\n'
        print lin

def Dot (U,V): return (U.T*V)[0]
def Dyad(U,V):
    m = U.shape[0]
    n = V.shape[0]
    M = Matrix(zeros(shape=(m,n)))
    for i in range(m):
        for j in range(n):
            M[i,j] = U[i]*V[j]
    return M

if __name__=="__main__":

    K = Matrix([[1., 2., 3., 4.],
                [5., 6., 7., 8.]])

    v = Vector([1., 2., 3., 4.])

    vdyv = Dyad(v,v)

    print 'K ='; K.write()
    print 'v ='; v.write()
    print 'v dyad v ='; vdyv.write()
