from numpy import *
from pylab import *
r = 2.5
A = linspace(0.,pi/2.,200)
X = r*cos(A)
Y = r*sin(A)
plot (X,Y,'r-',lw=2)
r = 3.5
X = r*cos(A)
Y = r*sin(A)
plot (X,Y,'r-',lw=2)

