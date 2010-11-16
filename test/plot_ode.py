from numpy import exp, log, linspace
from pylab import plot, show, xlabel, ylabel, grid, legend, axis, subplot
from msys_readdata import *

fsz = {'fontsize':20}
dat = read_table("ode1.dat")
res = read_table("test_ode1.dat")

subplot(2,1,1)
plot   (dat['t'], dat['y0'],   'b-', label='GSL')
plot   (dat['t'], dat['y0me'], 'r.', label='ME')
plot   (res['t'], res['y0'],   'k-', label='test_ode1')
xlabel (r'$t$',  fsz)
ylabel (r'$y_0$',fsz)
legend ()
grid   ()

subplot(2,1,2)
plot   (dat['y0'],   dat['y1'],   'b-', label='GSL')
plot   (dat['y0me'], dat['y1me'], 'r.', label='ME')
plot   (res['y0'],   res['y1'],   'k-', label='test_ode1')
xlabel (r'$y_0$',  fsz)
ylabel (r'$y_1$',fsz)
legend ()
#axis   ('scaled')
grid   ()

show   ()
