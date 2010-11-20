from msys_fig import *
from pylab import loglog

res = read_table("test_wrcs.res")

subplot(3,1,1)
plot   (res['pc'], res['Sw'], lw=2)
axis   ([0,max(res['pc']), 0, 1])
xlabel ('pc [kPa]')
ylabel ('Sw')
Grid   ()

subplot(3,1,2)
plot   (res['pc'], res['kw'], lw=2)
xlabel ('pc [kPa]')
ylabel ('kwr')
Grid   ()

subplot(3,1,3)
loglog (res['pc'], res['kw'], lw=2)
xlabel ('pc [kPa]')
ylabel ('kwr')
Grid   ()
show   ()
