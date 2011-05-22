from msys_fig import *
from pylab import loglog

res = read_table("test_wrcs.res")

subplot(2,2,1)
plot   (res['pc'], res['Sw'], lw=2, clip_on=False)
axis   ([0,max(res['pc']), 0, 1])
xlabel ('pc [kPa]')
ylabel ('Sw')
Grid   ()

subplot(2,2,2)
plot   (log(1.0+res['pc']), res['Sw'], lw=2, clip_on=False)
axis   ([0,max(log(1.0+res['pc'])), 0, 1])
xlabel ('log(1+pc) [kPa]')
ylabel ('Sw')
Grid   ()

subplot(2,2,3)
plot   (res['pc'], res['kw'], lw=2, clip_on=False)
xlabel ('pc [kPa]')
ylabel ('kwr')
Grid   ()

subplot(2,2,4)
loglog (res['pc'], res['kw'], lw=2, clip_on=False)
xlabel ('pc [kPa]')
ylabel ('kwr')
Grid   ()
show   ()
