from msys_fig import *

res = read_table("test_wrcs.res")
subplot(2,1,1)
plot   (res['pc'], res['Sw'])
xlabel ('pc')
ylabel ('Sw')
Grid   ()
subplot(2,1,2)
plot   (res['pc'], res['kw'])
xlabel ('pc')
ylabel ('kw')
Grid   ()
show   ()
