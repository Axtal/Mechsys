from msys_fig import *

r = read_table("test_timesteps.res")
plot   (r['i'], r['dt_sch0'], 'b-', marker='o', label='sch:0')
plot   (r['i'], r['dt_sch1'], 'r-', marker='.', label='sch:1')
xlabel ('i')
ylabel ('dt')
legend (loc='best')
Grid   ()
show   ()
