import optparse
from msys_fig import *

# input
op = optparse.OptionParser()
op.add_option('--tst', '-t', dest='tst', default='0', help='test number')
opts, args = op.parse_args()

if opts.tst=='0':
    dat = read_table("fig_11_01_nod_1.res")

    tsw = 1.0
    def calc_U(t):
        if (t<tsw): return  0.441*sin(pi*t/tsw)-0.216*sin(2.0*pi*t/tsw);
        else:       return -0.432*sin(6.284*(t-tsw))

    T = linspace(0.0,1.8,100)
    U = zeros(len(T))
    for i in range(len(T)): U[i] = calc_U(T[i])

    plot(T,U,'b-',label='Analytical')
    plot(dat['Time'],dat['uy'],'ro',lw=2,label='FEM')
    xlabel('Time')
    ylabel('uy')
    legend(loc='best')
    grid()
    show()

if opts.tst=='1':
    res   = read_table("fig_11_04_nod_17.res")
    p113d = read_table("sg_11_07_p113.dat")
    p113  = read_table("sg_11_07_p113.sim")
    p114  = read_table("sg_11_07_p114.sim")
    subplot(2,1,1)
    plot   (res  ['Time'],res  ['uy'],'r-',lw=2,      label='MechSys')
    plot   (p113d['Time'],p113d['uy'],'bo',           label='SG:prog11.3(scan)')
    plot   (p113 ['Time'],p113 ['uy'],'b-',marker='+',label='SG:prog11.3(sim)')
    plot   (p114 ['Time'],p114 ['uy'],'g-',lw=1,      label='SG:prog11.4')
    legend (loc='best')
    xlabel ('Time')
    ylabel ('uy')
    Grid   ()
    subplot(2,1,2)
    plot   (res['Time'],res['sy'],'r-',lw=2, label='MechSys')
    xlabel ('Time')
    ylabel ('sy')
    Grid   ()
    show   ()

if opts.tst=='2':
    res  = read_table("fig_11_19_nod_30.res")
    p117 = read_table("sg_11_19.sim")

    plot(res ['Time'],res ['uy'],'r-',lw=2)
    plot(p117['Time'],p117['uy'],'b-')
    plot(p117['Time'],p117['uy'],'b+')

    #legend(['MechSys','SG:p113d','SG:p113','SG:p114'],loc='upper left')
    #legend(['MechSys','SG:p113d','SG:p113'],loc='upper left')
    grid()
    show()

if opts.tst=='3':
    nod = 17
    res = read_table("fig_11_04.out")

    # displacements
    subplot(2,2,1)
    plot(res["Time"],res["ux"],'r-',lw=2)
    plot(res["Time"],res["uy"],'b-',lw=2)
    xlabel("Time");  ylabel("ux, uy")
    legend(["ux","uy"],loc="upper left")
    title("Displacements at node %d"%nod)
    grid()

    # velocities
    subplot(2,2,2)
    plot(res["Time"],res["vx"],'r-',lw=2)
    plot(res["Time"],res["vy"],'b-',lw=2)
    xlabel("Time");  ylabel("vx, vy")
    legend(["vx","vy"],loc="upper left")
    title("Velocities at node %d"%nod)
    grid()

    # accelerations
    subplot(2,2,3)
    plot(res["Time"],res["ax"],'r-',lw=2)
    plot(res["Time"],res["ay"],'b-',lw=2)
    xlabel("Time");  ylabel("ax, ay")
    legend(["ax","ay"],loc="upper left")
    title("Accelerations at node %d"%nod)
    grid()

    # disp-veloc
    subplot(2,2,4)
    plot(res["ux"],res["vx"],'r-',lw=2)
    plot(res["uy"],res["vy"],'b-',lw=2)
    xlabel("ux, uy");  ylabel("vx, vy")
    legend(["vx-ux","vy-uy"],loc="upper left")
    title("Disp-Veloc at node %d"%nod)
    grid()

    show()
