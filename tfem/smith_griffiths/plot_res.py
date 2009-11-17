from numpy import *
from pylab import *
from msys_readdata import *

if False:
    dat = read_table("fig_11_01_nod_1_-200.res")

    tsw = 1.0
    def calc_U(t):
        if (t<tsw): return  0.441*sin(pi*t/tsw)-0.216*sin(2.0*pi*t/tsw);
        else:       return -0.432*sin(6.284*(t-tsw))

    T = linspace(0.0,1.8,100)
    U = zeros(len(T))
    for i in range(len(T)): U[i] = calc_U(T[i])

    plot(T,U,'b-')
    plot(dat['Time'],dat['uy'],'ro',lw=2)
    grid()
    show()

if False:
    res   = read_table("fig_11_04_nod_17_-100.res")
    p113d = read_table("sg_11_07_p113.dat")
    p113  = read_table("sg_11_07_p113.res")
    p114  = read_table("sg_11_07_p114.res")

    plot(res  ['Time'],res  ['uy'],'r-',lw=2)
    plot(p113d['Time'],p113d['uy'],'bo')
    plot(p113 ['Time'],p113 ['uy'],'b-')
    plot(p113 ['Time'],p113 ['uy'],'b+')
    #plot(p114 ['Time'],p114 ['uy'],'b-',lw=1)
    #plot(p114 ['Time'],p114 ['uy'],'b+',lw=1)

    #legend(['MechSys','SG:p113d','SG:p113','SG:p114'],loc='upper left')
    legend(['MechSys','SG:p113d','SG:p113'],loc='upper left')
    grid()
    show()

if False:
    res  = read_table("fig_11_19_nod_30_-200.res")
    p117 = read_table("sg_11_19.res")

    plot(res ['Time'],res ['uy'],'r-',lw=2)
    plot(p117['Time'],p117['uy'],'b-')
    plot(p117['Time'],p117['uy'],'b+')

    #legend(['MechSys','SG:p113d','SG:p113','SG:p114'],loc='upper left')
    #legend(['MechSys','SG:p113d','SG:p113'],loc='upper left')
    grid()
    show()

if True:
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
