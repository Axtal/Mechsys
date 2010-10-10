from numpy import *
from pylab import *
from msys_readdata import *

plt = 1
ref = 0

quad8 = 1
both  = False
ninc  = 10
STOL  = 1.0e-7
NR    = 1
FE    = 0

if plt==0: fkey = "abbo_sloan_01"
else:      fkey = "abbo_sloan_02"
if not both:
    if quad8: fkey += "_quad8"
    else:     fkey += "_tri15"
if NR or FE: fkey += "_%d"%ninc
else:        fkey += "_%d_%g"%(ninc,STOL)

def plot_one(data,txt,clr,marker,lwd,hline=False,dutxt=0.0,a=1.0,c=10.0,t=0.2):
    fx = []
    nn = len(data)
    for i in range(nn): fx.append (array(data[i]['fx']))
    np = len(fx[0])
    q  = zeros(np)
    u  = array(data[0]['ux'])
    for i in range(np):
        f = zeros(nn)
        for j in range(nn): f[j] = fx[j][i]
        norm_f = sqrt(sum(f*f))
        if   nn==5: q[i] = 90.0*norm_f/(sqrt(2290.0)*t)
        elif nn==3: q[i] = sqrt(2.0)*norm_f/t
        else: raise Exception('nn=%d is wrong'%nn)
    plot (u*1000.0/a, q/c, color=clr, marker=marker, lw=lwd, label=txt)
    if hline:
        max_q = max(q)
        max_u = max(u)
        axhline (max_q/c, lw=lwd, color=clr)
        text    (max_u+dutxt,max_q/c,r'$%g$'%(max_q/c), color=clr)

if plt==0:
    # loaded nodes
    if quad8:
        quad8 = [read_table("%s_nod_0.res" %fkey),
                 read_table("%s_nod_6.res" %fkey),
                 read_table("%s_nod_22.res"%fkey)]
    else:
        tri15 = [read_table("%s_nod_0.res" %fkey),
                 read_table("%s_nod_6.res" %fkey),
                 read_table("%s_nod_22.res"%fkey),
                 read_table("%s_nod_44.res"%fkey),
                 read_table("%s_nod_45.res"%fkey)]

    # with reference solution
    if ref: passs

    if quad8: plot_one (quad8,'Quad8','red', 'None',2,True,0.05)
    else:     plot_one (tri15,'Tri15','blue','None',1,True,0.1)

    axhline (1.0174, color='k')
    text    (0,1.0174,r'$1.0174$')
    xlabel  (r'$\frac{u}{a}\,\times\,10^3$')
    ylabel  (r'$\frac{q}{c}$')
    legend  (loc='lower right')
    grid()
    show()

def plot_one2(L,B,data,txt,clr,marker,lwd,hline=False,dutxt=0.0,c=10.0):
    fy = []
    nn = len(data)
    for i in range(nn): fy.append (array(data[i]['fy']))
    np = len(fy[0])
    q  = zeros(np)
    u  = -array(data[0]['uy'])
    for i in range(np):
        f = zeros(nn)
        for j in range(nn): f[j] = fy[j][i]
        norm_f = sqrt(sum(f*f))
        if   nn==5: q[i] = 90.0*norm_f/(sqrt(2290.0)*L)
        elif nn==3: q[i] = sqrt(2.0)*norm_f/L
        else: raise Exception('nn=%d is wrong'%nn)
    plot (u/B, q/c, color=clr, marker=marker, lw=lwd, label=txt)
    if hline:
        max_q = max(q)
        max_u = max(u)
        axhline (max_q/c, lw=lwd, color=clr)
        text    (dutxt,max_q/c,r'$%g$'%(max_q/c), color=clr)

if plt==1:
    if quad:
        quad8 = [read_table("abbo_sloan_02_quad8_nod_12.res"),
                 read_table("abbo_sloan_02_quad8_nod_13.res"),
                 read_table("abbo_sloan_02_quad8_nod_43.res")]
    else:
        tri15 = [read_table("abbo_sloan_02_tri15_nod_0.res"  ),
                 read_table("abbo_sloan_02_tri15_nod_207.res"),
                 read_table("abbo_sloan_02_tri15_nod_60.res" ),
                 read_table("abbo_sloan_02_tri15_nod_206.res"),
                 read_table("abbo_sloan_02_tri15_nod_13.res" )]

    if ref:
        tri15_ref = [read_table("abbo_sloan_01_tri15_ref_nod_0.res"),
                     read_table("abbo_sloan_01_tri15_ref_nod_6.res"),
                     read_table("abbo_sloan_01_tri15_ref_nod_22.res"),
                     read_table("abbo_sloan_01_tri15_ref_nod_44.res"),
                     read_table("abbo_sloan_01_tri15_ref_nod_45.res")]

        quad8_ref = [read_table("abbo_sloan_01_quad8_ref_nod_0.res"),
                     read_table("abbo_sloan_01_quad8_ref_nod_6.res"),
                     read_table("abbo_sloan_01_quad8_ref_nod_22.res")]

    L = 0.564
    B = 2.0
    if quad: plot_one2 (L,B,quad8,'Quad8','red', 'None',2,True,0.05)
    else:    plot_one2 (L,B,tri15,'Tri15','blue','None',1,True,dutxt=0.02)

    axhline (30.1396, color='k')
    text    (0,30.1396,r'$30.1396$')
    xlabel  (r'-$\frac{u}{B}$')
    ylabel  (r'$\frac{q}{c}$')
    legend  (loc='lower right')
    grid()
    show()
