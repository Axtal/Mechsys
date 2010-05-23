from numpy import *
from pylab import *
from msys_readdata import *

plt = 0

def plot_one(tri15,quad8,txt,tclr,qclr):
    tri15_fx  = [array(tri15[0]['fx']),
                 array(tri15[1]['fx']),
                 array(tri15[2]['fx']),
                 array(tri15[3]['fx']),
                 array(tri15[4]['fx'])]
    quad8_fx  = [array(quad8[0]['fx']),
                 array(quad8[1]['fx']),
                 array(quad8[2]['fx'])]
    a = 1.0
    c = 10.0
    L = 0.2
    h = 1.0
    n = len(tri15_fx[0])
    tri15_q = zeros(n)
    quad8_q = zeros(n)
    for i in range(n):
        tri15_f = zeros(5)
        quad8_f = zeros(3)
        for j in range(5): tri15_f[j] = tri15_fx[j][i]
        for j in range(3): quad8_f[j] = quad8_fx[j][i]
        tri15_norm_f = sqrt(sum(tri15_f*tri15_f))
        quad8_norm_f = sqrt(sum(quad8_f*quad8_f))
        tri15_q[i]   = 90.0*tri15_norm_f/(sqrt(2290.0)*L*h)
        quad8_q[i]   = sqrt(2.0)*quad8_norm_f/(L*h)

    tri15_u = array(tri15[0]['ux'])
    quad8_u = array(quad8[0]['ux'])

    plot (tri15_u*1000.0/a, tri15_q/c, color=tclr,             label='Tri15'+txt)
    plot (quad8_u*1000.0/a, quad8_q/c, color=qclr, marker='+', label='Quad8'+txt)


if plt==0:
    tri15 = [read_table("abbo_sloan_01_tri15_nod_0_0.res"),
             read_table("abbo_sloan_01_tri15_nod_6_0.res"),
             read_table("abbo_sloan_01_tri15_nod_22_0.res"),
             read_table("abbo_sloan_01_tri15_nod_44_0.res"),
             read_table("abbo_sloan_01_tri15_nod_45_0.res")]
    quad8 = [read_table("abbo_sloan_01_quad8_nod_0_0.res"),
             read_table("abbo_sloan_01_quad8_nod_6_0.res"),
             read_table("abbo_sloan_01_quad8_nod_22_0.res")]

    tri15_ref = [read_table("abbo_sloan_01_tri15_ref_nod_0_0.res"),
                 read_table("abbo_sloan_01_tri15_ref_nod_6_0.res"),
                 read_table("abbo_sloan_01_tri15_ref_nod_22_0.res"),
                 read_table("abbo_sloan_01_tri15_ref_nod_44_0.res"),
                 read_table("abbo_sloan_01_tri15_ref_nod_45_0.res")]
    quad8_ref = [read_table("abbo_sloan_01_quad8_ref_nod_0_0.res"),
                 read_table("abbo_sloan_01_quad8_ref_nod_6_0.res"),
                 read_table("abbo_sloan_01_quad8_ref_nod_22_0.res")]

    plot_one (tri15_ref, quad8_ref, ' ref', 'red',    'blue')
    plot_one (tri15,     quad8,     '',     'orange', 'cyan')

    axhline (1.0174, color='k')
    text    (0,1.0174,r'$1.0174$')

    xlabel (r'$\frac{u}{a}\,\times\,10^3$')
    ylabel (r'$\frac{q}{c}$')
    legend (loc='lower right')
    grid()
    show()
