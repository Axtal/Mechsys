from numpy import pi, sin, cos, tan, arcsin, arccos, arctan, log, log10, exp
from numpy import array, linspace, insert, repeat, zeros
from pylab import rcParams, gca, gcf, clf, savefig
from pylab import plot, xlabel, ylabel, show, grid, legend, subplot, axis, text, axhline, axvline, title
from matplotlib.transforms import offset_copy
from matplotlib.patches    import FancyArrowPatch, PathPatch
from matplotlib.patches    import Arc  as MPLArc
from matplotlib.path       import Path as MPLPath

def SetForEps (proport=0.75, fig_width_pt=455.24):
    # fig_width_pt = 455.24411                  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27                   # Convert pt to inch
    fig_width     = fig_width_pt*inches_per_pt  # width in inches
    fig_height    = fig_width*proport           # height in inches
    fig_size      = [fig_width,fig_height]
    params = {'backend':         'ps',
              'axes.labelsize':  10,
              'text.fontsize':   10,
              'legend.fontsize': 9,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex':     False,
              'figure.figsize': fig_size}
    rcParams.update(params)

def Save (filename): savefig (filename, bbox_inches='tight')

def Grid (color='grey', zorder=-100): grid (color=color, zorder=zorder)

def Text (x, y, txt, x_offset=0, y_offset=0, units='points', va='bottom', ha='left', color='black', fontsize=10):
    trans = offset_copy(gca().transData, fig=gcf(), x=x_offset, y=y_offset, units=units)
    text(x, y, txt, transform=trans, va=va, ha=ha, color=color, fontsize=fontsize)

def Arc (xc,yc,R, alp_min=0.0, alp_max=pi, ec='red', fc='None', lw=2, ls='solid', label=None, useArc=True, res=200, zorder=0):
    if useArc:
        gca().add_patch(MPLArc([xc,yc], 2.*R,2.*R, angle=0, theta1=alp_min, theta2=alp_max*180.0/pi, ls=ls, color=ec, lw=lw, label=label, zorder=zorder))
    else:
        A   = linspace(alp_min,alp_max,res)
        dat = [(MPLPath.MOVETO, (xc+R*cos(alp_min),yc+R*sin(alp_min)))]
        for a in A[1:]:
            x = xc + R*cos(a)
            y = yc + R*sin(a)
            dat.append ((MPLPath.LINETO, (x,y)))
        cmd,vert = zip(*dat)
        ph = MPLPath (vert, cmd)
        gca().add_patch(PathPatch(ph, fc=fc, ec=ec, linewidth=lw, linestyle=ls, label=label, zorder=zorder))

def Arrow (xi,yi, xf,yf, scale=20, fc='#a2e3a2', ec='black', zorder=0):
    gca().add_patch(FancyArrowPatch((xi,yi), (xf,yf), arrowstyle='simple', mutation_scale=scale, ec=ec, fc=fc, zorder=zorder))

if __name__=='__main__':
    SetForEps ()
    x = linspace (0, 10, 100)
    y = x**1.5
    plot    (x,y, 'b-', label='sim')
    Arc     (0,0,10, useArc=False)
    Arc     (0,0,20, ec='magenta', useArc=True)
    Arrow   (-10,0,max(x),max(y))
    Text    (-20,25,r'$\sigma$')
    Text    (-20,25,r'$\sigma$',y_offset=-10)
    axvline (0,color='black',zorder=-1)
    axhline (0,color='black',zorder=-1)
    Grid    ()
    axis    ('equal')
    legend  (loc='upper left')
    xlabel  (r'$x$')
    ylabel  (r'$y$')
    #show    ()
    Save    ('test_msys_fig.eps')
