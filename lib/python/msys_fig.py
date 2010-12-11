########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Dorival M. Pedroso                                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

import os.path
from os.path import basename
from numpy import pi, sin, cos, tan, arcsin, arccos, arctan, log, log10, exp, sqrt
from numpy import array, linspace, insert, repeat, zeros, matrix
from pylab import rcParams, gca, gcf, clf, savefig
from pylab import plot, xlabel, ylabel, show, grid, legend, subplot, axis, text, axhline, axvline, title
from pylab import contour, contourf, colorbar, clabel
from pylab import cm as MPLcm
from matplotlib.transforms   import offset_copy
from matplotlib.patches      import FancyArrowPatch, PathPatch
from matplotlib.patches      import Arc  as MPLArc
from matplotlib.path         import Path as MPLPath
from matplotlib.font_manager import FontProperties

def SetForEps (proport=0.75, fig_width_pt=455.24):
    # fig_width_pt = 455.24411                  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27                   # Convert pt to inch
    fig_width     = fig_width_pt*inches_per_pt  # width in inches
    fig_height    = fig_width*proport           # height in inches
    fig_size      = [fig_width,fig_height]
    params = {'backend':         'ps',
              'axes.labelsize':  10,
              'text.fontsize':   10,
              'legend.fontsize': 8,
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

def Contour (X,Y,Z, label, nlevels=16, cmap=None, fmt='%g'):
    c1 = contourf (X,Y,Z, cmap=cmap)
    c2 = contour  (X,Y,Z, nlevels=nlevels, colors=('k'))
    cb = colorbar (c1, format=fmt)
    cb.ax.set_ylabel (label)
    clabel (c2, inline=0)

def GetClr (idx=0): # color
    C = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', '#de9700', '#89009d', '#7ad473', '#737ad4', '#d473ce', '#7e6322', '#462222', '#98ac9d', '#37a3e8']
    return C[idx % len(C)]

def GetLst (idx=0): # linestyle
    L = ['solid', 'dashed', 'dash_dot', 'dotted']
    return L[idx % len(L)]

# Read file with table
# ====================
# dat: dictionary with the following content:
#   dat = {'sx':[1,2,3],'ex':[0,1,2]}
def read_table(filename):
    if not os.path.isfile(filename): raise Exception("[1;31mread_table: could not find file <[1;34m%s[0m[1;31m>[0m"%filename)
    file   = open(filename,'r')
    header = file.readline().split()
    dat    = {}
    for key in header: dat[key] = []
    for lin in file:
        res = lin.split()
        for i, key in enumerate(header):
            dat[key].append(float(res[i]))
    file.close()
    for k, v in dat.iteritems():
        dat[k] = array(v)
    return dat


# Read many files with tables
# ===========================
# filenames: list with file names
# dat: dictionary with the following content:
#   dat = {'fkey1':{'sx':[1,2,3],'ex':[0,1,2]}}
def read_tables(filenames):
    dat = {}
    for fn in filenames:
        if not os.path.isfile(fn): raise Exception("[1;31mread_tables: could not find file <[1;34m%s[0m[1;31m>[0m"%filename)
        fkey      = basename(fn).replace('.dat','')
        file      = open(fn,'r')
        header    = file.readline().split()
        dat[fkey] = {}
        for key in header: dat[fkey][key] = []
        for lin in file:
            res = lin.split()
            for i, key in enumerate(header):
                dat[fkey][key].append(float(res[i]))
        file.close()
    return dat


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