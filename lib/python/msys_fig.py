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
from numpy.linalg import norm, eig
from numpy import pi, sin, cos, tan, arcsin, arccos, arctan2, log, log10, exp, sqrt
from numpy import array, linspace, insert, repeat, zeros, matrix, ones, eye, arange, diag, dot
from numpy import logical_or, logical_and
from pylab import rcParams, gca, gcf, clf, savefig
from pylab import plot, xlabel, ylabel, show, grid, legend, subplot, axis, text, axhline, axvline, title, xticks
from pylab import contour, contourf, colorbar, clabel, xlim
from pylab import cm as MPLcm
from matplotlib.transforms   import offset_copy
from matplotlib.patches      import FancyArrowPatch, PathPatch, Polygon
from matplotlib.patches      import Arc  as MPLArc
from matplotlib.path         import Path as MPLPath
from matplotlib.font_manager import FontProperties
from matplotlib.ticker       import FuncFormatter
import mpl_toolkits.mplot3d.axes3d as A3D
from matplotlib.ticker import MaxNLocator

# to stop clipping:
# plot(..., clip_on=0)


# Set figure proportions
# ======================
def SetForEps (proport=0.75, fig_width_pt=455.24):
    # fig_width_pt = 455.24411                  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27                   # Convert pt to inch
    fig_width     = fig_width_pt*inches_per_pt  # width in inches
    fig_height    = fig_width*proport           # height in inches
    fig_size      = [fig_width,fig_height]
    params = {'mathtext.fontset':'stix', # 'cm', 'stix', 'stixsans', 'custom'
              'backend':         'ps',
              'axes.labelsize':  10,
              'text.fontsize':   10,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex':     False,
              'figure.figsize': fig_size}
    rcParams.update(params)


# Save fig with extra artists
# ===========================
def Save (filename, extra_artists=[]):
    savefig (filename, bbox_inches='tight', bbox_extra_artists=extra_artists)
    print 'File <[1;34m%s[0m> written'%filename
# As a workaround, savefig can take bbox_extra_artists keyword (this may
# only be in the svn version though), which is a list artist that needs
# to be accounted for the bounding box calculation. So in your case, the
# below code will work.
# t1 = ax.text(-0.2,0.5,'text',transform=ax.transAxes)
# fig.savefig('test.png', bbox_inches='tight', bbox_extra_artists=[t1])


# Set number of x ticks
# =====================
# ex.: ax = subplot(2,2,1); SetXnticks(ax, 4)
def SetXnticks(subplot_fig, num):  subplot_fig.xaxis.set_major_locator(MaxNLocator(num))


# Legend
# ======
def Leg (fsz=8, ncol=1): legend (loc='best',prop={'size':fsz},ncol=ncol)


# Grid
# ====
def Grid (color='grey', zorder=-100): grid (color=color, zorder=zorder)


# Grid, labels and legend
# =======================
def Gll (xl, yl, leg=True, grd=True):
    xlabel(xl)
    ylabel(yl)
    if leg: Leg()
    if grd: Grid()


# Draw cross through zero
# =======================
def Cross (x0=0.0, y0=0.0, clr='black', ls='dashed', lw=1):
    axvline(x0, color=clr, linestyle=ls, linewidth=lw)
    axhline(y0, color=clr, linestyle=ls, linewidth=lw)


# Add text
# ========
def Text (x, y, txt, x_offset=0, y_offset=0, units='points', va='bottom', ha='left', color='black', fontsize=10):
    trans = offset_copy(gca().transData, fig=gcf(), x=x_offset, y=y_offset, units=units)
    text(x, y, txt, transform=trans, va=va, ha=ha, color=color, fontsize=fontsize)


# Draw arc
# ========
def Arc (xc,yc,R, alp_min=0.0, alp_max=pi, ec='red', fc='None', lw=2, ls='solid', label=None, useArc=True, res=200, zorder=0):
    if useArc:
        gca().add_patch(MPLArc([xc,yc], 2.*R,2.*R, clip_on=False, angle=0, theta1=alp_min*180.0/pi, theta2=alp_max*180.0/pi, ls=ls, color=ec, lw=lw, label=label, zorder=zorder))
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


# Draw arrow
# ==========
def Arrow (xi,yi, xf,yf, scale=20, fc='#a2e3a2', ec='black', zorder=0):
    gca().add_patch(FancyArrowPatch((xi,yi), (xf,yf), arrowstyle='simple', mutation_scale=scale, ec=ec, fc=fc, zorder=zorder))


# Draw quad = tetragon
# ====================
def Quad (x0,y0, x1,y1, x2,y2, x3,y3, fc='#a2e3a2', ec='black', zorder=0, alpha=1.0):
    gca().add_patch(Polygon(array([[x0,y0],[x1,y1],[x2,y2],[x3,y3]]), ec=ec, fc=fc, zorder=zorder, alpha=alpha))


# Plot contour
# ============
def Contour (X,Y,Z, label, nlevels=16, cmap=None, fmt='%g', wire=True):
    c1 = contourf (X,Y,Z, cmap=cmap)
    if wire: c2 = contour (X,Y,Z, nlevels=nlevels, colors=('k'))
    cb = colorbar (c1, format=fmt)
    cb.ax.set_ylabel (label)
    if wire: clabel (c2, inline=0)


# Get ordered color
# =================
def GetClr (idx=0, scheme=1): # color
    if scheme==1:
        C = ['blue', 'green', 'cyan', 'black', 'magenta', 'orange', 'red', '#de9700', '#89009d', '#7ad473', '#737ad4', '#d473ce', '#7e6322', '#462222', '#98ac9d', '#37a3e8', 'yellow']
    else:
        C = ['blue', 'green', 'red', 'cyan', 'magenta', 'orange', 'black', '#de9700', '#89009d', '#7ad473', '#737ad4', '#d473ce', '#7e6322', '#462222', '#98ac9d', '#37a3e8', 'yellow']
    return C[idx % len(C)]


# Get ordered line style
# ======================
def GetLst (idx=0): # linestyle
    L = ['solid', 'dashed', 'dotted']
    #L = ['solid', 'dashed', 'dash_dot', 'dotted']
    return L[idx % len(L)]


# Read file with table
# ====================
# dat: dictionary with the following content:
#   dat = {'sx':[1,2,3],'ex':[0,1,2]}
#   int_cols = ['idx','num']
def read_table(filename, int_cols=[], make_maps=True):
    if not os.path.isfile(filename): raise Exception("[1;31mread_table: could not find file <[1;34m%s[0m[1;31m>[0m"%filename)
    file   = open(filename,'r')
    header = file.readline().split()
    while len(header) == 0:
        header = file.readline().split()
    dat = {}
    im  = ['id','Id','tag','Tag'] # int keys to be mapped
    fm  = ['time','Time']         # float keys to be mapped
    k2m = []                      # all keys to be mapped
    k2m.extend(im)
    k2m.extend(fm)
    int_cols.extend(im)
    for key in header:
        dat[key] = []
        if make_maps:
            if key in k2m: dat['%s2row'%key] = {}
    row = 0
    for lin in file:
        res = lin.split()
        if len(res) == 0: continue
        for i, key in enumerate(header):
            if key in int_cols: dat[key].append(int  (res[i]))
            else:               dat[key].append(float(res[i]))
            if make_maps:
                if key in k2m:
                    if dat[key][row] in dat['%s2row'%key]:
                        if type(dat['%s2row'%key][dat[key][row]]) == int:
                            dat['%s2row'%key][dat[key][row]] = [dat['%s2row'%key][dat[key][row]], row]
                        else:
                            dat['%s2row'%key][dat[key][row]].append(row)
                    else:
                        dat['%s2row'%key][dat[key][row]] = row
        row += 1
    file.close()
    mkeys = ['%s2row'%k for k in k2m]
    for key, val in dat.iteritems(): # convert lists to arrays (floats only)
        if key in int_cols or key in mkeys: continue
        dat[key] = array(val)
    return dat


# Read many files with tables
# ===========================
# filenames: list with file names
# dat: dictionary with the following content:
#   dat = {'fkey1':{'sx':[1,2,3],'ex':[0,1,2]}}
def read_tables(filenames, num_int_columns=0):
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
                if i<num_int_columns: dat[fkey][key].append(int  (res[i]))
                else:                 dat[fkey][key].append(float(res[i]))
        file.close()
    return dat


# Read temporal data
# ==================
def read_tdata (fnkey, ids, arc_lens=[]):
    r0  = read_table ('%s_%d.res'%(fnkey,ids[0]))
    str_time = 'Time'
    if r0.has_key('time'): str_time = 'time'
    if r0.has_key('t'):    str_time = 't'
    nt  = len(r0[str_time])
    np  = len(ids)
    if len(arc_lens)>0: dat = {'arc_len':zeros((np,nt))}
    else:               dat = {}
    for k, v in r0.iteritems(): dat[k] = zeros((np,nt))
    for i, n in enumerate(ids):
        r = read_table ('%s_%d.res'%(fnkey,n), make_maps=False)
        for k, v in r.iteritems(): dat[k][i,:] = v
        if len(arc_lens)>0:
            dat['arc_len'][i,:] = arc_lens[i]
    return dat


# Radians formatting
# ==================
def rad_formatting (x, pos=None):
  n = int((x/(pi/6.0))+pi/12.0)
  if n== 0: return r'$0$'
  if n== 1: return r'$\frac{\pi}{6}$'
  if n== 2: return r'$\frac{\pi}{3}$'
  if n== 3: return r'$\frac{\pi}{2}$'
  if n== 4: return r'$2\frac{\pi}{3}$'
  if n== 6: return r'$\pi$'
  if n== 8: return r'$4\frac{\pi}{3}$'
  if n== 9: return r'$3\frac{\pi}{2}$'
  if n==10: return r'$5\frac{\pi}{3}$'
  if n==12: return r'$2\pi$'
  return r'$%d\frac{\pi}{6}$'%n


# Radians and degrees formatting
# ==============================
def rad_deg_formatting (x, pos=None):
  n = int((x/(pi/6.0))+pi/12.0)
  if n== 0: return r'$0$'
  if n== 1: return r'$\frac{\pi}{6}$''\n$30^\circ$'
  if n== 2: return r'$\frac{\pi}{3}$''\n$60^\circ$'
  if n== 3: return r'$\frac{\pi}{2}$''\n$90^\circ$'
  if n== 4: return r'$2\frac{\pi}{3}$''\n$120^\circ$'
  if n== 6: return r'$\pi$''\n$180^\circ$'
  if n== 8: return r'$4\frac{\pi}{3}$''\n$240^\circ$'
  if n== 9: return r'$3\frac{\pi}{2}$''\n$270^\circ$'
  if n==10: return r'$5\frac{\pi}{3}$''\n$300^\circ$'
  if n==12: return r'$2\pi$''\n$360^\circ$'
  return r'$%d\frac{\pi}{6}$''\n$%g^\circ$'%(n,x*180.0/pi)


# Set radians formatting
# ======================
def xRadFmt (gen_ticks=False, rad_and_deg=True):
    if gen_ticks: xticks (linspace(0,2*pi,13))
    if rad_and_deg: gca().xaxis.set_major_formatter(FuncFormatter(rad_deg_formatting))
    else:           gca().xaxis.set_major_formatter(FuncFormatter(rad_formatting))


# test
# ====
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


# Column nodes
# ============
# nc: number of cells along y
def column_nodes (nc=10, o2=True):
    ny = nc + 1 # number of rows along y
    l  = arange(ny) * 2
    r  = arange(ny) * 2 + 1
    if o2:
        c = 2 * ny + arange(ny)
        m = 3 * ny + arange(nc) * 2
        L = zeros(ny + nc, dtype=int)
        i = arange(ny + nc)
        L[i%2==0] = l
        L[i%2==1] = m
        R = L + 1
        return l, c, r, L, R
    return l, r
