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

from numpy import array, sqrt, linspace, pi, cos, sin
from pylab import figure, text, show, axis
from pylab import matplotlib as MPL

class DrawMesh:

    # Constructor
    # ============
    #
    # V: vertices
    # C: cells
    #
    # V=[[ 0, -1,   0.0,  0.0],
    #    [ 1, -2,   1.0,  0.0],
    #    [ 2, -4,   0.0,  1.0],
    #    [ 3, -3,   1.0,  1.0]]
    #
    # C=[[ 0, -1, [0,1,3,2], {0:-10,1:-20,2:-30,3:-40}, {}]]
    #
    # pct: percentage of drawing limits to use for icons
    def __init__(self, V,C, Pins={}, Shares={}, pct=0.001, fsz1=8, fsz2=6):
        # mesh
        self.V      = V
        self.C      = C
        self.Pins   = Pins
        self.Shares = Shares
        self.ndim   = len(V[0])-2
        if self.ndim>2: raise Exception("Drawing: NDim=%d is invalid"%self.ndim)

        # constants
        self.pct  = pct
        self.fsz1 = fsz1
        self.fsz2 = fsz2

        # colors
        self.pink    = (250/255.0,204/255.0,228/255.0)
        self.lblue   = (217/255.0,228/255.0,255/255.0)
        self.lgreen  = (100/255.0,241/255.0,193/255.0)
        self.dblue   = ( 45/255.0,  0/255.0,160/255.0)
        self.orange  = (241/255.0,125/255.0,  0/255.0)
        self.lyellow = (234/255.0,228/255.0,179/255.0)
        #self.purple  = '#9b8de3'
        self.purple  = '#c5a9f3'
        self.dred    = '#b30000'

        # assign colors
        self.celledgeclr = self.dblue
        self.lineedgeclr = self.dred

        # drawing limits (bounding box)
        allx      = [v[2] for v in self.V]
        if self.ndim==1: ally = [0.0, 1.0]
        else:            ally = [v[3] for v in self.V]
        self.lims = array([min(allx),max(allx),min(ally),max(ally)])
        self.diag = sqrt((self.lims[1]-self.lims[0])**2.0+(self.lims[3]-self.lims[2])**2.0)

        # noise to move tags and ids
        self.yidnoise = 0.001*self.diag

        # icon's length
        self.l = self.pct*self.diag

        # matplotlib's structures
        self.PH = MPL.path.Path
        self.PP = MPL.patches
        self.PC = MPL.patches.PathPatch

    # Draw mesh
    #==========
    def draw(self, with_tags=True, with_ids=True, with_shares=True, only_lin_cells=False):
        # create figure
        fig = figure()
        ax  = fig.add_subplot(111)
        ax.grid()

        # draw points at bounding box
        dlim  = sqrt((self.lims[1]-self.lims[0])**2.0+(self.lims[3]-self.lims[2])**2.0)
        limsx = self.lims + array([-self.pct*dlim,+self.pct*dlim,-self.pct*dlim,+self.pct*dlim]) # extended limits
        ax.plot(limsx[:2],limsx[2:],'o',marker='None')

        # draw solid cells
        dat = []
        if not only_lin_cells:
            for c in self.C:
                con  = c[2] # connectivity
                nnod = len(con)
                if nnod>2:
                    x0 = self.V[con[0]][2]
                    y0 = self.V[con[0]][3]
                    dat.append((self.PH.MOVETO, (x0,y0)))
                    # edges for 3,4
                    if nnod<=4:
                        for j in range(1,nnod):
                            xj = self.V[con[j]][2]
                            yj = self.V[con[j]][3]
                            dat.append((self.PH.LINETO, (xj,yj)))
                    # edges for 6,8
                    if nnod==6:
                        dat.append((self.PH.LINETO, (self.V[con[3]][2], self.V[con[3]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[1]][2], self.V[con[1]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[4]][2], self.V[con[4]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[2]][2], self.V[con[2]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[5]][2], self.V[con[5]][3])))
                    if nnod==8:
                        dat.append((self.PH.LINETO, (self.V[con[4]][2], self.V[con[4]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[1]][2], self.V[con[1]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[5]][2], self.V[con[5]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[2]][2], self.V[con[2]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[6]][2], self.V[con[6]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[3]][2], self.V[con[3]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[7]][2], self.V[con[7]][3])))
                    if nnod==15:
                        dat.append((self.PH.LINETO, (self.V[con[ 6]][2], self.V[con[ 6]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 3]][2], self.V[con[ 3]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 7]][2], self.V[con[ 7]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 1]][2], self.V[con[ 1]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 8]][2], self.V[con[ 8]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 4]][2], self.V[con[ 4]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 9]][2], self.V[con[ 9]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 2]][2], self.V[con[ 2]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[10]][2], self.V[con[10]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 5]][2], self.V[con[ 5]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[11]][2], self.V[con[11]][3])))
                    dat.append((self.PH.CLOSEPOLY, (0,0)))
            if len(dat)>0:
                cmd,vert = zip(*dat)
                ph0 = self.PH (vert, cmd)
                pc0 = self.PC (ph0, facecolor=self.lblue, edgecolor=self.celledgeclr, linewidth=2)
                ax.add_patch  (pc0)

        # draw linear cells
        if self.ndim==1:
            for c in self.C:
                con = c[2] # connectivity
                if len(con)==2:
                    x0 = self.V[con[0]][2]
                    y0 = 0.0
                    x1 = self.V[con[1]][2]
                    y1 = 0.0
                    XY = array([[x0,y0],[x1,y1]])
                    ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=self.celledgeclr, lw=2))
        else:
            for c in self.C:
                con = c[2] # connectivity
                if len(con)==2:
                    x0 = self.V[con[0]][2]
                    y0 = self.V[con[0]][3]
                    x1 = self.V[con[1]][2]
                    y1 = self.V[con[1]][3]
                    XY = array([[x0,y0],[x1,y1]])
                    ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=self.lineedgeclr, lw=4))

        # text
        if with_ids or with_tags:
            for c in self.C:
                # centre
                con = c[2] # connectivity
                nnod = len(con)
                if only_lin_cells and nnod>2: continue
                if nnod==2:
                    if self.ndim==1:
                        x0 = self.V[con[0]][2]
                        y0 = 0.0
                        x1 = self.V[con[1]][2]
                        y1 = 0.0
                        xc = x0 + 0.5*(x1-x0)
                        yc = y0 + 0.5*(y1-y0)
                    else:
                        x0 = self.V[con[0]][2]
                        y0 = self.V[con[0]][3]
                        x1 = self.V[con[1]][2]
                        y1 = self.V[con[1]][3]
                        cf = 0.5 if only_lin_cells else 0.3
                        xc = x0 + cf*(x1-x0)
                        yc = y0 + cf*(y1-y0)
                else:
                    xc = self.V[con[0]][2]
                    yc = self.V[con[0]][3]
                    for j in range(1,nnod):
                        xc += self.V[con[j]][2]
                        yc += self.V[con[j]][3]
                    xc = xc/nnod
                    yc = yc/nnod
                # text
                txt = '   '
                if with_tags: txt = '%s%d' % (txt,c[1])
                if len(c)>4:
                    if len(c[4])>0: txt += '\n   '
                    for brykey, side_id in c[4].iteritems():
                        # neighbours
                        side = side_id[0]
                        eid  = side_id[1]
                        txt  = '%s %d' % (txt,eid)
                if len(txt)>0: ax.text(xc,yc, txt,  ha='left',  backgroundcolor=self.lgreen,fontsize=self.fsz2)
                if with_ids:   ax.text(xc,yc, c[0], ha='right', backgroundcolor=self.purple,fontsize=self.fsz1)
                # edge tags
                if with_tags and self.ndim>1 and len(c)>3:
                    nnod = len(con)
                    for side, tag in c[3].iteritems():
                        if tag<0: # has tag
                            if nnod==3:
                                if   side==0: na, nb = con[0], con[1]
                                elif side==1: na, nb = con[1], con[2]
                                elif side==2: na, nb = con[2], con[0]
                            elif nnod==6:
                                if   side==0: na, nb,  nc, nd = con[0], con[3],  con[3], con[1]
                                elif side==1: na, nb,  nc, nd = con[1], con[4],  con[4], con[2]
                                elif side==2: na, nb,  nc, nd = con[2], con[5],  con[5], con[0]
                            elif nnod==4:
                                if   side==0: na, nb = con[0], con[1]
                                elif side==1: na, nb = con[1], con[2]
                                elif side==2: na, nb = con[2], con[3]
                                elif side==3: na, nb = con[3], con[0]
                            elif nnod==8:
                                if   side==0: na, nb,  nc, nd = con[0], con[4],  con[4], con[1]
                                elif side==1: na, nb,  nc, nd = con[1], con[5],  con[5], con[2]
                                elif side==2: na, nb,  nc, nd = con[2], con[6],  con[6], con[3]
                                elif side==3: na, nb,  nc, nd = con[3], con[7],  con[7], con[0]
                            elif nnod==15:
                                if   side==0: na,nb, nc,nd, ne,nf, ng,nh = con[0],con[ 6], con[ 6],con[3], con[3],con[ 7], con[ 7],con[1]
                                elif side==1: na,nb, nc,nd, ne,nf, ng,nh = con[1],con[ 8], con[ 8],con[4], con[4],con[ 9], con[ 9],con[2]
                                elif side==2: na,nb, nc,nd, ne,nf, ng,nh = con[2],con[10], con[10],con[5], con[5],con[11], con[11],con[0]
                            xa = self.V[na][2]
                            ya = self.V[na][3]
                            xb = self.V[nb][2]
                            yb = self.V[nb][3]
                            xc = (xa+xb)/2.0
                            yc = (ya+yb)/2.0
                            ax.text(xc,yc, '(%d)'%tag, ha='center', va='center', fontsize=self.fsz2, backgroundcolor=self.pink)
                            if nnod==6 or nnod==8 or nnod==15:
                                xa = self.V[nc][2]
                                ya = self.V[nc][3]
                                xb = self.V[nd][2]
                                yb = self.V[nd][3]
                                xc = (xa+xb)/2.0
                                yc = (ya+yb)/2.0
                                ax.text(xc,yc, '(%d)'%tag, ha='center', va='center', fontsize=self.fsz2, backgroundcolor=self.pink)
                            if nnod==15:
                                xa = self.V[ne][2]
                                ya = self.V[ne][3]
                                xb = self.V[nf][2]
                                yb = self.V[nf][3]
                                xc = (xa+xb)/2.0
                                yc = (ya+yb)/2.0
                                ax.text(xc,yc, '(%d)'%tag, ha='center', va='center', fontsize=self.fsz2, backgroundcolor=self.pink)
                                xa = self.V[ng][2]
                                ya = self.V[ng][3]
                                xb = self.V[nh][2]
                                yb = self.V[nh][3]
                                xc = (xa+xb)/2.0
                                yc = (ya+yb)/2.0
                                ax.text(xc,yc, '(%d)'%tag, ha='center', va='center', fontsize=self.fsz2, backgroundcolor=self.pink)

        #self.circle(0,0, 5.0,pi/2.0,ax)
        #self.circle(0,0,15.0,pi/2.0,ax)

        # draw nodes

        #prob = [ 455, 464, 468, 471, 478, 480, 483, 492, 494, 498, 2446, 2461, 2568, 2610, 2954, 3368, 3371, 3381, 3384]

        if with_ids:
            for v in self.V:

                #if not v[0] in prob: continue

                # skip if not connected to lin cell
                if only_lin_cells and with_shares:
                    con_lin_cell = False # connected to lin cell ?
                    if len(self.Shares)<1: continue
                    for e in self.Shares[v[0]]:
                        if len(self.C[e][2])==2: con_lin_cell = True
                    if not con_lin_cell: continue
                # shares
                if v[0] in self.Shares and with_shares:
                    s    = '('
                    eles = self.Shares[v[0]]
                    for i, e in enumerate(eles):
                        if i==len(eles)-1: s += '%d) .' % e
                        else:              s += '%d,' % e
                    text(v[2], v[3], s, va='bottom', ha='right', color='black', backgroundcolor='white', fontsize=self.fsz2)
                # ids
                s = '%d' % v[0]
                yval = v[3]+self.yidnoise if self.ndim>1 else self.yidnoise
                text(v[2], yval, s, va='bottom', ha='right', color='black', backgroundcolor=self.lyellow, fontsize=self.fsz1)
                tag = v[1]
                # tag
                if tag<0 and with_tags:
                    yval = v[3]-self.yidnoise if self.ndim>1 else -self.yidnoise
                    text(v[2], yval, '%d'%tag, va='top', ha='left', color='black', backgroundcolor=self.orange, fontsize=self.fsz2)

        # pins
        for key, val in self.Pins.iteritems():
            # ids
            x = self.V[key][2]
            y = self.V[key][3]
            s = '(%d,' % key
            for i in range(len(val)):
                if i==len(val)-1: s += '%d)' % val[i]
                else:             s += '%d,' % val[i]
            text(x, y, s, va='bottom', ha='left', backgroundcolor=self.lyellow, fontsize=self.fsz1)
            # tag
            tag = self.V[key][1]
            if tag<0 and with_tags:
                text(x, y, '%d'%tag, va='top', color='black', backgroundcolor='none', fontsize=self.fsz2)

    # Show figure
    # ===========
    def show(self):
        axis('equal')
        show()

    # Add circle
    # ==========
    def circle(self,xc,yc,R,alp_max,ax):
        A   = linspace(0.0,alp_max,200)
        dat = [(self.PH.MOVETO, (R,0.0))]
        for a in A:
            x = xc + R*cos(a)
            y = yc + R*sin(a)
            dat.append((self.PH.LINETO, (x,y)))
        cmd,vert = zip(*dat)
        ph = self.PH (vert, cmd)
        pc = self.PC (ph, facecolor=self.lyellow, edgecolor="red", linewidth=4)
        ax.add_patch (pc)
