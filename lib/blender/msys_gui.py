########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raul D. D. Farfan             #
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

import Blender
from   Blender import Draw, BGL

def caption1_clr():
    #BGL.glColor3f (0.4, 0.4, 0.4)
    BGL.glColor3f (0.573, 0.512, 0.585)

def caption2_clr():
    #BGL.glColor3f (0.431, 0.443, 0.514)
    BGL.glColor3f (0.573, 0.512, 0.585)

def caption3_clr():
    #BGL.glColor3f (0.53, 0.54, 0.6)
    BGL.glColor3f (0.573, 0.512, 0.585)

def box1_clr():
    #BGL.glColor3f (0.85, 0.85, 0.85)
    #BGL.glColor3f (0.553, 0.553, 0.718)
    BGL.glColor3f (0.656, 0.696, 0.617)

def box2_clr():
    #BGL.glColor3f (0.72, 0.72, 0.8)
    BGL.glColor3f (0.667, 0.667, 0.729)

def box3_clr():
    #BGL.glColor3f (0.82, 0.82, 0.9)
    BGL.glColor3f (0.7, 0.75, 0.7)


def label(txt, c,r,w,h):
    BGL.glColor3f     (0.663, 0.663, 0.663)
    BGL.glRecti       (c, r, c+w, r+h)
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         (txt)
    # frame
    BGL.glColor3f  (0.3, 0.3, 0.3)
    BGL.glBegin    (BGL.GL_LINE_LOOP)
    BGL.glVertex2i (c,   r)
    BGL.glVertex2i (c+w, r)
    BGL.glVertex2i (c+w, r+h)
    BGL.glVertex2i (c,   r+h)
    BGL.glEnd      ()

def label_(txt, c,r,w,h):
    l = len(txt)*10
    BGL.glColor3f     (0.663, 0.663, 0.663)
    BGL.glRecti       (c, r, c+w, r+h)
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (c+5+(w-l)/2, r+2+(h-12)/2)
    Draw.Text         (txt)
    # frame
    BGL.glColor3f  (0.3, 0.3, 0.3)
    BGL.glBegin    (BGL.GL_LINE_LOOP)
    BGL.glVertex2i (c,   r)
    BGL.glVertex2i (c+w, r)
    BGL.glVertex2i (c+w, r+h)
    BGL.glVertex2i (c,   r+h)
    BGL.glEnd      ()

def text(c,r, txt):
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         (txt)

def background():
    #BGL.glClearColor (0.531, 0.543, 0.614, 0.0)
    BGL.glClearColor (0.725, 0.725, 0.725, 0.0)
    BGL.glClear      (Blender.BGL.GL_COLOR_BUFFER_BIT)

def caption1(c,r,w,rh, lab, evt_refresh,evt_hideall,evt_showonly,evt_showhide):
    caption1_clr      ()
    BGL.glRecti       (c, r, c+w, r+int(1.3*rh))
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         (lab)
    Draw.PushButton   ('Refresh',   evt_refresh,  c+w-5-70-70-60-60, r+3, 60, rh-4, 'Redraw all windows')
    Draw.PushButton   ('Hide All',  evt_hideall,  c+w-5-70-70-60,    r+3, 60, rh-4, 'Hide all boxes')
    Draw.PushButton   ('Show Only', evt_showonly, c+w-5-70-70,       r+3, 70, rh-4, 'Show only this box')
    Draw.PushButton   ('Show/Hide', evt_showhide, c+w-5-70,          r+3, 70, rh-4, 'Show/Hide this box')

def box1_in(W,cg,rh, c,r,w,h):
    box1_clr      ()
    BGL.glRecti   (c,r-h,c+w,r)
    return r-2*rh, c+cg, W-2*(c+cg)

def box1_out(W,cg,rh, c,r):
    return r-rh, c-cg, W-2*(c-cg)

def caption2(c,r,w,rh, lab):
    caption2_clr      ()
    BGL.glRecti       (c, r, c+w, r+rh)
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         (lab)

def caption2_(c,r,w,rh, lab, evt_add,evt_delall):
    caption2_clr      ()
    BGL.glRecti       (c, r, c+w, r+rh)
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         (lab)
    Draw.PushButton   ('Add',        evt_add,    c+w-5-60-80, r+2, 60, rh-4, 'Add ')
    Draw.PushButton   ('Delete all', evt_delall, c+w-5-80,    r+2, 80, rh-4, 'Delete all ')

def caption2__(c,r,w,rh, lab, evt_add,evt_del,evt_delall):
    caption2_clr      ()
    BGL.glRecti       (c, r, c+w, r+rh)
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         (lab)
    Draw.PushButton   ('Add',        evt_add,    c+w-5-60-60-80, r+2, 60, rh-4, 'Add ')
    Draw.PushButton   ('Delete',     evt_del,    c+w-5-60-80,    r+2, 60, rh-4, 'Delete')
    Draw.PushButton   ('Delete all', evt_delall, c+w-5-80,       r+2, 80, rh-4, 'Delete all ')

def box2_in(W,cg,rh, c,r,w,h):
    box2_clr      ()
    BGL.glRecti   (c,r-h,c+w,r)
    return r-2*rh, c+cg, W-2*(c+cg)

def box2__in(W,cg,rh, c,r,w,h):
    box2_clr      ()
    BGL.glRecti   (c,r-h,c+w,r)
    return r-rh, c+cg, W-2*(c+cg)

def box2_out(W,cg,rh, c,r):
    return r-2*rh, c-cg, W-2*(c-cg)

def box2__out(W,cg,rh, c,r):
    return r-rh, c-cg, W-2*(c-cg)

def caption3(c,r,w,rh, lab, evt_add,evt_delall):
    caption3_clr      ()
    BGL.glRecti       (c, r, c+w, r+rh)
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         (lab)
    Draw.PushButton   ('Add',        evt_add,    c+w-5-60-80, r+2, 60, rh-4, 'Add ')
    Draw.PushButton   ('Delete all', evt_delall, c+w-5-80,    r+2, 80, rh-4, 'Delete all ')

def caption3_(c,r,w,rh, lab):
    caption3_clr      ()
    BGL.glRecti       (c, r, c+w, r+rh)
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         (lab)

def box3_in(W,cg,rh, c,r,w,h):
    box3_clr      ()
    BGL.glRecti   (c,r-h,c+w,r)
    return r-rh, c+cg, W-2*(c+cg)

def box3_out(W,cg,rh, c,r):
    return r-rh, c-cg, W-2*(c-cg)
