# Copyright (C) 2009 Dorival M Pedroso
# Copyright (C) 2009 Sergio G Torres
# ------------------------------------
# Discrete Element Method (DEM)

from numpy        import *
from pylab        import *
from data_handler import *

class vector(matrix):
    def __getitem__(self, key):
        return matrix.__getitem__(self, (key,0))
    def norm(self):
        return sqrt((self.T*self)[0])

class Generator:
    # Constructor
    # ===========
    def __init__(self, filename):
        dat     = read_table(filename)
        np      = len(dat['Xc']) # number of particles
        self.X  = []
        self.R  = []
        for i in range(np):
            x = vector([[dat['Xc'][i]],[dat['Yc'][i]]])
            r = dat['R'][i]
            self.X.append(x)
            self.R.append(r)

    # Draw
    # ====
    def draw(self, scale=0.1, show_veloc=False):
        # colors
        lblue = (217/255.0,228/255.0,255/255.0)

        # create figure
        fig = figure()
        ax  = fig.add_subplot(111)

        # calculate bounding box
        xmin, ymin, xmax, ymax = self.X[0][0], self.X[0][1], self.X[0][0], self.X[0][1]
        for i in range(len(self.X)):
            if self.X[i][0]-self.R[i]<xmin: xmin = self.X[i][0]-self.R[i]
            if self.X[i][1]-self.R[i]<ymin: ymin = self.X[i][1]-self.R[i]
            if self.X[i][0]+self.R[i]>xmax: xmax = self.X[i][0]+self.R[i]
            if self.X[i][1]+self.R[i]>ymax: ymax = self.X[i][1]+self.R[i]

        # draw particles
        for i in range(len(self.X)):
            #r, x, v = self.R[i], self.X[i], self.V[i]
            r, x = self.R[i], self.X[i]
            ax.add_patch(matplotlib.patches.CirclePolygon(array([x[0],x[1]]),r,facecolor=lblue,edgecolor='black',linewidth=1))
            if show_veloc:
                normv = v.norm()
                dx    = scale*v[0]/normv
                dy    = scale*v[1]/normv
                arrow(x[0],x[1],dx,dy,linewidth=1,edgecolor='red')

        # scale plot
        plot ([xmin,xmax],[ymin,ymax],'o',marker='None')

# Test program
# ============
if __name__=='__main__':
    #for i in range(100):
        #print 'circles_%04d.out'%i

    for i in range(100):
        fn_in  = 'circles_%04d.out'%i
        fn_out = 'circles_%04d.png'%i
        g = Generator(fn_in)
        g.draw()
        axis('scaled')
        #show()
        savefig(fn_out)
