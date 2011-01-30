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

from mechsys  import *
from msys_fig import *

# Solving:
#             x''(t) + \mu x'(t) (x(t)^2 - 1) + x(t) = 0
# Using:
#             x' = y
#             y' = -x + \mu y (1-x^2)

class ODE:
    def __init__ (self, mu): self.mu = mu

    def fun (self, t, Y, dYdt):
        dYdt[0] =  Y[1]
        dYdt[1] = -Y[0] - self.mu*Y[1]*(Y[0]*Y[0] - 1)

# data
mu    = 2.0
y0ini = 0.0
y1ini = 3.0
tf    = 10.209022
dtOut = 0.1

# solve
ode = ODE(mu)
sol = ODESolver(ode, 'ODE', 'fun')
y   = [y0ini, y1ini]
sol.Init (0.0, y, "RK12", 1.e-2)
sol.EvolveOut (tf, dtOut, "py_test_ode1.dat")
