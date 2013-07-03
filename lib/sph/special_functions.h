/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_SPH_SPECIAL_H
#define MECHSYS_SPH_SPECIAL_H

// Std lib
#include <iostream>

// MechSys
#include <mechsys/linalg/matvec.h>

namespace SPH {

inline double Kernel(double r,double h)
{
    double C = 1.0/(h*h*h*M_PI);
    double q = r/h;
    if ((q>=0.0)&&(q<1)) return C*(1-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q);
    else if (q<=2)       return C*((1.0/4.0)*(2-q)*(2-q)*(2-q));
    else                 return 0.0;
}

inline double GradKernel(double r, double h)
{
    double C = 1.0/(h*h*h*M_PI);
    double q = r/h;
    if ((q>=0.0)&&(q<1)) return C*(3.0*q+(9.0/4.0)*q*q);
    else if (q<=2)       return C*(-(3.0/4.0)*(2-q)*(2-q));
    else                 return 0.0;
}

inline double Kernel(Vec3_t & x, Vec3_t & xp, double h)
{
    double r = norm(xp-x);
    return Kernel(r,h);
}

inline double Pressure(double rho)
{
    double P0 = 10.0;
    double rho0 = 10.0;
    return P0*(pow(rho/rho0,7)-1);
}

inline double SoundSpeed(double rho)
{
    double P0 = 10.0;
    double rho0 = 10.0;
    return sqrt(7*P0*(pow(rho/rho0,6)/rho0));
}

}; // namespace SPH

#endif // MECHSYS_SPH_SPECIAL_H
