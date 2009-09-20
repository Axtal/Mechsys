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

// MechSys
#include "dem/domain.h"
#include "util/fatal.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // test the varius Domain and Particle constructors.
    Domain d;

    // add cube
    d.AddCube (Vec3_t(20,0,0),1.,15,1.,0.,NULL);
    double cube_vol = (4./3.)*M_PI-3.*M_PI*15.-pow(15,3.0)-6*pow(15,2.0);

    // add rice
    d.AddRice (Vec3_t(0,0,0),1.,10.,1.,0.,NULL);
    double rice_vol = (4./3.)*M_PI-M_PI*10.0;
    Vec3_t rice_I   = ((1./3.)*M_PI*100+(1./12.)*M_PI*1000+0.75*M_PI*10+(8./15.)*M_PI,
                       (1./3.)*M_PI*100+(1./12.)*M_PI*1000+0.75*M_PI*10+(8./15.)*M_PI,
                       0.5*M_PI*10+(8./15.)*M_PI);

    // initialize
    d.Initialize (/*dt*/0.0);

    // check
    double cube_tol_vol = 10.0;
    double rice_tol_vol = 0.1;
    double rice_tol_I   = 2.0;
    double cube_err_vol = fabs(d.Particles[0]->V - cube_vol);
    double rice_err_vol = fabs(d.Particles[1]->V - rice_vol);
    double rice_err_I   = fabs(d.Particles[1]->I(0) - rice_I(0)) +
                          fabs(d.Particles[1]->I(1) - rice_I(1)) +
                          fabs(d.Particles[1]->I(2) - rice_I(2));

    // output
    std::cout << "[1;33m\n--- Results ----------------------------------------------------[0m\n";
    cout << "  Cube Volume            = " << d.Particles[0]->V    << "  [1;31mError = " << cube_err_vol << "[0m\n";
    cout << "  Cube Center of mass    = " << d.Particles[0]->x(0) << ", " << d.Particles[0]->x(1) << ", " << d.Particles[0]->x(2) << "\n";
    cout << "  Cube Moment of inertia = " << d.Particles[0]->I(0) << ", " << d.Particles[0]->I(1) << ", " << d.Particles[0]->I(2) << "\n";
    cout << "  Cube Quaternion        = " << d.Particles[0]->Q(0) << ", " << d.Particles[0]->Q(1) << ", " << d.Particles[0]->Q(2) << ", " << d.Particles[0]->Q(3) << "\n";
    cout << endl;
    cout << "  Rice Volume            = " << d.Particles[1]->V    << "  [1;31mError = " << rice_err_vol << "[0m\n";
    cout << "  Rice Center of mass    = " << d.Particles[1]->x(0) << ", " << d.Particles[1]->x(1) << ", " << d.Particles[1]->x(2) << endl;
    cout << "  Rice Moment of inertia = " << d.Particles[1]->I(0) << ", " << d.Particles[1]->I(1) << ", " << d.Particles[1]->I(2) << "  [1;31mError =" << rice_err_I << "[0m\n";
    cout << "  Rice Quaternion        = " << d.Particles[1]->Q(0) << ", " << d.Particles[1]->Q(1) << ", " << d.Particles[1]->Q(2) << ", " << d.Particles[1]->Q(3) << "\n";

    // draw
    d.WriteBPY ("test_domain");

    // results
    if ((rice_err_vol>rice_tol_vol) || (rice_err_I>rice_tol_I) || (cube_err_vol>cube_tol_vol)) return 1;
    else return 0;
}
MECHSYS_CATCH
