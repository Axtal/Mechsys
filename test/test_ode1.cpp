/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

// Std Lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// MechSys
#include <mechsys/numerical/odesolver.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;

/* Solving:
 *             x''(t) + \mu x'(t) (x(t)^2 - 1) + x(t) = 0
 * Using:
 *             x' = y
 *             y' = -x + \mu y (1-x^2)
 */


int func (double t, const double y[], double f[], void *params)
{
    double mu = *(double *)params;
    f[0] = y[1];
    f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat 
     = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix; 
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

class ODE
{
public:
    ODE (double Mu) : mu(Mu) {}
    int Fun (double t, double const Y[], double dYdt[]) const
    {
        dYdt[0] =  Y[1];
        dYdt[1] = -Y[0] - mu*Y[1]*(Y[0]*Y[0] - 1);
        return GSL_SUCCESS;
    }
    double mu;
};

int main(int argc, char **argv) try
{
    // directly calling GSL
    int nincs = 50;
    Vec_t GSL_T(nincs+1);
    Mat_t GSL_Y(nincs+1, 2);
    {
        gsl_odeiv_step_type const * T = gsl_odeiv_step_rk8pd;
        gsl_odeiv_step            * s = gsl_odeiv_step_alloc    (T, /*neq*/2);
        gsl_odeiv_control         * c = gsl_odeiv_control_y_new (/*stol*/1e-6, /*eps*/0.0);
        gsl_odeiv_evolve          * e = gsl_odeiv_evolve_alloc  (/*neq*/2);

        //printf ("step method is '%s'\n",    gsl_odeiv_step_name    (s));
        //printf ("control method is '%s'\n", gsl_odeiv_control_name (c));

        double mu = 10;
        gsl_odeiv_system sys = {func, jac, /*neq*/2, &mu};

        double t = 0.0, t1 = 100.0;
        double h = 1e-6; /*initial step-size*/
        double y[2] = { 1.0, 0.0 };
        GSL_T(0) = t;
        GSL_Y(0, 0) = 1.0;
        GSL_Y(0, 1) = 0.0;

        for (int i=1; i<=nincs; i++)
        {
            double ti = i * t1 / (double)nincs;
            while (t < ti)
            {
               int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, ti, &h, y);
               if (status!=GSL_SUCCESS) throw new Fatal("gsl_odeiv_evolve_apply failed. status=%d",status);
               //if (h<hmin) { h = hmin; }
            }
            GSL_T(i)   = t;
            GSL_Y(i,0) = y[0];
            GSL_Y(i,1) = y[1];
        }

        gsl_odeiv_evolve_free  (e);
        gsl_odeiv_control_free (c);
        gsl_odeiv_step_free    (s);
    }

    // wrapper
    Vec_t T(nincs+1);
    Mat_t Y(nincs+1, 2);
    {
        ODE ode(/*mu*/10.0);
        Numerical::ODESolver<ODE> Solver(&ode);
        //Solver.Scheme = "RKF45";
        T(0)    = 0.0;
        Y(0, 0) = 1.0;
        Y(0, 1) = 0.0;
        Solver.Evolve (&ODE::Fun, T, Y, /*tf*/100.0);
    }

    // check
    for (int i=0; i<nincs+1; ++i)
    {
        printf("t=%5.2e, %5.2e, d=%5.2e    y0=% 10.4e, % 10.4e, d=% 10.4e    y1=% 10.4e, % 10.4e, d=% 10.4e\n", 
                T(i),   GSL_T(i),   T(i)  -GSL_T(i),
                Y(i,0), GSL_Y(i,0), Y(i,0)-GSL_Y(i,0),
                Y(i,1), GSL_Y(i,1), Y(i,1)-GSL_Y(i,1));
    }

    return 0;
}
MECHSYS_CATCH
