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
#include <mechsys/util/numstreams.h>

using std::cout;
using std::endl;
using Util::_8s;
using Util::PI;

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
    int Fun (double t, double const Y[], double dYdt[])
    {
        dYdt[0] =  Y[1];
        dYdt[1] = -Y[0] - mu*Y[1]*(Y[0]*Y[0] - 1);
        return GSL_SUCCESS;
    }
    int TgIncs (double t, double const Y[], double dt, double dY[])
    {
        double dYdt0 =  Y[1];
        double dYdt1 = -Y[0] - mu*Y[1]*(Y[0]*Y[0] - 1);
        dY[0] = dYdt0*dt;
        dY[1] = dYdt1*dt;
        return GSL_SUCCESS;
    }
    double mu;
};

int main(int argc, char **argv) try
{
    // number of output increments
    int    nincs = 100;
    double mu    = 2.0;
    double tf    = 10.209022;
    double y0ini = 0.0;
    double y1ini = 3.0;

    // GSL
    Vec_t T(nincs+1);
    Mat_t Y(nincs+1, 2);
    {
        ODE ode(mu);
        Numerical::ODESolver<ODE> sol(&ode);
        //sol.Scheme = "RKF45";
        T(0)    = 0.0;
        Y(0, 0) = y0ini;
        Y(0, 1) = y1ini;
        sol.Evolve (&ODE::Fun, T, Y, tf);
    }

    // ME
    Vec_t Tme(nincs+1);
    Mat_t Yme(nincs+1, 2);
    {
        ODE ode(mu);
        Numerical::ODESolver<ODE> sol(&ode);
        //sol.STOL  = 1.0e-8;
        //sol.MaxSS = 100000;
        Tme(0)    = 0.0;
        Yme(0, 0) = y0ini;
        Yme(0, 1) = y1ini;
        sol.Evolve (&ODE::TgIncs, Tme, Yme, tf);
    }

    // check
    Array<double> err_t (nincs+1);
    Array<double> err_y0(nincs+1);
    Array<double> err_y1(nincs+1);
    for (int i=0; i<nincs+1; ++i)
    {
        err_t [i] = fabs(T(i)  -Tme(i));
        err_y0[i] = fabs(Y(i,0)-Yme(i,0));
        err_y1[i] = fabs(Y(i,1)-Yme(i,1));
    }
    double tol_t  = 1.0e-13;
    double tol_y0 = 1.0e-5;
    double tol_y1 = 1.0e-4;
    printf("\n%s--- Error ------------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    cout << TERM_CLR2 << Util::_4<< "   " << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm" << TERM_RST << endl;
    cout << Util::_4<< "t"  << Util::_8s<<err_t .TheMin() << Util::_8s<<err_t .Mean() << (err_t .TheMax()>tol_t  ? TERM_RED : TERM_GREEN) << Util::_8s<<err_t .TheMax() << TERM_RST << Util::_8s<<err_t .Norm() << "\n";
    cout << Util::_4<< "y0" << Util::_8s<<err_y0.TheMin() << Util::_8s<<err_y0.Mean() << (err_y0.TheMax()>tol_y0 ? TERM_RED : TERM_GREEN) << Util::_8s<<err_y0.TheMax() << TERM_RST << Util::_8s<<err_y0.Norm() << "\n";
    cout << Util::_4<< "y1" << Util::_8s<<err_y1.TheMin() << Util::_8s<<err_y1.Mean() << (err_y1.TheMax()>tol_y1 ? TERM_RED : TERM_GREEN) << Util::_8s<<err_y1.TheMax() << TERM_RST << Util::_8s<<err_y1.Norm() << "\n";
    cout << endl;

    // output file
    std::ofstream of("test_ode1.dat", std::ios::out);
    of << _8s<<"t" << _8s<<"y0" << _8s<<"y1" << _8s<<"y0me" << _8s<<"y1me" << endl;
    for (int i=0; i<nincs+1; ++i) of << _8s<<T(i) << _8s<<Y(i,0) << _8s<<Y(i,1) << _8s<<Yme(i,0) << _8s<<Yme(i,1) << endl;
    of.close();
    cout << "File <" << TERM_CLR_BLUE << "test_ode1.dat" << TERM_RST << "> written" << endl;
    cout << endl;

    // end
    return 0;
}
MECHSYS_CATCH
