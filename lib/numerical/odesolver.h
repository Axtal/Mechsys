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

#ifndef MECHSYS_ODE_H
#define MECHSYS_ODE_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON
#include <stdio.h>

// GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/string.h>
#include <mechsys/util/fatal.h>

namespace Numerical
{

enum ODEScheme
{
    RK23_t,  
    RKF45_t, 
    RKDP89_t,
};

template<typename Instance>
class ODESolver
{
public:
    // Typedefs
    typedef int (Instance::*pFun) (double t, double const Y[], double dYdt[]) const; ///< Callback function

    // Constructor
    ODESolver (Instance const * p2Inst);

    // Methods
    void Evolve (pFun p2Fun, Vec_t & T, Mat_t & Y, double tf); ///< Evolve from T[0] to tf

    // Internal methods
    int CallFun (double t, double const Y[], double dYdt[]) { return (_p2inst->*_p2fun)(t, Y, dYdt); }

    // Data
    String Scheme; ///< Scheme: RK23, RKF45, RKDP89 => Classical Runge-Kutta, Runge-Kutta-Fehlberg, Runge-Kutta-Dormand-Prince
    double STOL;   ///< Local error estimate
    double EPSREL; ///< Relative error
    double DtIni;  ///< Initial step-size

private:
    Instance const * _p2inst;  ///< Pointer to an instance
    pFun             _p2fun;   ///< Pointer to instance function
};

/** Trick to pass pointers to member functions to GSL.
 * This will use (*params) to store the pointer to an instance of ODESolver, therefore (*params) will no longer be available. */
template<typename Instance>
int __ode_call_fun__ (double t, double const Y[], double dYdt[], void * not_used_for_params)
{
    ODESolver<Instance> * p2inst = static_cast<ODESolver<Instance>*>(not_used_for_params);
    return p2inst->CallFun (t, Y, dYdt);
}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Instance>
inline ODESolver<Instance>::ODESolver (Instance const * p2Inst)
    : Scheme  ("RKDP89"),
      STOL    (1.0e-6),
      EPSREL  (DBL_EPSILON),
      DtIni   (1.0e-6),
      _p2inst (p2Inst),
      _p2fun  (NULL)
{}

template<typename Instance>
inline void ODESolver<Instance>::Evolve (pFun p2Fun, Vec_t & T, Mat_t & Y, double tf)
{
    // check
    if (size(T)<2)             throw new Fatal("ODESolver::Evolve: The size of vector T must be at least equal to 2. The first item corresponds to the initial time.");
    if (Y.num_rows()<2)        throw new Fatal("ODESolver::Evolve: The number of rows in matrix Y must be at least equal to 2. The first row corresponds to the initial values.");
    if (Y.num_cols()<1)        throw new Fatal("ODESolver::Evolve: The number of columns in matrix Y must be at least equal to 1, corresponding to the number of equations.");
    if (size(T)!=Y.num_rows()) throw new Fatal("ODESolver::Evolve: The size of vector T must be equal to the number of rows in matrix Y");

    // data
    int nincs = size(T)-1;
    int neq   = Y.num_cols();

    // scheme
    gsl_odeiv_step_type const * scheme;
    if      (Scheme=="RK23")   scheme = gsl_odeiv_step_rk2;
    else if (Scheme=="RKF45")  scheme = gsl_odeiv_step_rkf45;
    else if (Scheme=="RKDP89") scheme = gsl_odeiv_step_rk8pd;
    else throw new Fatal("Scheme == %s is invalid. Valid ones are: RK23, RKF45, and RKDP89",Scheme.CStr());

    // set pointer to function
    _p2fun = p2Fun;

    // auxiliar structures
    gsl_odeiv_step    * s = gsl_odeiv_step_alloc    (scheme, neq);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (STOL, EPSREL);
    gsl_odeiv_evolve  * e = gsl_odeiv_evolve_alloc  (neq);
    gsl_odeiv_system  sys = {__ode_call_fun__<Instance>, /*jac*/NULL, neq, this};

    // initial values
    double h = DtIni;
    double t = T(0);
    double * y = new double [neq];
    for (int j=0; j<neq; ++j) y[j] = Y(0, j);

    // solve
    for (int i=1; i<=nincs; i++)
    {
        // evolve
        double ti = i * tf / static_cast<double>(nincs);
        while (t < ti)
        {
           int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, ti, &h, y);
           if (status!=GSL_SUCCESS) throw new Fatal("ODESolver::Evolve: gsl_odeiv_evolve_apply failed. status=%d",status);
           //if (h<hmin) { h = hmin; }
        }

        // results
        T(i) = t;
        for (int j=0; j<neq; ++j) Y(i, j) = y[j];
    }

    // clean up
    delete [] y;
    gsl_odeiv_evolve_free  (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free    (s);
}

}; // namespace Numerical

#endif // MECHSYS_ODE_H
