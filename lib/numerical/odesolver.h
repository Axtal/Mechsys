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

template<typename Instance>
class ODESolver
{
public:
    // Typedefs
    typedef int (Instance::*pFun) (double t, double const Y[], double dYdt[]) const; ///< Callback function
    typedef int (Instance::*pTgIncs) (double t, double const Y[], double dt, double dY[]) const; ///< Callback function

    // Constructor
    ODESolver (Instance const * p2Inst);

    // Methods
    void Evolve (pFun p2Fun, Vec_t & T, Mat_t & Y, double tf); ///< Evolve from T[0] to tf
    void Evolve (pTgIncs p2Incs, Vec_t & T, Mat_t & Y, double tf); ///< Evolve from T[0] to tf

    // Internal methods
    int CallFun (double t, double const Y[], double dYdt[]) { return (_p2inst->*_p2fun)(t, Y, dYdt); }

    // Data
    String Scheme; ///< Scheme: RK23, RKF45, RKDP89 => Classical Runge-Kutta, Runge-Kutta-Fehlberg, Runge-Kutta-Dormand-Prince
    double STOL;   ///< Local error estimate
    double EPSREL; ///< Relative error
    double DtIni;  ///< Initial step-size
    double dTini;  ///< Initial pseudo time step-size
    size_t MaxSS;  ///< Max number of sub-steps
    double mMin;   ///< Minimum change of substep size
    double mMax;   ///< Maximum change of substep size

private:
    Instance const * _p2inst;  ///< Pointer to an instance
    pFun             _p2fun;   ///< Pointer to instance function
    pTgIncs          _p2incs;  ///< Pointer to instance function
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
      dTini   (1.0e-3),
      MaxSS   (2000),
      mMin    (0.1),
      mMax    (10.0),
      _p2inst (p2Inst),
      _p2fun  (NULL),
      _p2incs (NULL)
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
    double Dt = (tf - T(0)) / static_cast<double>(nincs);
    for (int i=1; i<=nincs; i++)
    {
        // evolve
        double ti = T(0) + i*Dt;
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

template<typename Instance>
inline void ODESolver<Instance>::Evolve (pTgIncs p2Incs, Vec_t & tVec, Mat_t & Y, double tf)
{
    // check
    if (size(tVec)<2)             throw new Fatal("ODESolver::Evolve: The size of vector T must be at least equal to 2. The first item corresponds to the initial time.");
    if (Y.num_rows()<2)           throw new Fatal("ODESolver::Evolve: The number of rows in matrix Y must be at least equal to 2. The first row corresponds to the initial values.");
    if (Y.num_cols()<1)           throw new Fatal("ODESolver::Evolve: The number of columns in matrix Y must be at least equal to 1, corresponding to the number of equations.");
    if (size(tVec)!=Y.num_rows()) throw new Fatal("ODESolver::Evolve: The size of vector T must be equal to the number of rows in matrix Y");

    // data
    int nincs = size(tVec)-1;
    int neq   = Y.num_cols();

    // set pointer to function
    _p2incs = p2Incs;

    // initial values
    double   t   = tVec(0);
    double * y   = new double [neq];
    double * yFE = new double [neq];
    double * yME = new double [neq];
    double * dy1 = new double [neq];
    double * dy2 = new double [neq];
    for (int j=0; j<neq; ++j) y[j] = Y(0, j);

    // solve
    double Dt = (tf - tVec(0)) / static_cast<double>(nincs);
    for (int i=1; i<=nincs; i++)
    {
        // for each pseudo time T
        double ti  = tVec(0) + i*Dt;
        double Dti = ti - t;
        double T   = 0.0;
        double dT  = dTini;
        size_t k   = 0;
        for (k=0; k<MaxSS; ++k)
        {
            // exit point
            if (T>=1.0) break;

            // FE increments
            double dt = dT*Dti;
            (_p2inst->*_p2incs) (t, y, dt, dy1);
            for (int j=0; j<neq; ++j) yFE[j] = y[j] + dy1[j];

            // ME increments and local error
            (_p2inst->*_p2incs) (t+dt, yFE, dt, dy2);
            double norm_yME = 0.0;
            double norm_err = 0.0;
            for (int j=0; j<neq; ++j)
            {
                yME[j]    = y[j] + 0.5*(dy1[j]+dy2[j]);
                norm_yME += yME[j]*yME[j];
                norm_err += (yME[j]-yFE[j])*(yME[j]-yFE[j]);
            }

            // local error estimate
            double err = sqrt(norm_err)/(1.0+sqrt(norm_yME));

            // step multiplier
            double m = (err>0.0 ? 0.9*sqrt(STOL/err) : mMax);

            // update
            if (err<STOL)
            {
                // update state
                T += dT;
                t += dt;
                for (int j=0; j<neq; ++j) y[j] = yME[j];

                // limit change on stepsize
                if (m>mMax) m = mMax;
            }
            else if (m<mMin) m = mMin;

            // change next step size
            dT = m * dT;

            // check for last increment
            if (dT>1.0-T) dT = 1.0-T;
        }
        if (k>=MaxSS) throw new Fatal("ODESolver::Evolve: Modified-Euler did not converge after %d substeps",k);

        // results
        tVec(i) = t;
        for (int j=0; j<neq; ++j) Y(i, j) = y[j];
    }

    // clean up
    delete [] y;
    delete [] yFE;
    delete [] yME;
    delete [] dy1;
    delete [] dy2;
}

}; // namespace Numerical

#endif // MECHSYS_ODE_H
