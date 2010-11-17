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
#include <cmath>   // for sqrt
#include <cfloat>  // for DBL_EPSILON
#include <cstring> // for strcmp

// GSL
#include <gsl/gsl_errno.h> // for GSL_SUCCESS
#include <gsl/gsl_odeiv.h>

// MechSys
#include <mechsys/util/fatal.h>

namespace Numerical
{

template<typename Instance>
class ODESolver
{
public:
    // Typedefs
    typedef int (Instance::*pFun) (double t, double const Y[], double dYdt[]); ///< Callback function

    // Constructor
    ODESolver (Instance * p2Inst, pFun p2Fun, size_t NEq, char const * Scheme="RKDP89",
               double STOL=1.0e-6, double h=1.0e-6, double EPSREL=DBL_EPSILON);

    // Destructor
    ~ODESolver ();

    // Methods
    void Evolve (double tf); ///< Evolve from t to tf

    // Data
    double   t; ///< Current time
    double * Y; ///< Current vector of state variables

    // Internal methods
    int _call_fun (double time, double const y[], double dydt[]) { return (_p2inst->*_p2fun)(time, y, dydt); }

private:
    // Auxiliary variables
    Instance          * _p2inst;
    pFun                _p2fun;
    double              _h;
    bool                _FE;
    bool                _ME;
    gsl_odeiv_step    * _s;
    gsl_odeiv_control * _c;
    gsl_odeiv_evolve  * _e;
    gsl_odeiv_system    _sys;

    // for ME
    double * yFE;
    double * yME;
    double * dy1;
    double * dy2;
    double   stol;
};

/** Trick to pass pointers to member functions to GSL.
 * This will use (*params) to store the pointer to an instance of ODESolver, therefore (*params) will no longer be available. */
template<typename Instance>
int __ode_call_fun__ (double time, double const y[], double dydt[], void * not_used_for_params)
{
    ODESolver<Instance> * p2inst = static_cast<ODESolver<Instance>*>(not_used_for_params);
    return p2inst->_call_fun (time, y, dydt);
}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Instance>
inline ODESolver<Instance>::ODESolver (Instance * p2Inst, pFun p2Fun, size_t NEq, char const * Scheme, double STOL, double h, double EPSREL)
    : t(-1.0), Y(NULL), _p2inst(p2Inst), _p2fun(p2Fun), _h(h), _FE(false), _ME(false), yFE(NULL), yME(NULL), dy1(NULL), dy2(NULL), stol(STOL)
{
    // set scheme
    gsl_odeiv_step_type const * scheme = gsl_odeiv_step_rk8pd;
    if (Scheme!=NULL)
    {
        if      (strcmp(Scheme,"FE"    )==0) { scheme = gsl_odeiv_step_rk2; _FE=true; }
        else if (strcmp(Scheme,"RK12"  )==0) { scheme = gsl_odeiv_step_rk2; _ME=true; }
        else if (strcmp(Scheme,"RK23"  )==0) { scheme = gsl_odeiv_step_rk2;   }
        else if (strcmp(Scheme,"RKF45" )==0) { scheme = gsl_odeiv_step_rkf45; }
        else if (strcmp(Scheme,"RKDP89")==0) { scheme = gsl_odeiv_step_rk8pd; }
        else throw new Fatal("ODESolver::ODESolver: Scheme %s is invalid. Valid ones are: FE, RK12, RK23, RKF45, and RKDP89",Scheme);
    }

    // allocate auxiliary variables
    _s   = gsl_odeiv_step_alloc    (scheme, NEq);
    _c   = gsl_odeiv_control_y_new (STOL, EPSREL);
    _e   = gsl_odeiv_evolve_alloc  (NEq);
    _sys.function  = __ode_call_fun__<Instance>;
    _sys.jacobian  = NULL;
    _sys.dimension = NEq;
    _sys.params    = this;

    // allocate array
    Y = new double [NEq];

    // FE or ME
    if (_FE) dy1 = new double [NEq];
    if (_ME)
    {
        yFE = new double [NEq];
        yME = new double [NEq];
        dy1 = new double [NEq];
        dy2 = new double [NEq];
    }
}

template<typename Instance>
inline ODESolver<Instance>::~ODESolver ()
{
    if (Y  !=NULL) delete [] Y;
    if (yFE!=NULL) delete [] yFE;
    if (yME!=NULL) delete [] yME;
    if (dy1!=NULL) delete [] dy1;
    if (dy2!=NULL) delete [] dy2;
    gsl_odeiv_evolve_free  (_e);
    gsl_odeiv_control_free (_c);
    gsl_odeiv_step_free    (_s);
}

template<typename Instance>
inline void ODESolver<Instance>::Evolve (double tf)
{
    // check
    if (t<0.0) throw new Fatal("ODESolver::Evolve: Initial values (t, Y) must be set first");

    // solve
    int neq = _sys.dimension;
    if (_FE)
    {
        while (t<tf)
        {
            double dt = (t+_h>tf ? tf-t : _h);
            (_p2inst->*_p2fun) (t, Y, dy1);
            for (int j=0; j<neq; ++j) Y[j] += dy1[j]*dt;
            t += dt;
        }
    }
    else if (_ME)
    {
        // for each pseudo time T
        size_t MaxSS = 2000;
        double mMin  = 0.1;
        double mMax  = 10.0;
        double Dtf   = tf - t;
        double T     = 0.0;
        double dT    = _h;
        size_t k     = 0;
        for (k=0; k<MaxSS; ++k)
        {
            // exit point
            if (T>=1.0) break;

            // FE increments
            double dt = dT*Dtf;
            (_p2inst->*_p2fun) (t, Y, dy1);
            for (int j=0; j<neq; ++j) yFE[j] = Y[j] + dy1[j]*dt;

            // ME increments and local error
            (_p2inst->*_p2fun) (t+dt, yFE, dy2);
            double norm_yME = 0.0;
            double norm_err = 0.0;
            for (int j=0; j<neq; ++j)
            {
                yME[j]    = Y[j] + 0.5*(dy1[j]+dy2[j])*dt;
                norm_yME += yME[j]*yME[j];
                norm_err += (yME[j]-yFE[j])*(yME[j]-yFE[j]);
            }

            // local error estimate
            double err = sqrt(norm_err)/(1.0+sqrt(norm_yME));

            // step multiplier
            double m = (err>0.0 ? 0.9*sqrt(stol/err) : mMax);

            // update
            if (err<stol)
            {
                // update state
                T += dT;
                t += dt;
                for (int j=0; j<neq; ++j) Y[j] = yME[j];

                // limit change on stepsize
                if (m>mMax) m = mMax;
            }
            else if (m<mMin) m = mMin;

            // change next step size
            dT = m * dT;

            // check for last increment
            if (dT>1.0-T) dT = 1.0-T;
        }
        if (k>=MaxSS) throw new Fatal("ODESolver::Evolve: RK12 (Modified-Euler) did not converge after %d substeps",k);
    }
    else
    {
        while (t<tf)
        {
           int status = gsl_odeiv_evolve_apply (_e, _c, _s, &_sys, &t, tf, &_h, Y);
           if (status!=GSL_SUCCESS) throw new Fatal("ODESolver::Evolve: gsl_odeiv_evolve_apply failed. status=%d",status);
           //if (h<hmin) { h = hmin; }
        }
    }
}

}; // namespace Numerical

#endif // MECHSYS_ODE_H
