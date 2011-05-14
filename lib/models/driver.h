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

#ifndef MECHSYS_DRIVER_H
#define MECHSYS_DRIVER_H

// Std Lib
#include <cmath>
#include <sstream>

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/numerical/odesolver.h>

using std::cout;
using std::endl;

class Driver
{
public:
    // Constructor & Destructor
     Driver (String const & ModelName, SDPair const & Prms, SDPair const & Inis);
    ~Driver ();

    // Methods
    void Update ();

    // Data to be set before calling update methods
    Array<bool> PrescDEps;     ///< Prescribed strains
    Vec_t       DEps;          ///< Strain increments
    Vec_t       DSig;          ///< Stress increments

    // Data
    Model        * Mdl;              ///< Model
    EquilibState * Sta;              ///< State
    double         Tol;              ///< Tolerance to check for division by zero
    bool           CheckModelTgIncs; ///< Check model TgIncs ?
    String         RKScheme;         ///< Scheme
    double         RKStol;           ///< STOL
    Vec_t          DIvs;             ///< Increments of internal values
    int            GivenIncsCase;    ///< 3 components cases: 8 cases
    bool           Linear;           ///< Linear model ?

private:
    // Callbacks
    int _RK_3Comps (double t, double const Y[], double dYdt[]);
    EquilibState * _RK_state;
    Vec_t          _RK_EpsIni;
    Vec_t          _RK_SigIni;
    double         _tmp_dYdt[3];
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation


inline Driver::Driver (String const & ModelName, SDPair const & Prms, SDPair const & Inis)
    : Tol(DBL_EPSILON), CheckModelTgIncs(true), RKScheme("RKF45"), RKStol(1.0e-5), GivenIncsCase(-1),
      Linear(false)
{
    Mdl = AllocModel (ModelName, 3, Prms, /*AnotherMdl*/NULL);
    Sta = new EquilibState (3);
    Mdl->InitIvs (Inis, Sta);

    PrescDEps.Resize (6);
    DEps.change_dim  (6);
    DSig.change_dim  (6);
    if (Mdl->NIvs>0) DIvs.change_dim (Mdl->NIvs);

    _RK_state = new EquilibState (3);
    _RK_EpsIni.change_dim (6);
    _RK_SigIni.change_dim (6);
}

inline Driver::~Driver ()
{
    delete Mdl;
    delete Sta;
    delete _RK_state;
}

inline void Driver::Update ()
{
    // ODE solver
    Numerical::ODESolver<Driver> ode(this, &Driver::_RK_3Comps, /*neq*/3+Mdl->NIvs, RKScheme.CStr(), RKStol);

    // initial values
    ode.t = 0.0;
    if ( PrescDEps[0] &&  PrescDEps[1] &&  PrescDEps[2]) // 1) all DEps are given
    {
        GivenIncsCase = 1;
        ode.Y[0] = Sta->Sig(0);
        ode.Y[1] = Sta->Sig(1);
        ode.Y[2] = Sta->Sig(2);
    }
    else if (!PrescDEps[0] && !PrescDEps[1] && !PrescDEps[2]) // 2) all DSig are given
    {
        GivenIncsCase = 2;
        ode.Y[0] = Sta->Eps(0);
        ode.Y[1] = Sta->Eps(1);
        ode.Y[2] = Sta->Eps(2);
    }
    else if ( PrescDEps[0] && !PrescDEps[1] && !PrescDEps[2]) // 3) DEps[0], DSig[1], and DSig[2] are given
    {
        GivenIncsCase = 3;
        ode.Y[0] = Sta->Sig(0);
        ode.Y[1] = Sta->Eps(1);
        ode.Y[2] = Sta->Eps(2);
    }
    else if (!PrescDEps[0] &&  PrescDEps[1] && !PrescDEps[2]) // 4) DSig[0], DEps[1], and DSig[2] are given
    {
        GivenIncsCase = 4;
        ode.Y[0] = Sta->Eps(0);
        ode.Y[1] = Sta->Sig(1);
        ode.Y[2] = Sta->Eps(2);
    }
    else if (!PrescDEps[0] && !PrescDEps[1] &&  PrescDEps[2]) // 5) DSig[0], DSig[1], and DEps[2] are given
    {
        GivenIncsCase = 5;
        ode.Y[0] = Sta->Eps(0);
        ode.Y[1] = Sta->Eps(1);
        ode.Y[2] = Sta->Sig(2);
    }
    else if ( PrescDEps[0] &&  PrescDEps[1] && !PrescDEps[2]) // 6) DEps[0], DEps[1], and DSig[2] are given
    {
        GivenIncsCase = 6;
        ode.Y[0] = Sta->Sig(0);
        ode.Y[1] = Sta->Sig(1);
        ode.Y[2] = Sta->Eps(2);
    }
    else if ( PrescDEps[0] && !PrescDEps[1] &&  PrescDEps[2]) // 7) DEps[0], DSig[1], and DEps[2] are given
    {
        GivenIncsCase = 7;
        ode.Y[0] = Sta->Sig(0);
        ode.Y[1] = Sta->Eps(1);
        ode.Y[2] = Sta->Sig(2);
    }
    else if (!PrescDEps[0] &&  PrescDEps[1] &&  PrescDEps[2]) // 8) DSig[0], DEps[1], and DEps[2] are given
    {
        GivenIncsCase = 8;
        ode.Y[0] = Sta->Eps(0);
        ode.Y[1] = Sta->Sig(1);
        ode.Y[2] = Sta->Sig(2);
    }
    else throw new Fatal("Driver::Update: __internal_error__: combination of prescribed DEps is invalid");

    // initial internal values
    for (size_t i=0; i<Mdl->NIvs; ++i) ode.Y[3+i] = Sta->Ivs(i);

    // set variables for _RK_3Comps
    _RK_EpsIni        = Sta->Eps;
    _RK_SigIni        = Sta->Sig;
    _RK_state->Eps(3) = 0.0;
    _RK_state->Eps(4) = 0.0;
    _RK_state->Eps(5) = 0.0;
    _RK_state->Sig(3) = 0.0;
    _RK_state->Sig(4) = 0.0;
    _RK_state->Sig(5) = 0.0;
    if (Mdl->NIvs>0) _RK_state->Ivs.change_dim (Mdl->NIvs);

    // loading condition
    _RK_3Comps (0.0, ode.Y, _tmp_dYdt);
    if (GivenIncsCase==1) // 1) all DEps are given
    {
        DSig(0) = _tmp_dYdt[0]; // * 1.0
        DSig(1) = _tmp_dYdt[1];
        DSig(2) = _tmp_dYdt[2];
    }
    else if (GivenIncsCase==2) // 2) all DSig are given
    {
        DEps(0) = _tmp_dYdt[0];
        DEps(1) = _tmp_dYdt[1];
        DEps(2) = _tmp_dYdt[2];
    }
    else if (GivenIncsCase==3) // 3) DEps[0], DSig[1], and DSig[2] are given
    {
        DSig(0) = _tmp_dYdt[0];
        DEps(1) = _tmp_dYdt[1];
        DEps(2) = _tmp_dYdt[2];
    }
    else if (GivenIncsCase==4) // 4) DSig[0], DEps[1], and DSig[2] are given
    {
        DEps(0) = _tmp_dYdt[0];
        DSig(1) = _tmp_dYdt[1];
        DEps(2) = _tmp_dYdt[2];
    }
    else if (GivenIncsCase==5) // 5) DSig[0], DSig[1], and DEps[2] are given
    {
        DEps(0) = _tmp_dYdt[0];
        DEps(1) = _tmp_dYdt[1];
        DSig(2) = _tmp_dYdt[2];
    }
    else if (GivenIncsCase==6) // 6) DEps[0], DEps[1], and DSig[2] are given
    {
        DSig(0) = _tmp_dYdt[0];
        DSig(1) = _tmp_dYdt[1];
        DEps(2) = _tmp_dYdt[2];
    }
    else if (GivenIncsCase==7) // 7) DEps[0], DSig[1], and DEps[2] are given
    {
        DSig(0) = _tmp_dYdt[0];
        DEps(1) = _tmp_dYdt[1];
        DSig(2) = _tmp_dYdt[2];
    }
    else if (GivenIncsCase==6) // 8) DSig[0], DEps[1], and DEps[2] are given
    {
        DEps(0) = _tmp_dYdt[0];
        DSig(1) = _tmp_dYdt[1];
        DSig(2) = _tmp_dYdt[2];
    }
    double aint = -1.0; // no intersection
    Sta->Ldg  = Mdl->LoadCond (Sta, DEps, aint); // returns true if there is loading (also when there is intersection)
    _RK_state->Ldg = Sta->Ldg;


    cout << "case = " << GivenIncsCase << endl;


    // solve
    if (Linear)
    {
        Sta->Sig += DSig;
        Sta->Eps += DEps;
        return;
    }
    else ode.Evolve (/*tf*/1.0);

    // results
    if (GivenIncsCase==1) // 1) all DEps are given
    {
        Sta->Sig(0) = ode.Y[0];
        Sta->Sig(1) = ode.Y[1];
        Sta->Sig(2) = ode.Y[2];
    }
    else if (GivenIncsCase==2) // 2) all DSig are given
    {
        Sta->Eps(0) = ode.Y[0];
        Sta->Eps(1) = ode.Y[1];
        Sta->Eps(2) = ode.Y[2];
    }
    else if (GivenIncsCase==3) // 3) DEps[0], DSig[1], and DSig[2] are given
    {
        Sta->Sig(0) = ode.Y[0];
        Sta->Eps(1) = ode.Y[1];
        Sta->Eps(2) = ode.Y[2];
    }
    else if (GivenIncsCase==4) // 4) DSig[0], DEps[1], and DSig[2] are given
    {
        Sta->Eps(0) = ode.Y[0];
        Sta->Sig(1) = ode.Y[1];
        Sta->Eps(2) = ode.Y[2];
    }
    else if (GivenIncsCase==5) // 5) DSig[0], DSig[1], and DEps[2] are given
    {
        Sta->Eps(0) = ode.Y[0];
        Sta->Eps(1) = ode.Y[1];
        Sta->Sig(2) = ode.Y[2];
    }
    else if (GivenIncsCase==6) // 6) DEps[0], DEps[1], and DSig[2] are given
    {
        Sta->Sig(0) = ode.Y[0];
        Sta->Sig(1) = ode.Y[1];
        Sta->Eps(2) = ode.Y[2];
    }
    else if (GivenIncsCase==7) // 7) DEps[0], DSig[1], and DEps[2] are given
    {
        Sta->Sig(0) = ode.Y[0];
        Sta->Eps(1) = ode.Y[1];
        Sta->Sig(2) = ode.Y[2];
    }
    else if (GivenIncsCase==6) // 8) DSig[0], DEps[1], and DEps[2] are given
    {
        Sta->Eps(0) = ode.Y[0];
        Sta->Sig(1) = ode.Y[1];
        Sta->Sig(2) = ode.Y[2];
    }

    // internal values at time t
    for (size_t i=0; i<Mdl->NIvs; ++i) Sta->Ivs(i) = ode.Y[3+i];
}

inline int Driver::_RK_3Comps (double t, double const Y[], double dYdt[])
{
    // state at time t
    if (GivenIncsCase==1) // 1) all DEps are given
    {
        _RK_state->Eps(0) = _RK_EpsIni(0) + t * DEps(0);
        _RK_state->Eps(1) = _RK_EpsIni(1) + t * DEps(1);
        _RK_state->Eps(2) = _RK_EpsIni(2) + t * DEps(2);
        _RK_state->Sig(0) = Y[0];
        _RK_state->Sig(1) = Y[1];
        _RK_state->Sig(2) = Y[2];
    }
    else if (GivenIncsCase==2) // 2) all DSig are given
    {
        _RK_state->Eps(0) = Y[0];
        _RK_state->Eps(1) = Y[1];
        _RK_state->Eps(2) = Y[2];
        _RK_state->Sig(0) = _RK_SigIni(0) + t * DSig(0);
        _RK_state->Sig(1) = _RK_SigIni(1) + t * DSig(1);
        _RK_state->Sig(2) = _RK_SigIni(2) + t * DSig(2);
    }
    else if (GivenIncsCase==3) // 3) DEps[0], DSig[1], and DSig[2] are given
    {
        _RK_state->Eps(0) = _RK_EpsIni(0) + t * DEps(0);
        _RK_state->Eps(1) = Y[1];
        _RK_state->Eps(2) = Y[2];
        _RK_state->Sig(0) = Y[0];
        _RK_state->Sig(1) = _RK_SigIni(1) + t * DSig(1);
        _RK_state->Sig(2) = _RK_SigIni(2) + t * DSig(2);
    }
    else if (GivenIncsCase==4) // 4) DSig[0], DEps[1], and DSig[2] are given
    {
        _RK_state->Eps(0) = Y[0];
        _RK_state->Eps(1) = _RK_EpsIni(1) + t * DEps(1);
        _RK_state->Eps(2) = Y[2];
        _RK_state->Sig(0) = _RK_SigIni(0) + t * DSig(0);
        _RK_state->Sig(1) = Y[1];
        _RK_state->Sig(2) = _RK_SigIni(2) + t * DSig(2);
    }
    else if (GivenIncsCase==5) // 5) DSig[0], DSig[1], and DEps[2] are given
    {
        _RK_state->Eps(0) = Y[0];
        _RK_state->Eps(1) = Y[1];
        _RK_state->Eps(2) = _RK_EpsIni(2) + t * DEps(2);
        _RK_state->Sig(0) = _RK_SigIni(0) + t * DSig(0);
        _RK_state->Sig(1) = _RK_SigIni(1) + t * DSig(1);
        _RK_state->Sig(2) = Y[2];
    }
    else if (GivenIncsCase==6) // 6) DEps[0], DEps[1], and DSig[2] are given
    {
        _RK_state->Eps(0) = _RK_EpsIni(0) + t * DEps(0);
        _RK_state->Eps(1) = _RK_EpsIni(1) + t * DEps(1);
        _RK_state->Eps(2) = Y[2];
        _RK_state->Sig(0) = Y[0];
        _RK_state->Sig(1) = Y[1];
        _RK_state->Sig(2) = _RK_SigIni(2) + t * DSig(2);
    }
    else if (GivenIncsCase==7) // 7) DEps[0], DSig[1], and DEps[2] are given
    {
        _RK_state->Eps(0) = _RK_EpsIni(0) + t * DEps(0);
        _RK_state->Eps(1) = Y[1];
        _RK_state->Eps(2) = _RK_EpsIni(2) + t * DEps(2);
        _RK_state->Sig(0) = Y[0];
        _RK_state->Sig(1) = _RK_SigIni(1) + t * DSig(1);
        _RK_state->Sig(2) = Y[2];
    }
    else if (GivenIncsCase==6) // 8) DSig[0], DEps[1], and DEps[2] are given
    {
        _RK_state->Eps(0) = Y[0];
        _RK_state->Eps(1) = _RK_EpsIni(1) + t * DEps(1);
        _RK_state->Eps(2) = _RK_EpsIni(2) + t * DEps(2);
        _RK_state->Sig(0) = _RK_SigIni(0) + t * DSig(0);
        _RK_state->Sig(1) = Y[1];
        _RK_state->Sig(2) = Y[2];
    }

    // internal values at time t
    for (size_t i=0; i<Mdl->NIvs; ++i) _RK_state->Ivs(i) = Y[3+i];

    // stiffness
    Mat_t D;
    Mdl->Stiffness (_RK_state, D);

    // stress/strain increments
    if (GivenIncsCase==1) // 1) all DEps are given
    {
        dYdt[0] = DSig(0) = D(0,0)*DEps(0) + D(0,1)*DEps(1) + D(0,2)*DEps(2);
        dYdt[1] = DSig(1) = D(1,0)*DEps(0) + D(1,1)*DEps(1) + D(1,2)*DEps(2);
        dYdt[2] = DSig(2) = D(2,0)*DEps(0) + D(2,1)*DEps(1) + D(2,2)*DEps(2);
    }
    else if (GivenIncsCase==2) // 2) all DSig are given
    {
        double den2 = D(0,0)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+D(0,1)*(D(1,0)*D(2,2)-D(1,2)*D(2,0))+D(0,2)*(D(1,1)*D(2,0)-D(1,0)*D(2,1));
        if (fabs(den2)<Tol) throw new Fatal("Driver::TgIncs3Comps: division by null denominator (%g)",den2);
        dYdt[0] = DEps(0) =  (DSig(0)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+D(0,1)*(DSig(1)*D(2,2)-D(1,2)*DSig(2))+D(0,2)*(D(1,1)*DSig(2)-DSig(1)*D(2,1)))/den2;
        dYdt[1] = DEps(1) = -(DSig(0)*(D(1,2)*D(2,0)-D(1,0)*D(2,2))+D(0,0)*(DSig(1)*D(2,2)-D(1,2)*DSig(2))+D(0,2)*(D(1,0)*DSig(2)-DSig(1)*D(2,0)))/den2;
        dYdt[2] = DEps(2) =  (DSig(0)*(D(1,1)*D(2,0)-D(1,0)*D(2,1))+D(0,0)*(DSig(1)*D(2,1)-D(1,1)*DSig(2))+D(0,1)*(D(1,0)*DSig(2)-DSig(1)*D(2,0)))/den2;
    }
    else if (GivenIncsCase==3) // 3) DEps[0], DSig[1], and DSig[2] are given
    {
        double den3 = D(1,2)*D(2,1)-D(1,1)*D(2,2);
        if (fabs(den3)<Tol) throw new Fatal("Driver::TgIncs3Comps: division by null denominator (%g)",den3);
        dYdt[0] = DSig(0) =  (DEps(0)*D(0,0)*den3+D(0,1)*(DEps(0)*(D(1,0)*D(2,2)-D(1,2)*D(2,0))-DSig(1)*D(2,2)+D(1,2)*DSig(2))+D(0,2)*(DEps(0)*(D(1,1)*D(2,0)-D(1,0)*D(2,1))+DSig(1)*D(2,1)-D(1,1)*DSig(2)))/den3;
        dYdt[1] = DEps(1) = -(DEps(0)*(D(1,2)*D(2,0)-D(1,0)*D(2,2))+DSig(1)*D(2,2)-D(1,2)*DSig(2))/den3;
        dYdt[2] = DEps(2) =  (DEps(0)*(D(1,1)*D(2,0)-D(1,0)*D(2,1))+DSig(1)*D(2,1)-D(1,1)*DSig(2))/den3;
    }
    else if (GivenIncsCase==4) // 4) DSig[0], DEps[1], and DSig[2] are given
    {
        double den4 = D(0,2)*D(2,0)-D(0,0)*D(2,2);
        if (fabs(den4)<Tol) throw new Fatal("Driver::TgIncs3Comps: division by null denominator (%g)",den4);
        dYdt[0] = DEps(0) = -(DEps(1)*(D(0,2)*D(2,1)-D(0,1)*D(2,2))+DSig(0)*D(2,2)-D(0,2)*DSig(2))/den4;
        dYdt[1] = DSig(1) =  (DEps(1)*(D(0,0)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+D(0,1)*(D(1,0)*D(2,2)-D(1,2)*D(2,0))+D(0,2)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))+DSig(0)*(D(1,2)*D(2,0)-D(1,0)*D(2,2))-D(0,0)*D(1,2)*DSig(2)+D(0,2)*D(1,0)*DSig(2))/den4;
        dYdt[2] = DEps(2) = -(DEps(1)*(D(0,1)*D(2,0)-D(0,0)*D(2,1))-DSig(0)*D(2,0)+D(0,0)*DSig(2))/den4;
    }
    else if (GivenIncsCase==5) // 5) DSig[0], DSig[1], and DEps[2] are given
    {
        double den5 = D(0,1)*D(1,0)-D(0,0)*D(1,1);
        if (fabs(den5)<Tol) throw new Fatal("Driver::TgIncs3Comps: division by null denominator (%g)",den5);
        dYdt[0] = DEps(0) =  ((D(0,2)*D(1,1)-D(0,1)*D(1,2))*DEps(2)-DSig(0)*D(1,1)+D(0,1)*DSig(1))/den5;
        dYdt[1] = DEps(1) = -((D(0,2)*D(1,0)-D(0,0)*D(1,2))*DEps(2)-DSig(0)*D(1,0)+D(0,0)*DSig(1))/den5;
        dYdt[2] = DSig(2) =  (DEps(2)*(D(0,0)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+D(0,1)*(D(1,0)*D(2,2)-D(1,2)*D(2,0))+D(0,2)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))+DSig(0)*(D(1,0)*D(2,1)-D(1,1)*D(2,0))-D(0,0)*DSig(1)*D(2,1)+D(0,1)*DSig(1)*D(2,0))/den5;
    }
    else if (GivenIncsCase==6) // 6) DEps[0], DEps[1], and DSig[2] are given
    {
        if (fabs(D(2,2))<Tol) throw new Fatal("Driver::TgIncs3Comps: division by null denominator (%g)",D(2,2));
        dYdt[0] = DSig(0) = -(DEps(1)*(D(0,2)*D(2,1)-D(0,1)*D(2,2))-DEps(0)*D(0,0)*D(2,2)+D(0,2)*(DEps(0)*D(2,0)-DSig(2)))/D(2,2);
        dYdt[1] = DSig(1) = -(DEps(1)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+DEps(0)*(D(1,2)*D(2,0)-D(1,0)*D(2,2))-D(1,2)*DSig(2))/D(2,2);
        dYdt[2] = DEps(2) = -(DEps(1)*D(2,1)+DEps(0)*D(2,0)-DSig(2))/D(2,2);
    }
    else if (GivenIncsCase==7) // 7) DEps[0], DSig[1], and DEps[2] are given
    {
        if (fabs(D(1,1))<Tol) throw new Fatal("Driver::TgIncs3Comps: division by null denominator (%g)",D(1,1));
        dYdt[0] = DSig(0) = ((D(0,2)*D(1,1)-D(0,1)*D(1,2))*DEps(2)+DEps(0)*D(0,0)*D(1,1)+D(0,1)*(DSig(1)-DEps(0)*D(1,0)))/D(1,1);
        dYdt[1] = DEps(1) = -(D(1,2)*DEps(2)+DEps(0)*D(1,0)-DSig(1))/D(1,1);
        dYdt[2] = DSig(2) = -(DEps(2)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+DEps(0)*(D(1,0)*D(2,1)-D(1,1)*D(2,0))-DSig(1)*D(2,1))/D(1,1);
    }
    else if (GivenIncsCase==8) // 8) DSig[0], DEps[1], and DEps[2] are given
    {
        if (fabs(D(0,0))<Tol) throw new Fatal("Driver::TgIncs3Comps: division by null denominator (%g)",D(0,0));
        dYdt[0] = DEps(0) = -(D(0,2)*DEps(2)+D(0,1)*DEps(1)-DSig(0))/D(0,0);
        dYdt[1] = DSig(1) = -((D(0,2)*D(1,0)-D(0,0)*D(1,2))*DEps(2)+DEps(1)*(D(0,1)*D(1,0)-D(0,0)*D(1,1))-DSig(0)*D(1,0))/D(0,0);
        dYdt[2] = DSig(2) = -(DEps(2)*(D(0,2)*D(2,0)-D(0,0)*D(2,2))+DEps(1)*(D(0,1)*D(2,0)-D(0,0)*D(2,1))-DSig(0)*D(2,0))/D(0,0);
    }

    // increments of internal values
    if (Mdl->NIvs>0)
    {
        Vec_t dsig_tmp;
        Mdl->TgIncs (_RK_state, DEps, dsig_tmp, DIvs);
        for (size_t i=0; i<Mdl->NIvs; ++i) dYdt[3+i] = DIvs(i);
        if (CheckModelTgIncs)
        {
            Vec_t err(DSig - dsig_tmp);
            if (Norm(err)>1.0e-8)
            {
                std::ostringstream oss;
                oss << "  dsig_tmp  = " << PrintVector(dsig_tmp);
                oss << "  DSig      = " << PrintVector(DSig);
                oss << "  Norm(err) = " << Norm(err) << std::endl;
                throw new Fatal("Driver::_RK_3Comps: __internal_error__: DSig differs from dsig_tmp calculated by model\n%s",oss.str().c_str());
            }
        }
    }

    return GSL_SUCCESS;
}

#endif
