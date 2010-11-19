/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
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

#ifndef MECHSYS_UNSATFLOW_H
#define MECHSYS_UNSATFLOW_H

// MechSys
#include<mechsys/models/model.h>
#include<mechsys/linalg/matvec.h>
#include<mechsys/numerical/odesolver.h>

class UnsatFlowState : public State
{
public:
    // static
    static Mat_t Im; // identity 3x3

    // Constructor
    UnsatFlowState (int NDim);

    // Methods
    void Init    (SDPair const & Ini, size_t NIvs=0);
    void Backup  () {}
    void Restore () {}

    // Data
    double n, Sw; // porosity, saturation
    double pc;    // capillary pressure (pc=pa-pw=-pw)
    Mat_t  kwb;   // current conductivity divided by gammaW
};

class UnsatFlow : public Model
{
public:
    // Constructor
    UnsatFlow (int NDim, SDPair const & Prms);

    // Methods
    void   InitIvs (SDPair const & Ini, State * Sta) const;
    void   Update  (double Dpw, double DEv, UnsatFlowState * Sta);
    void   TgVars  (UnsatFlowState const * Sta) const;
    double FindSw  (double pc); ///< Find Sw corresponding to pc by integrating from (pc,Sw)=(0,1) to pc (disregarding Dev)

    // Data
    double gamW;  // water unit weight
    double kwsat; // saturated (isotropic) conductivity
    double Mkw;   // expoent of kw model

    // SWRC parameters
    double bc_lam, bc_sb, bc_wr; ///< Brooks & Corey model parameters

    // Read-write variables (scratchpad)
    mutable double c, C, chi;
    mutable Mat_t  kwb;

private:
    double _Cpc  (double pc, double Sw) const;
    double _Ceps (double pc, double Sw) const;

    double _RK_pc, _RK_Dev, _RK_Dpc;
    int _RK_func (double t, double const Sw[], double dSwdt[]);
};

Mat_t UnsatFlowState::Im;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline UnsatFlowState::UnsatFlowState (int NDim)
    : State(NDim), n(0.2), Sw(1.0), pc(0.0)
{
    if (num_rows(Im)!=(size_t)NDim)
    {
        Im.change_dim (NDim, NDim);
        if (NDim==3)
        {
            Im = 1., 0., 0.,
                 0., 1., 0.,
                 0., 0., 1.;
        }
        else
        {
            Im = 1., 0.,
                 0., 1.;
        }
    }
    double kk = 1.0e-2/9.81;
    kwb = kk*Im;
}

inline void UnsatFlowState::Init (SDPair const & Ini, size_t NIvs)
{
    n   = Ini("n");
    pc  = Ini("pc");
    Sw  = Ini("Sw");
    kwb = Ini("kwb")*Im;
}

inline UnsatFlow::UnsatFlow (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"UnsatFlow")
{
    gamW   = Prms("gamW");
    kwsat  = Prms("kwsat");
    Mkw    = Prms("Mkw");
    bc_lam = 0.8;
    bc_sb  = 1.8;
    bc_wr  = 0.01;
    if (Prms.HasKey("bc_lam")) bc_lam = Prms("bc_lam");
    if (Prms.HasKey("bc_sb"))  bc_sb  = Prms("bc_sb");
    if (Prms.HasKey("bc_sbb")) bc_sb  = exp(Prms("bc_sbb"))-1.0;
    if (Prms.HasKey("bc_wr"))  bc_wr  = Prms("bc_wr");
}

inline void UnsatFlow::InitIvs (SDPair const & Ini, State * Sta) const
{
    SDPair ini(Ini);
    UnsatFlowState * sta = static_cast<UnsatFlowState*>(Sta);
    double pc = -Ini("pw");
    double Sw =  Ini("Sw");
    double f;
    if (pc>bc_sb) f = Sw - bc_wr - (1.0-bc_wr)*pow(bc_sb/pc,bc_lam);
    else          f = Sw - 1.0;
    if (fabs(f)>1.0e-7) throw new Fatal("UnsatFlow::InitIvs: Sw=%g and pc=%g are not in WR curve",Sw,pc);
    ini.Set   ("pc",  pc);
    ini.Set   ("kwb", kwsat*pow(Sw,Mkw)/gamW);
    sta->Init (ini);
}

inline void UnsatFlow::Update (double Dpw, double DEv, UnsatFlowState * Sta)
{
    // set variables for _RK_func
    _RK_pc  =  Sta->pc;
    _RK_Dev =  DEv;
    _RK_Dpc = -Dpw;

    // solve ODE
    Numerical::ODESolver<UnsatFlow> ode(this, &UnsatFlow::_RK_func, /*neq*/1, "RKF45", /*stol*/1.e-6);
    ode.t    = 0.0;
    ode.Y[0] = Sta->Sw;
    ode.Evolve (/*tf*/1.0);

    // set new state
    Sta->Sw  = ode.Y[0];
    Sta->pc += (-Dpw);
    Sta->n  += DEv;
    for (size_t i=0; i<num_rows(Sta->kwb); ++i) Sta->kwb(i,i) = kwsat*pow(Sta->Sw,Mkw)/gamW;
}

inline void UnsatFlow::TgVars (UnsatFlowState const * Sta) const
{
    double Cpc  = _Cpc  (Sta->pc, Sta->Sw);
    double Ceps = _Ceps (Sta->pc, Sta->Sw);
    c   = Sta->Sw + Sta->n * Ceps;
    C   = -Sta->n * Cpc;
    chi = Sta->Sw;
    kwb = Sta->kwb;
}

inline double UnsatFlow::FindSw (double pc)
{
    if (pc>0.0) // water unsaturated
    {
        // set variables for _RK_func
        _RK_Dev = 0.0;
        _RK_pc  = 0.0;
        _RK_Dpc = pc;

        // solve ODE
        Numerical::ODESolver<UnsatFlow> ode(this, &UnsatFlow::_RK_func, /*neq*/1, "RKF45", /*stol*/1.e-6);
        ode.t    = 0.0;
        ode.Y[0] = 1.0; // initial Sw
        ode.Evolve (/*tf*/1.0);

        // final Sw
        return ode.Y[0];
    }
    else return 1.0; // water saturated
}

inline double UnsatFlow::_Cpc (double pc, double Sw) const
{
    if (Sw<=bc_wr) return 0.0;
    if (pc<bc_sb || pc<1.0e-7) return 0.0;
    else return -bc_lam*pow(bc_sb/pc,bc_lam)*(1.0-bc_wr)/pc;
}

inline double UnsatFlow::_Ceps (double pc, double Sw) const
{
    return 0.0;
}

inline int UnsatFlow::_RK_func (double t, double const Y[], double dYdt[])
{
    double pc   = _RK_pc + t*_RK_Dpc;
    double Sw   = Y[0];
    double Cpc  = _Cpc  (pc, Sw);
    double Ceps = _Ceps (pc, Sw);
    dYdt[0] = Ceps*_RK_Dev + Cpc*_RK_Dpc;
    return GSL_SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * UnsatFlowMaker(int NDim, SDPair const & Prms) { return new UnsatFlow(NDim,Prms); }

int UnsatFlowRegister()
{
    ModelFactory["UnsatFlow"] = UnsatFlowMaker;
    MODEL.Set ("UnsatFlow", (double)MODEL.Keys.Size());
    return 0;
}

int __UnsatFlow_dummy_int = UnsatFlowRegister();


#endif // MECHSYS_UNSATFLOW_H
