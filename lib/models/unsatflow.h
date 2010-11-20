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
    double Findpc  (double Sw); ///< Find pc corresponding to Sw by integrating from (pc,Sw)=(0,1) to Sw (disregarding Dev)

    // Data
    double gamW;  ///< water unit weight
    double kwsat; ///< saturated (isotropic) conductivity
    double Mkw;   ///< expoent of kw model
    int    WRC;   ///< water retention curve model: 0:BC, 1:HZ, 2:ZI

    // SWRC parameters
    double bc_lam, bc_sb,  bc_wr;                      ///< Brooks & Corey model parameters
    double hz_a,   hz_b,   hz_A,   hz_B;               ///< Huang and Zienkiewicz (Gawin, Schrefler ...)
    double zi_del, zi_bet, zi_gam, zi_a, zi_b, zi_alp; ///< Zienkiewicz et al. 1990

    // Read-write variables (scratchpad)
    mutable double c, C, chi;
    mutable Mat_t  kwb;

private:
    double _kw_mult (double Sw)            const;
    double _Cpc     (double pc, double Sw) const;
    double _Ceps    (double pc, double Sw) const;

    double _RK_pc, _RK_Dev, _RK_Dpc, _RK_Sw, _RK_DSw;
    int _RK_func_Sw (double t, double const Sw[], double dSwdt[]);
    int _RK_func_pc (double t, double const pc[], double dpcdt[]);
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
    gamW    = Prms("gamW");
    kwsat   = Prms("kwsat");
    Mkw     = Prms("Mkw");
    bc_lam  = 0.8;
    bc_sb   = 1.8;
    bc_wr   = 0.01;
    hz_a    = 0.10152;
    hz_b    = 2.4279;
    hz_A    = 2.207;
    hz_B    = 1.0121;
    zi_del  = 0.0842;
    zi_bet  = 0.7; // m
    zi_gam  = 2.0;
    zi_a    = 5.0; // m
    zi_b    = 4.0;
    zi_alp  = 0.9;
    WRC     = (Prms.HasKey("WRC") ? static_cast<int>(Prms("WRC")) : 0);
    if (Prms.HasKey("bc_lam")) bc_lam = Prms("bc_lam");
    if (Prms.HasKey("bc_sb"))  bc_sb  = Prms("bc_sb");
    if (Prms.HasKey("bc_sbb")) bc_sb  = exp(Prms("bc_sbb"))-1.0;
    if (Prms.HasKey("bc_wr"))  bc_wr  = Prms("bc_wr");
    if (WRC<0 || WRC>2) throw new Fatal("UnsatFlow::UnsatFlow: WRC must be 0(BC), 1(HZ) or 2(ZI)");
}

inline void UnsatFlow::InitIvs (SDPair const & Ini, State * Sta) const
{
    SDPair ini(Ini);
    UnsatFlowState * sta = static_cast<UnsatFlowState*>(Sta);
    double pc, Sw;
    if (Ini.HasKey("pw"))
    {
        pc = -Ini("pw");
        Sw = const_cast<UnsatFlow*>(this)->FindSw(pc);
    }
    else if (Ini.HasKey("Sw"))
    {
        Sw = Ini("Sw");
        pc = const_cast<UnsatFlow*>(this)->Findpc(Sw);
    }
    else throw new Fatal("UnsatFlow::InitIvs: Either 'pw' or 'Sw' must be given in 'inis' dictionary in order to initialize WRC");
    ini.Set   ("pc",  pc);
    ini.Set   ("Sw",  Sw);
    ini.Set   ("kwb", kwsat*_kw_mult(Sw)/gamW);
    sta->Init (ini);
}

inline void UnsatFlow::Update (double Dpw, double DEv, UnsatFlowState * Sta)
{
    // set variables for _RK_func
    _RK_pc  =  Sta->pc;
    _RK_Dev =  DEv;
    _RK_Dpc = -Dpw;

    // solve ODE
    Numerical::ODESolver<UnsatFlow> ode(this, &UnsatFlow::_RK_func_Sw, /*neq*/1, "RKF45", /*stol*/1.e-6);
    ode.t    = 0.0;
    ode.Y[0] = Sta->Sw;
    ode.Evolve (/*tf*/1.0);

    // set new state
    Sta->Sw  = ode.Y[0];
    Sta->pc += (-Dpw);
    Sta->n  += DEv;
    double kwm = _kw_mult(Sta->Sw);
    for (size_t i=0; i<num_rows(Sta->kwb); ++i) Sta->kwb(i,i) = kwsat*kwm/gamW;
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
        Numerical::ODESolver<UnsatFlow> ode(this, &UnsatFlow::_RK_func_Sw, /*neq*/1, "RKF45", /*stol*/1.e-6);
        ode.t    = 0.0;
        ode.Y[0] = 1.0; // initial Sw
        ode.Evolve (/*tf*/1.0);

        // final Sw
        return ode.Y[0];
    }
    else return 1.0; // water saturated
}

inline double UnsatFlow::Findpc (double Sw)
{
    if (Sw<1.0) // water unsaturated
    {
        // set variables for _RK_func_pc
        _RK_Dev = 0.0;
        _RK_Sw  = 1.0;
        _RK_DSw = Sw-1.0;

        // solve ODE
        Numerical::ODESolver<UnsatFlow> ode(this, &UnsatFlow::_RK_func_pc, /*neq*/1, "RKF45", /*stol*/1.e-6);
        ode.t    = 0.0;
        ode.Y[0] = 0.0; // initial pc
        ode.Evolve (/*tf*/1.0);

        // final pc
        return ode.Y[0];
    }
    else return 0.0; // water saturated
}

inline double UnsatFlow::_kw_mult (double Sw) const
{
    if (Sw<1.0)
    {
        if (Sw>0.0)
        {
            if (WRC==0) // BC
            {
                return pow(Sw,Mkw);
            }
            else if (WRC==1) // HZ
            {
                double kwm = 1.0 - hz_A*pow(1.0-Sw, hz_B);
                if (kwm<0.0) return 0.0;
                else         return kwm;
            }
            else // ZI
            {
                if (Sw>zi_del)
                {
                    double hc = pow((1.0-Sw)/(Sw-zi_del), 1.0/zi_gam)/zi_bet; // m
                    return pow(1.0+pow(zi_a*hc, zi_b), -zi_alp);
                }
                else return 0.0;
            }
        }
        else return 0.0;
    }
    else return 1.0;
}

inline double UnsatFlow::_Cpc (double pc, double Sw) const
{
    if (pc>0.0 && Sw>0.0)
    {
        if (WRC==0) // BC
        {
            if (pc>bc_sb && Sw>bc_wr) return -bc_lam*pow(bc_sb/pc,bc_lam)*(1.0-bc_wr)/pc;
            else return 0.0;
        }
        else if (WRC==1) // HZ
        {
            return -(hz_a*hz_b*pow(pc/gamW,hz_b))/pc;
        }
        else // ZI
        {
            double hc = pc/gamW; // m
            double K  = pow(zi_bet*hc, zi_gam);
            return -(1.0/gamW)*(zi_gam*K*(1.0-zi_del))/(hc*pow(K+1.0,2.0));
        }
    }
    else return 0.0;
}

inline double UnsatFlow::_Ceps (double pc, double Sw) const
{
    return 0.0;
}

inline int UnsatFlow::_RK_func_Sw (double t, double const Y[], double dYdt[])
{
    double pc   = _RK_pc + t*_RK_Dpc;
    double Sw   = Y[0];
    double Cpc  = _Cpc  (pc, Sw);
    double Ceps = _Ceps (pc, Sw);
    dYdt[0] = Ceps*_RK_Dev + Cpc*_RK_Dpc; // dSwdt
    //std::cout << Y[0] << "   " << dYdt[0] << std::endl;
    return GSL_SUCCESS;
}

inline int UnsatFlow::_RK_func_pc (double t, double const Y[], double dYdt[])
{
    double Sw   = _RK_Sw + t*_RK_DSw;
    double pc   = Y[0];
    double Cpc  = _Cpc  (pc, Sw);
    double Ceps = _Ceps (pc, Sw);
    if (fabs(Cpc)>1.0e-7) dYdt[0] = (1.0/Cpc) * (_RK_DSw - Ceps*_RK_Dev); // dpcdt
    else                  dYdt[0] = 1.0e+8;
    //std::cout << "_RK_func_pc: pc = " << Y[0] << "    dpcdt = " << dYdt[0] << std::endl;
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
