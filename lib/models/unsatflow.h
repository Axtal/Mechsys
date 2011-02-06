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

// Std Lib
#include <fstream>
#include <sstream>

// MechSys
#include<mechsys/util/numstreams.h>
#include<mechsys/models/model.h>
#include<mechsys/linalg/matvec.h>
#include<mechsys/numerical/odesolver.h>

class UnsatFlowState : public State
{
public:
    // Constructor
    UnsatFlowState (int NDim) : State(NDim) {}

    // Methods
    void   Init    (SDPair const & Ini, size_t NIvs=0);
    void   Backup  () { throw new Fatal("UnsatFlowState::Backup: this method is not available yet"); }
    void   Restore () { throw new Fatal("UnsatFlowState::Restore: this method is not available yet"); }
    size_t PckSize () const { return 3; }
    void   Pack    (Array<double>       & V) const;
    void   Unpack  (Array<double> const & V);

    // Data
    double n;  ///< Porosity
    double Sw; ///< Saturation
    double pc; ///< Capillary pressure (pc = pa-pw = -pw)
};

class UnsatFlow : public Model
{
public:
    // enums
    enum WRC_t { BC_t, ZI_t }; ///< WRC type

    // Constructor
    UnsatFlow (int NDim, SDPair const & Prms);

    // Methods
    void   InitIvs  (SDPair const & Ini, State * Sta) const;                  ///< Initialize internal values
    void   Update   (double Dpw, double DEv, UnsatFlowState * Sta);           ///< Update state
    void   TgVars   (UnsatFlowState const * Sta) const;                       ///< Calculate c, C, chi, and kwb
    double FindSw   (double pc);                                              ///< Find Sw corresponding to pc by integrating from (pc,Sw)=(0,1) to pc (disregarding Dev)
    double Findpc   (double Sw);                                              ///< Find pc corresponding to Sw by integrating from (pc,Sw)=(0,1) to Sw (disregarding Dev)
    double mkw      (double Sw) const;                                        ///< kw multiplier
    void   GenCurve (Array<double> pcs, char const * Filekey, size_t Np=100); ///< Generate WRC

    // Data
    Mat_t  kwb_sat; ///< Saturated kw bar
    double kwsat;   ///< Saturated (isotropic) conductivity
    double akw;     ///< Expoent of kw model
    WRC_t  wrc;     ///< Water retention curve model

    // SWRC parameters
    double bc_lam, bc_sb,  bc_wr;                      ///< Brooks & Corey model parameters
    double zi_del, zi_bet, zi_gam, zi_a, zi_b, zi_alp; ///< Zienkiewicz et al. 1990 (zi_bet and zi_a => [m])

    // Read-write variables (scratchpad)
    mutable double c, C, chi;
    mutable Mat_t  kwb;

private:
    double _Cpc  (double pc, double Sw) const;
    double _Ceps (double pc, double Sw) const;

    double _RK_pc, _RK_Dev, _RK_Dpc, _RK_Sw, _RK_DSw;
    int _RK_func_Sw (double t, double const Sw[], double dSwdt[]);
    int _RK_func_pc (double t, double const pc[], double dpcdt[]);
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void UnsatFlowState::Init (SDPair const & Ini, size_t NIvs)
{
    n  = Ini("n");
    pc = Ini("pc");
    Sw = Ini("Sw");
}

inline void UnsatFlowState::Pack (Array<double> & V) const
{
    V.Resize (PckSize());
    V[0] = n;
    V[1] = Sw;
    V[2] = pc;
}

inline void UnsatFlowState::Unpack (Array<double> const & V)
{
    if (V.Size()!=PckSize()) throw new Fatal("UnsatFlowState::Unpack: Size of given vector (%zd) is different of correct size of Pack (%zd)",V.Size(),PckSize());
    n  = V[0];
    Sw = V[1];
    pc = V[2];
}

inline UnsatFlow::UnsatFlow (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,/*niv*/0,"UnsatFlow")
{
    // parameters
    kwsat   = Prms("kwsat");
    akw     = Prms("akw");
    bc_lam  = (Prms.HasKey("bc_lam") ? Prms("bc_lam") : 0.8    );
    bc_sb   = (Prms.HasKey("bc_sb" ) ? Prms("bc_sb" ) : 1.8    );
    bc_wr   = (Prms.HasKey("bc_wr" ) ? Prms("bc_wr" ) : 0.01   );
    zi_del  = (Prms.HasKey("zi_del") ? Prms("zi_del") : 0.0842 );
    zi_bet  = (Prms.HasKey("zi_bet") ? Prms("zi_bet") : 0.7    );
    zi_gam  = (Prms.HasKey("zi_gam") ? Prms("zi_gam") : 2.0    );
    zi_a    = (Prms.HasKey("zi_a"  ) ? Prms("zi_a"  ) : 5.0    );
    zi_b    = (Prms.HasKey("zi_b"  ) ? Prms("zi_b"  ) : 4.0    );
    zi_alp  = (Prms.HasKey("zi_alp") ? Prms("zi_alp") : 0.9    );
    wrc     = BC_t;
    if (Prms.HasKey("bc_sbb")) bc_sb = exp(Prms("bc_sbb"))-1.0;
    if (Prms.HasKey("ZI"))     wrc   = ZI_t;

    // saturated conductivity matrix
    kwb    .change_dim (NDim, NDim);
    kwb_sat.change_dim (NDim, NDim);
    double m = kwsat / Prms("gamW");
    if (NDim==3) kwb_sat =  m, 0., 0.,   0., m, 0.,   0., 0., m;
    else         kwb_sat =  m, 0.,   0., m;
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
    ini.Set   ("n",   Prms("por"));
    ini.Set   ("pc",  pc);
    ini.Set   ("Sw",  Sw);
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
}

inline void UnsatFlow::TgVars (UnsatFlowState const * Sta) const
{
    double Cpc  = _Cpc  (Sta->pc, Sta->Sw);
    double Ceps = _Ceps (Sta->pc, Sta->Sw);
    c   = Sta->Sw + Sta->n * Ceps;
    C   = -Sta->n * Cpc;
    chi = Sta->Sw;
    kwb = mkw(Sta->Sw) * kwb_sat;
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

inline double UnsatFlow::mkw (double Sw) const
{
    if (Sw<1.0)
    {
        if (Sw>0.0)
        {
            switch (wrc)
            {
                case BC_t: { return pow(Sw,akw); }
                case ZI_t:
                {
                    if (Sw>zi_del)
                    {
                        double hc = pow((1.0-Sw)/(Sw-zi_del), 1.0/zi_gam)/zi_bet; // m
                        return pow(1.0+pow(zi_a*hc, zi_b), -zi_alp);
                    }
                    else return 0.0;
                }
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
        switch (wrc)
        {
            case BC_t:
            {
                if (pc>bc_sb && Sw>bc_wr) return -bc_lam*pow(bc_sb/pc,bc_lam)*(1.0-bc_wr)/pc;
                else return 0.0;
            }
            case ZI_t:
            {
                double hc = pc/GamW; // m
                double K  = pow(zi_bet*hc, zi_gam);
                return -(1.0/GamW)*(zi_gam*K*(1.0-zi_del))/(hc*pow(K+1.0,2.0));
            }
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

inline void UnsatFlow::GenCurve (Array<double> pcs, char const * Filekey, size_t Np)
{
    // initial state
    UnsatFlowState sta(NDim);
    SDPair ini;
    ini.Set ("n",   Por);
    ini.Set ("pw", -pcs[0]);
    InitIvs (ini, &sta);

    // header and initial output
    std::ostringstream oss;
    oss << Util::_8s<<"pc"   << Util::_8s<<"Sw"   << Util::_8s<<"n"   << std::endl;
    oss << Util::_8s<<sta.pc << Util::_8s<<sta.Sw << Util::_8s<<sta.n << std::endl;

    // update
    for (size_t i=1; i<pcs.Size(); ++i)
    {
        double dpc = (pcs[i]-sta.pc)/Np;
        for (size_t j=0; j<Np; ++j)
        {
            Update (-dpc, 0.0, &sta);
            oss << Util::_8s<<sta.pc << Util::_8s<<sta.Sw << Util::_8s<<sta.n << std::endl;
        }
    }

    // save file
    String buf(Filekey);
    buf += ".res";
    std::ofstream of(buf.CStr(), std::ios::out);
    of << oss.str();
    of.close();
    printf("\nFile <%s%s%s> written\n",TERM_CLR_BLUE_H,buf.CStr(),TERM_RST);
}

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * UnsatFlowMaker(int NDim, SDPair const & Prms) { return new UnsatFlow(NDim,Prms); }

int UnsatFlowRegister()
{
    ModelFactory   ["UnsatFlow"] = UnsatFlowMaker;
    MODEL.Set      ("UnsatFlow", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["UnsatFlow"].Resize(13);
    MODEL_PRM_NAMES["UnsatFlow"] = "kwsat",  "akw",
                                   "bc_lam", "bc_sb",  "bc_wr",
                                   "zi_del", "zi_bet", "zi_gam", "zi_a", "zi_b", "zi_alp",
                                   "BC",     "ZI";
    MODEL_IVS_NAMES["UnsatFlow"].Resize(2);
    MODEL_IVS_NAMES["UnsatFlow"] = "pw", "Sw";
    return 0;
}

int __UnsatFlow_dummy_int = UnsatFlowRegister();


#endif // MECHSYS_UNSATFLOW_H
