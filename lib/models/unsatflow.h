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
#include<mechsys/models/unsatflowstate.h>
#include<mechsys/linalg/matvec.h>
#include<mechsys/numerical/odesolver.h>

class UnsatFlow : public Model
{
public:
    // static
    static Vec_t Iv; ///< Identity vector (NCo)

    // enums
    enum WRC_t { BC_t, ZI_t, HZ_t, PW_t }; ///< WRC type

    // Constructor
    UnsatFlow (int NDim, SDPair const & Prms, Model const * EquilibMdl=NULL);

    // Methods
    void   InitIvs   (SDPair const & Ini, State * Sta) const;                  ///< Initialize internal values
    void   Update    (double Dpw, double DEv, UnsatFlowState * Sta);           ///< Update state
    void   TgVars    (UnsatFlowState const * Sta) const;                       ///< Calculate c, C, chi, and kwb
    void   TgVars    (UnsatFlowState const * Sta, Vec_t const & Ww, double & Cpw, double & Cvs, Vec_t & fwd) const;
    void   TgIncs    (UnsatFlowState const * Sta, double Dpw, double DEv, double & DSw, double & Dchi) const;
    double FindSw    (double pc);                                              ///< Find Sw corresponding to pc by integrating from (pc,Sw)=(0,1) to pc (disregarding Dev)
    double Findpc    (double Sw);                                              ///< Find pc corresponding to Sw by integrating from (pc,Sw)=(0,1) to Sw (disregarding Dev)
    double rw        (double Sw) const;                                        ///< kwsat multiplier
    void   GenCurve  (Array<double> pcs, char const * Filekey, size_t Np=100); ///< Generate WRC
    void   Stiffness (State const * Sta, UnsatFlowState const * FSta, Mat_t & D) const; ///< Coupled "effective" stiffness
    void   RateAndUpdate (size_t Idx, Vec_t const & d, double dpwdt, double dt, UnsatFlowState * FSta, EquilibState * Sta) const;

    // Pointer to equilib model
    Model const * EMdl;

    // Data
    Mat_t  kwsatI; ///< Inverse of saturated kw
    Mat_t  kwsatb; ///< Saturated kw bar
    double akw;    ///< Exponent of kw model
    WRC_t  WRC;    ///< Water retention curve model
    String Name;   ///< Model name
    double Cw;     ///< Compressibility of water

    // SWRC parameters
    double bc_lam, bc_sb,  bc_wr;                      ///< Brooks & Corey model parameters
    double zi_del, zi_bet, zi_gam, zi_a, zi_b, zi_alp; ///< Zienkiewicz et al. 1990 (zi_bet and zi_a => [m])
    double hz_a,   hz_b,   hz_A,   hz_B;               ///< Huang and Zienkiewicz (Gawin, Schrefler ...)
    double ld,  xRd, yR, xRw, bd, bw, b1;              ///< Pedroso and Williams (PW)

    // Derived constants
    double y0, c1d, c2d, c3d, c1w, c2w, c3w; // Pedroso and Williams (PW)

    // Read-write variables (scratchpad)
    mutable double c, C, chi, Cpc, Ceps;

#define HMSTRESSUPDATE_DECLARE
    #include <mechsys/models/hmstressupdate.h>
    mutable HMStressUpdate HMSUp;
#undef HMSTRESSUPDATE_DECLARE

private:
    double _Cpc  (bool Drying, double pc, double Sw) const;
    double _Ceps (bool Drying, double pc, double Sw) const;

    double _RK_pc, _RK_Dev, _RK_Dpc, _RK_Sw, _RK_DSw;
    int _RK_func_Sw (double t, double const Sw[], double dSwdt[]);
    int _RK_func_pc (double t, double const pc[], double dpcdt[]);
};

Vec_t UnsatFlow::Iv;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


#define HMSTRESSUPDATE_IMPLEMENT
  #include <mechsys/models/hmstressupdate.h>
#undef HMSTRESSUPDATE_IMPLEMENT


inline UnsatFlow::UnsatFlow (int NDim, SDPair const & Prms, Model const * EquilibMdl)
    : Model (NDim,Prms,/*niv*/0,"UnsatFlow"), EMdl(EquilibMdl), WRC(PW_t), Name("Pedroso-Williams")
{
    // parameters
    if (Prms.HasKey("BC")) WRC = BC_t;
    if (Prms.HasKey("ZI")) WRC = ZI_t;
    if (Prms.HasKey("HZ")) WRC = HZ_t;
    if (Prms.HasKey("PW")) WRC = PW_t;
    switch (WRC)
    {
        case BC_t:
        {
            Name   = "Brooks-Corey";
            akw    = Prms("akw");
            bc_lam = (Prms.HasKey("bc_lam") ? Prms("bc_lam") : 0.8 );
            bc_sb  = (Prms.HasKey("bc_sb" ) ? Prms("bc_sb" ) : 1.8 );
            bc_wr  = (Prms.HasKey("bc_wr" ) ? Prms("bc_wr" ) : 0.01);
            if (Prms.HasKey("bc_sbb")) bc_sb = exp(Prms("bc_sbb"))-1.0;
            break;
        }
        case ZI_t:
        {
            Name   = "Zienkiewicz";
            zi_del = (Prms.HasKey("zi_del") ? Prms("zi_del") : 0.0842);
            zi_bet = (Prms.HasKey("zi_bet") ? Prms("zi_bet") : 0.7   );
            zi_gam = (Prms.HasKey("zi_gam") ? Prms("zi_gam") : 2.0   );
            zi_a   = (Prms.HasKey("zi_a"  ) ? Prms("zi_a"  ) : 5.0   );
            zi_b   = (Prms.HasKey("zi_b"  ) ? Prms("zi_b"  ) : 4.0   );
            zi_alp = (Prms.HasKey("zi_alp") ? Prms("zi_alp") : 0.9   );
            break;
        }
        case HZ_t:
        {
            Name = "Huang-Zienkiewicz";
            hz_a = (Prms.HasKey("hz_a"  ) ? Prms("hz_a"  ) : 0.10152);
            hz_b = (Prms.HasKey("hz_b"  ) ? Prms("hz_b"  ) : 2.4279 );
            hz_A = (Prms.HasKey("hz_A"  ) ? Prms("hz_A"  ) : 2.207  );
            hz_B = (Prms.HasKey("hz_B"  ) ? Prms("hz_B"  ) : 1.0121 );
            break;
        }
        case PW_t:
        {
            Name = "Pedroso-Williams";
            akw  = Prms("akw");
            ld   = (Prms.HasKey("ld" ) ? Prms("ld" ) : 2.8);
            xRd  = (Prms.HasKey("xRd") ? Prms("xRd") : 0.88);
            yR   = (Prms.HasKey("yR" ) ? Prms("yR" ) : 0.04);
            xRw  = (Prms.HasKey("xRw") ? Prms("xRw") : 0.7);
            bd   = (Prms.HasKey("bd" ) ? Prms("bd" ) : 3.0);
            bw   = (Prms.HasKey("bw" ) ? Prms("bw" ) : 5.0);
            b1   = (Prms.HasKey("b1" ) ? Prms("b1" ) : 1.8);

            // derived constants
            y0  = 1.0; // saturation for pc=0
            c1d = bd * ld;
            c2d = exp(bd * yR);
            c3d = exp(bd * (y0 + ld * xRd)) - c2d * exp(c1d * xRd);
            c1w = -bw * ld;
            c2w = exp(-bw * y0);
            c3w = exp(-bw * ld * xRw) - c2w * exp(c1w * xRw);
        }
    }

    // saturated conductivity matrix
    Mat_t kwsat(NDim, NDim);
    kwsatb.change_dim (NDim, NDim);
    GamW = Prms("gamW");
    double m = Prms("kwsat") / GamW;
    if (NDim==3) kwsat  =  Prms("kwsat"), 0., 0.,   0., Prms("kwsat"), 0.,   0., 0., Prms("kwsat");
    else         kwsat  =  Prms("kwsat"), 0.,   0., Prms("kwsat");
    if (NDim==3) kwsatb =  m, 0., 0.,   0., m, 0.,   0., 0., m;
    else         kwsatb =  m, 0.,   0., m;
    Inv (kwsat, kwsatI);

    // compressibility of water
    Cw = (Prms.HasKey("Cw") ? Prms("Cw") : 1.0e-7);

    // stress update
    HMSUp.SetModel (EquilibMdl, this);

    // Iv
    if (size(Iv)==0)
    {
        Iv.change_dim(NDim==2 ? 4 : 6);
        Iv(0) = 1.;
        Iv(1) = 1.;
        Iv(2) = 1.;
    }
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
    ini.Set   ("RhoW", Prms("gamW")/Prms("grav"));
    ini.Set   ("RhoS", Prms("rhoS"));
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
    Sta->Sw     = ode.Y[0];
    Sta->pc    += (-Dpw);
    Sta->n     += DEv;
    Sta->Drying = (_RK_Dpc>0.0); // drying/wetting ?
}

inline void UnsatFlow::TgVars (UnsatFlowState const * Sta) const
{
    Cpc  = _Cpc  (Sta->Drying, Sta->pc, Sta->Sw);
    Ceps = _Ceps (Sta->Drying, Sta->pc, Sta->Sw);
    c   = Sta->Sw + Sta->n * Ceps;
    C   = -Sta->n * Cpc;
    chi = Sta->Sw;
}

inline void UnsatFlow::TgVars (UnsatFlowState const * Sta, Vec_t const & Ww, double & Cpw, double & Cvs, Vec_t & fwd) const
{
    double gamw = Grav * Sta->RhoW;
    double nw   = Sta->n * Sta->Sw;
    double Cn   = 0.0; // == dSwdn

    Cpw = nw * Cw  -  Sta->n * Sta->RhoW * _Cpc(Sta->Drying, Sta->pc, Sta->Sw);
    Cvs = Sta->n * Sta->RhoW * Cn * (1.0-Sta->n)  +  Sta->Sw * Sta->RhoW;
    fwd = (-nw*nw*gamw/rw(Sta->Sw)) * kwsatI * Ww;
}

inline void UnsatFlow::RateAndUpdate (size_t Idx, Vec_t const & d, double dpwdt, double dt, UnsatFlowState * FSta, EquilibState * Sta) const
{
    // alloc vectors
    if (FSta->dndt.Size()<1)
    {
        FSta->dndt   .Resize (2);
        FSta->dSwdt  .Resize (2);
        FSta->dRhoWdt.Resize (2);
    }
    if (Sta->dSigdt.Size()<1)
    {
        Sta->dSigdt.Resize (2);
        for (size_t i=0; i<2; ++i) Sta->dSigdt[i].change_dim (NCps);
    }

    // rate
    double Cn = 0.0; // == dSwdn
    double pw = -FSta->pc;
    Mat_t  D;
    EMdl->Stiffness (Sta, D);
    FSta->dndt   [Idx] = (1.0-FSta->n)*Tra(d);
    FSta->dSwdt  [Idx] = -Cpc*dpwdt + Cn*FSta->dndt[Idx];
    FSta->dRhoWdt[Idx] = Cw*dpwdt;
    Sta ->dSigdt [Idx] = D*d - (dpwdt*FSta->Sw + pw*FSta->dSwdt[Idx])*I;

    // update
    if (Idx==0) // Forward-Euler
    {
        FSta->n    += FSta->dndt   [0]*dt;
        FSta->Sw   += FSta->dSwdt  [0]*dt;
        FSta->RhoW += FSta->dRhoWdt[0]*dt;
        Sta ->Sig  += Sta ->dSigdt [0]*dt;
    }
    else if (Idx==1) // Modified-Euler
    {
        FSta->n    += (0.5*dt)*(FSta->dndt   [0] + FSta->dndt   [1]);
        FSta->Sw   += (0.5*dt)*(FSta->dSwdt  [0] + FSta->dSwdt  [1]);
        FSta->RhoW += (0.5*dt)*(FSta->dRhoWdt[0] + FSta->dRhoWdt[1]);
        Sta ->Sig  += (0.5*dt)*(Sta ->dSigdt [0] + Sta ->dSigdt [1]);
    }
    else throw new Fatal("UnsatFlow::RateAndUpdate: Idx==%zd is invalid (0=FE, 1=ME)",Idx);
    FSta->pc     += (-dpwdt)*dt;
    FSta->Drying  = (FSta->dSwdt[0]<0 ? true : false);
    Sta ->Eps    += d*dt;
}

inline void UnsatFlow::TgIncs (UnsatFlowState const * Sta, double Dpw, double DEv, double & DSw, double & Dchi) const
{
    TgVars (Sta);
    DSw  = Ceps*DEv + Cpc*(-Dpw);
    Dchi = DSw;
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

inline double UnsatFlow::rw (double Sw) const
{
    if (Sw<1.0)
    {
        if (Sw>0.0)
        {
            switch (WRC)
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
                case HZ_t:
                {
                    double res = 1.0 - hz_A*pow(1.0-Sw, hz_B);
                    return (res<0.0 ? 0.0 : res);
                }
                case PW_t: { return pow(Sw,akw); }
                default: throw new Fatal("UnsatFlow::rw: WRC is invalid");
            }
        }
        else return 0.0;
    }
    else return 1.0;
}

inline double UnsatFlow::_Cpc (bool Drying, double pc, double Sw) const
{
    if (pc>0.0 && Sw>0.0)
    {
        switch (WRC)
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
            case HZ_t:
            {
                return -(hz_a*hz_b*pow(pc/GamW,hz_b))/pc;
            }
            case PW_t:
            {
                double x   = log(1.0+pc);
                double y   = Sw;
                double lb  = 0.0;
                if (Drying)
                {
                    double Dd  = (y-yR>0.0 ? y-yR : 0.0);
                    double ldb = ld * (1.0 - exp(-bd * Dd));
                    double yd  = -ld * x + log(c3d + c2d * exp(c1d * x)) / bd;
                    double D   = (yd-y>0.0 ? yd-y : 0.0);
                    double yb  = y/y0;
                    double b2  = bw;
                    double b2b = b2 * sqrt((yb>0.0 ? yb : 0.0));  // b2b = b2
                           lb  = ldb * exp(-b2b * D);
                }
                else
                {
                    double Dw  = (y0-y>0.0 ? y0-y : 0.0);
                    double lwb = ld * (1.0 - exp(-bw * Dw));
                    double yw  = -ld * x - log(c3w + c2w * exp(c1w * x)) / bw;
                    double D   = (y-yw>0.0 ? y-yw : 0.0);
                           lb  = lwb * exp(-b1 * D);
                }
                return -lb/(1.0+pc);
            }
            default: throw new Fatal("UnsatFlow::_Cpc: WRC is invalid");
        }
    }
    else return 0.0;
}

inline double UnsatFlow::_Ceps (bool Drying, double pc, double Sw) const
{
    return 0.0;
}

inline int UnsatFlow::_RK_func_Sw (double t, double const Y[], double dYdt[])
{
    double pc  = _RK_pc + t*_RK_Dpc;
    double Sw  = Y[0];
    bool   dry = (_RK_Dpc>0.0);
    Cpc     = _Cpc  (dry, pc, Sw);
    Ceps    = _Ceps (dry, pc, Sw);
    dYdt[0] = Ceps*_RK_Dev + Cpc*_RK_Dpc; // dSwdt
    //std::cout << Y[0] << "   " << dYdt[0] << std::endl;
    return GSL_SUCCESS;
}

inline int UnsatFlow::_RK_func_pc (double t, double const Y[], double dYdt[])
{
    double Sw  = _RK_Sw + t*_RK_DSw;
    double pc  = Y[0];
    bool   dry = (_RK_DSw<0.0);
    Cpc  = _Cpc  (dry, pc, Sw);
    Ceps = _Ceps (dry, pc, Sw);
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
    oss << Util::_8s<<"pc"   << Util::_8s<<"Sw"   << Util::_6<<"Drying"   << Util::_8s<<"n"   << Util::_8s<<"rw"       << Util::_8s<<"kw"                        << std::endl;
    oss << Util::_8s<<sta.pc << Util::_8s<<sta.Sw << Util::_6<<sta.Drying << Util::_8s<<sta.n << Util::_8s<<rw(sta.Sw) << Util::_8s<<rw(sta.Sw)*kwsatb(0,0)*GamW << std::endl;

    // update
    for (size_t i=1; i<pcs.Size(); ++i)
    {
        double dpc = (pcs[i]-sta.pc)/Np;
        for (size_t j=0; j<Np; ++j)
        {
            Update (-dpc, 0.0, &sta);
            oss << Util::_8s<<sta.pc << Util::_8s<<sta.Sw << Util::_6<<sta.Drying << Util::_8s<<sta.n << Util::_8s<<rw(sta.Sw) << Util::_8s<<rw(sta.Sw)*kwsatb(0,0)*GamW << std::endl;
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

inline void UnsatFlow::Stiffness (State const * Sta, UnsatFlowState const * FSta, Mat_t & D) const
{
    // eq state
    EquilibState sta(NDim);
    sta = (*static_cast<EquilibState const *>(Sta));

    // constitutive stresses
    TgVars (FSta);
    sta.Sig += (chi*(-FSta->pc))*Iv;

    // stiffness
    EMdl->Stiffness (&sta, D);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * UnsatFlowMaker(int NDim, SDPair const & Prms, Model const * EquilibMdl) { return new UnsatFlow(NDim,Prms,EquilibMdl); }

int UnsatFlowRegister()
{
    ModelFactory   ["UnsatFlow"] = UnsatFlowMaker;
    MODEL.Set      ("UnsatFlow", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["UnsatFlow"].Resize(26);
    MODEL_PRM_NAMES["UnsatFlow"] = "kwsat",  "akw",
                                   "bc_lam", "bc_sb",  "bc_wr",
                                   "zi_del", "zi_bet", "zi_gam", "zi_a", "zi_b", "zi_alp",
                                   "hz_a",   "hz_b",   "hz_A",   "hz_B",
                                   "ld",     "xRd",    "yR",     "xRw",  "bd",   "bw",  "b1",
                                   "BC",     "ZI",     "HZ",     "PW";
    MODEL_IVS_NAMES["UnsatFlow"].Resize(2);
    MODEL_IVS_NAMES["UnsatFlow"] = "pw", "Sw";
    return 0;
}

int __UnsatFlow_dummy_int = UnsatFlowRegister();


#endif // MECHSYS_UNSATFLOW_H
