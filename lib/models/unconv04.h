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

#ifndef MECHSYS_UNCONV04_H
#define MECHSYS_UNCONV04_H

// MechSys
#include <mechsys/models/model.h>

using std::cout;
using std::endl;

class Unconv04 : public Model
{
public:
    // Constructor
    Unconv04 (int NDim, SDPair const & Prms);

    // Derived methods
    void InitIvs    (SDPair const & Ini, State * Sta)                             const;
    void TgIncs     (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const;
    void Stiffness  (State const * Sta, Mat_t & D)                                const;
    bool LoadCond   (State const * Sta, Vec_t const & DEps, double & alpInt)      const;
    void UpdatePath (State const * Sta, Vec_t const & DEps, Vec_t const & DSig)   const;

    // Internal methods
    void Ref (double x, double a, double b, double c, double A, double B, double bet, double x0, double y0, double & D, double & lam, double & y) const;

    // Parameters
    double lam0, lam1, lam2, x1, x2, bet0, bet1;
    double psi0, psi1, ev1, ev2, bet2, bet3;
    double g0, g1, Mcs, Mso, bet4, bet5;
    double K, G;

    // Auxiliar data
    Vec_t I;
    Mat_t IdyI, Psd;

    // Modifiable data
    mutable double alpha;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Unconv04::Unconv04 (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"Unconv04"), alpha(0.0)
{
    lam0 = Prms("lam0");
    lam1 = Prms("lam1");
    lam2 = Prms("lam2");
    x1   = Prms("x1");
    x2   = Prms("x2");
    bet0 = Prms("bet0");
    bet1 = Prms("bet1");

    psi0 = Prms("psi0");
    psi1 = Prms("psi1");
    ev1  = Prms("ev1");
    ev2  = Prms("ev2");
    bet2 = Prms("bet2");
    bet3 = Prms("bet3");

    g0   = Prms("g0");
    g1   = Prms("g1");
    Mcs  = Prms("Mcs");
    Mso  = Prms("Mso");
    bet4 = Prms("bet4");
    bet5 = Prms("bet5");

    K    = Prms("K");
    G    = Prms("G");

    Calc_I    (NCps, I);
    Calc_IdyI (NCps, IdyI);
    Calc_Psd  (NCps, Psd);
}

inline void Unconv04::InitIvs (SDPair const & Ini, State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini);
}

inline void Unconv04::TgIncs (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const
{
    Mat_t D;
    Stiffness (Sta, D);
    DSig = D*DEps;
}

inline void Unconv04::Stiffness (State const * Sta, Mat_t & D) const
{
    EquilibState const * sta = static_cast<EquilibState const*>(Sta);
    double p,q,t,ev,ed;
    OctInvs (sta->Sig, p,q,t);
    ev = Calc_evoct(sta->Eps)*100.;
    ed = Calc_edoct(sta->Eps)*100.;
    double x = log(1.0+p);

    double D1,D3,D5,r0,r1,r2,lr0,lr1,lr2;
    //   x  , a    , b    , c        , A     , B    , bet  , x0  , y0    , D  , lam , y
    Ref (x  , lam1 , 1.0  , -lam1*x1 , lam2  , lam1 , bet1 , x2  , 0.0   , D1 , lr0 , r0);
    Ref (ed , 0.0  , -1.0 , ev2      , -psi1 , 0.0  , bet3 , 0.0 , ev1   , D3 , lr1 , r1);
    //Ref (ed , 0.0  , 1.0  , -Mcs   , g1    , 0.0  , bet5 , 0.0 , Mso , D5 , lr2 , r2);
    Ref (ed , 0.0  , 1.0  , -Mcs*p   , g1    , 0.0  , bet5 , 0.0 , Mso*p , D5 , lr2 , r2);

    double D0  = r0 - ev;
    double D2  = ev - r1;
    double D4  = r2 - q;
    double lam = lam0 + ( lr0 - lam0)*exp(-bet0*D0);
    double psi = psi0 + ( lr1 - psi0)*exp(-bet2*D2);
    double g   = g0   + (-lr2 - g0  )*exp(-bet4*D4);

    printf("psi=%g\n",psi);
    //printf("D1 = %g\n",D1);
    //printf("alpha = %g\n",alpha*180.0/Util::PI);

    double a  = -lam/(3.0*(1.0+p))/100.;
    double A  = -a/Util::SQ3;

    Mat_t C(A*IdyI + (0.5/G)*Psd);
    if (q>1.0e-10)
    {
        double b = -psi/(3.0*g);
        double c = 1.0/g;
        Vec_t B, S;
        Mat_t BdyS;
        Dev (sta->Sig, S);
        B = (b/q)*I + (c/q)*S;
        Dyad (B,S,BdyS);
        C += BdyS;
    }
    //std::cout << "C =\n" << PrintMatrix(C);
    Inv (C,D);
    //std::cout << "D =\n" << PrintMatrix(D);

    //double Kb = K/(1.0+9.*K*A);
    //D = (2.0*G)*Psd + Kb*IdyI;

    //printf("q=%g\n",q);

    //if (q>1.0e-10)
    //{
        //double b = -psi/(3.0*g);
        //double c = 1.0/g;
        //printf("g0=%g, lr2=%g, D4=%g, bet4=%g, exp(-bet4*D4)=%g, g=%g\n",g0,lr2,D4,bet4,exp(-bet4*D4),g);
        //printf("psi=%g, alpha=%g, g=%g, b=%g, c=%g\n",psi,alpha*180./Util::PI,g,b,c);
        //Vec_t B, S;
        //Dev (sta->Sig, S);
        //B = (b/q)*I + (c/q)*S;
        //Vec_t DeB(D*B);
        //Vec_t SDe;
        //Mult (S,D, SDe);
        //double phi = 1.0 + dot(S,DeB);
        //for (size_t i=0; i<NCps; ++i)
        //for (size_t j=0; j<NCps; ++j)
            //D(i,j) -= DeB(i)*SDe(j)/phi;
        //std::cout << PrintMatrix(D);
    //}
}

inline bool Unconv04::LoadCond (State const * Sta, Vec_t const & DEps, double & alpInt) const
{
    return true;
}

inline void Unconv04::UpdatePath (State const * Sta, Vec_t const & DEps, Vec_t const & DSig) const
{
    double dq = Calc_qoct (DSig);
    double dp = Calc_poct (DSig);
    alpha = atan2(dq,dp);
    //printf("alpha = %g\n",alpha*180./Util::PI);
}

inline void Unconv04::Ref (double x, double a, double b, double c, double A, double B, double bet, double x0, double y0, double & D, double & lam, double & y) const
{
    double c1 = bet*(b*A-a);
    double c2 = (A-B)*exp(-c*bet)/(A-a/b);
    double c3 = exp(b*bet*(y0+A*x0)) - c2*exp(c1*x0);
    //printf("c1=%g, c2=%g, c3=%g  ",c1,c2,c3);
    y   = -A*x + log(c3+c2*exp(c1*x))/(b*bet);
    D   = a*x + b*y + c;
    lam = A + (B - A)*exp(-bet*D);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * Unconv04Maker(int NDim, SDPair const & Prms) { return new Unconv04(NDim,Prms); }

int Unconv04Register()
{
    ModelFactory["Unconv04"] = Unconv04Maker;
    MODEL.Set ("Unconv04", (double)MODEL.Keys.Size());
    return 0;
}

int __Unconv04_dummy_int = Unconv04Register();

#endif // MECHSYS_UNCONV04_H
