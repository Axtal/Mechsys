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

    Vec3_t Le, Ls;
    Vec_t P0,P1,P2, Q0,Q1,Q2;
    EigenProj (sta->Eps, Le, Q0,Q1,Q2, /*sort*/true);
    EigenProj (sta->Sig, Ls, P0,P1,P2, /*sort*/true);

    Mat3_t X, Y, Yi; // devedgamdL, dpqthdL, dLdpqth
    double ev,ed,te;
    double p, q, ts;
    OctDerivs (Le, ev, ed, te, X);
    OctDerivs (Ls, p,  q,  ts, Y);
    InvOctDerivs (Ls, p,  q,  ts, Yi);
    double x = log(1.0+p);

    X(0,0) *= -1.;
    X(0,1) *= -1.;
    X(0,2) *= -1.;

    ev *= 100.;
    ed *= 100.;

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
    double A   = 1.0;
    double B   = 1.0;

    double a = -100.*(1.0+p)/(lam*A);
    double b = -100.*psi*B*(1.0+p)/(lam*A);
    double d =  100.*g;
    double m = 1.0;

    printf("a=%g, b=%g, d=%g\n",a,b,d);


    Mat3_t Z;
    Z =  a,   b,   0., 
         0.,  d,   0., 
         0.,  0.,  m;

    Mat3_t tm, E;
    //Inv (Y, Yi);
    tm = product (Yi, Z); // tm = Yi*Z
    E  = product (tm, X); // E  = tm*X = Yi*Z*X

    std::cout << "E=\n"   << PrintMatrix(E)        << std::endl;

    if (q<1.0e-7)
    {
        E = K+(4.0*G)/3.0,  K-(2.0*G)/3.0,  K-(2.0*G)/3.0,
            K-(2.0*G)/3.0,  K+(4.0*G)/3.0,  K-(2.0*G)/3.0,
            K-(2.0*G)/3.0,  K-(2.0*G)/3.0,  K+(4.0*G)/3.0;
    }

    std::cout << "sig="   << PrintVector(sta->Sig) << std::endl;
    std::cout << "X=\n"   << PrintMatrix(X)        << std::endl;
    std::cout << "Y=\n"   << PrintMatrix(Y)        << std::endl;
    std::cout << "Yi=\n"  << PrintMatrix(Yi)       << std::endl;
    std::cout << "Z=\n"   << PrintMatrix(Z)        << std::endl;
    std::cout << "E=\n"   << PrintMatrix(E)        << std::endl;

    D.change_dim (NCps,NCps);
    set_to_zero  (D);
    for (size_t i=0; i<NCps; ++i)
    for (size_t j=0; j<NCps; ++j)
    {
        //if (i<3 && j<3) D(i,j) = E(i,j);
        D(i,j) = P0(i)*E(0,0)*P0(j) + P0(i)*E(0,1)*P1(j) + P0(i)*E(0,2)*P2(j) + 
                 P1(i)*E(1,0)*P0(j) + P1(i)*E(1,1)*P1(j) + P1(i)*E(1,2)*P2(j) + 
                 P2(i)*E(2,0)*P0(j) + P2(i)*E(2,1)*P1(j) + P2(i)*E(2,2)*P2(j);
        if (i>2 && j>2 && i==j) D(i,j) += 1.0;
    }

    //Mat_t De;
    //De = (2.0*G)*Psd + K*IdyI;
    //std::cout << "De=\n" << PrintMatrix(De) << std::endl;
    std::cout << "D=\n"  << PrintMatrix(D)  << std::endl;
    //D = De;
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
