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

// STL
#include <iostream>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/arrows.h>
#include <mechsys/vtk/cylinder.h>
#include <mechsys/vtk/sgrid.h>
#include <mechsys/vtk/isosurf.h>

using std::cout;
using std::endl;
using Util::PI;
using Util::SQ2;
using Util::SQ3;

const double PHI   = 30.0*PI/180.0;
const double SPHI  = sin(PHI);
const double ALPHA = (1.0+SPHI)/(1.0-SPHI);
const double kDP   = 2.0*SQ2*SPHI/(3.0-SPHI);
const double kMN   = (9.0-SPHI*SPHI)/(1.0-SPHI*SPHI);
const double kSMP  = 2.0*SQ2*tan(PHI)/3.0;
const double POCT0 = 0.0;

inline void NormalOCT (Vec3_t const & L, Vec3_t & n, Mat3_t & dndL)
{
    n = 1.0/SQ3, 1.0/SQ3, 1.0/SQ3;
    dndL = 0., 0., 0.,
           0., 0., 0.,
           0., 0., 0.;
}

inline void NormalSMP (Vec3_t const & L, Vec3_t & n, Mat3_t & dndL)
{
    // check
    if (!(L(0)<0.0 && L(1)<0.0 && L(2)<0.0)) throw new Fatal("SMP: This method only works for Li<0 (negative principal values of stress)");

    // characteristics invariants
    double I1 = L(0)+L(1)+L(2);
    double I2 = L(0)*L(1) + L(1)*L(2) + L(2)*L(0);
    double I3 = L(0)*L(1)*L(2);

    // normal to SMP
    n = sqrt(I3/(I2*L(0))),
        sqrt(I3/(I2*L(1))),
        sqrt(I3/(I2*L(2)));

    // derivative of n w.r.t L
    double c0 = I3/(2.0*n(0)*L(0)*I2);
    double c1 = I3/(2.0*n(1)*L(1)*I2);
    double c2 = I3/(2.0*n(2)*L(2)*I2);
    dndL = c0*(        -(I1-L(0))/I2), c0*(1.0/L(1)-(I1-L(1))/I2), c0*(1.0/L(2)-(I1-L(2))/I2),
           c1*(1.0/L(0)-(I1-L(0))/I2), c1*(        -(I1-L(1))/I2), c1*(1.0/L(2)-(I1-L(2))/I2),
           c2*(1.0/L(0)-(I1-L(0))/I2), c2*(1.0/L(1)-(I1-L(1))/I2), c2*(        -(I1-L(2))/I2);
}

inline void SigTauDerivs (Vec3_t const & L, Vec3_t const & n, Mat3_t const & dndL, double & sig, double & tau, Vec3_t & dsigdL, Vec3_t & dtaudL)
{
    // invariants
    double sb = L(0)*n(0)*n(0) + 
                L(1)*n(1)*n(1) + 
                L(2)*n(2)*n(2);
    sig = fabs(sb);
    tau = sqrt(pow(L(0)-sb,2.0)*n(0)*n(0) + 
               pow(L(1)-sb,2.0)*n(1)*n(1) + 
               pow(L(2)-sb,2.0)*n(2)*n(2));

    // derivative of sig_bar w.r.t L and of tau w.r.t L
    Vec3_t dsbdL(n(0)*n(0),
                 n(1)*n(1),
                 n(2)*n(2));
    dtaudL = n(0)*n(0)*(L(0)-sb)/tau,
             n(1)*n(1)*(L(1)-sb)/tau,
             n(2)*n(2)*(L(2)-sb)/tau; 
    for (int k=0; k<3; ++k)
    {
        for (int i=0; i<3; ++i) dsbdL (k) += 2.0*dndL(i,k)*n(i)*L(i);
        for (int i=0; i<3; ++i) dtaudL(k) += (dndL(i,k)*n(i)*pow(L(i)-sb,2.0) - dsbdL(k)*n(i)*n(i)*(L(i)-sb))/tau;
    }

    // derivative of sig w.r.t L
    if (sb>0.0) dsigdL =  dsbdL;
    else        dsigdL = -dsbdL;
}

struct Data
{
    Data () : Oct(false), Aniso(false), Obliq(true), Check(true), c(0.0), Tol(1.0e-6), a(0,0,1), k(-1) {}
    bool   Oct;
    bool   Aniso;
    bool   Obliq;
    bool   Check;
    double c;
    double Tol;
    Vec3_t a;
    double k;
    double CalcK ()
    {
        double M = 1.0/sqrt(2.0*ALPHA+1.0);
        double N = sqrt(ALPHA)*M;
        //double K = (2.0*N*N+(a(2)+a(1))*c*N+ALPHA*M*M+a(0)*ALPHA*c*M)/((a(2)+a(1))*c*N+a(0)*c*M+1.0);
        double K = (2.0*N*N+(a(1)+a(0))*c*N+ALPHA*M*M+a(2)*ALPHA*c*M)/((a(1)+a(0))*c*N+a(2)*c*M+1.0);

        //return sqrt(pow(ALPHA-K,2.0)*M*M + 2.0*pow(1.0-K,2.0)*N*N)/K;
        return sqrt(pow(ALPHA-K,2.0)*M*M + 2.0*pow(1.0-K,2.0)*N*N)/K;
    }
};

void FuncSMP (Vec3_t const & X, double & F, Vec3_t & V, void * UserData)
{
    V = 0,0,0;
    Vec3_t x(X-POCT0/SQ3);
    if (x(0)<0.0 && x(1)<0.0 && x(2)<0.0)
    {
        double I1   = x(0)+x(1)+x(2);
        double I2   = x(0)*x(1) + x(1)*x(2) + x(2)*x(0);
        double I3   = x(0)*x(1)*x(2);
        double ssmp = -3.0*I3/I2;
        double tsmp = sqrt(I1*I3/I2-pow(3.0*I3/I2,2.0));
        F = tsmp/ssmp - kSMP;
    }
    else F = 1.0e+10;
}

void Func (Vec3_t const & X, double & F, Vec3_t & V, void * UserData)
{
    // data
    Data * dat = static_cast<Data*>(UserData);

    // sigma modified
    Vec3_t x(X-POCT0/SQ3);

    // variables
    Vec3_t n;
    Mat3_t dndL;
    Vec3_t dsigdL, dtaudL;
    double sig, tau;
    double k = (dat->Oct ? sqrt(kMN/9.0-1.0) : sqrt(kMN/9.0-1.0));

    if (x(0)<0.0 && x(1)<0.0 && x(2)<0.0)
    {
        // normal and derivatives
        if (dat->Oct) NormalOCT (x, n, dndL);
        else          NormalSMP (x, n, dndL);

        if (dat->Aniso)
        {
            Vec3_t t(x(0)*n(0), x(1)*n(1), x(2)*n(2));
            Vec3_t b = dat->c * dat->a + n;
            Vec3_t p,q;
            if (dat->Obliq)
            {
                double A = dot(t,b);
                double B = dot(b,n);
                p = (A/B)*n;
            }
            else
            {
                n = b/norm(b);
                p = dot(t,b)*b;
            }
            q   = t - p;
            sig = norm(p);
            tau = norm(q);
            V   = 0,0,0;

            if (dat->k>0) k = dat->k;
            else          k = dat->CalcK();
        }
        else
        {
            // invariants and derivatives
            SigTauDerivs (x, n, dndL, sig, tau, dsigdL, dtaudL);
            V = (1.0/sig)*dtaudL - (tau/(sig*sig))*dsigdL;
            V /= norm(V);
            
            // check
            if (dat->Check)
            {
                // stresses on plane n
                Vec3_t t(x(0)*n(0), x(1)*n(1), x(2)*n(2));
                Vec3_t p(dot(t,n)*n);
                Vec3_t q(t-p);
                double sign = norm(p);
                double taun = norm(q);
                if ((fabs(sign-sig) > dat->Tol) || (fabs(taun-tau) > dat->Tol)) throw new Fatal("FuncSMP: Stress invariants are different\n\tsign != sig (%g != %g => Error = %g)\nor\ttaun != tau (%g != %g => Error = %g)",sign,sig,fabs(sign-sig),taun,tau,fabs(taun-tau));

                if (dat->Oct)
                {
                }
                else // SMP
                {
                    double I1   = x(0)+x(1)+x(2);
                    double I2   = x(0)*x(1) + x(1)*x(2) + x(2)*x(0);
                    double I3   = x(0)*x(1)*x(2);
                    double ssmp = -3.0*I3/I2;
                    double tsmp = sqrt(I1*I3/I2-pow(3.0*I3/I2,2.0));
                    if ((fabs(ssmp-sig) > dat->Tol) || (fabs(tsmp-tau) > dat->Tol)) throw new Fatal("FuncSMP: Stress invariants are different\n\tssmp != sig (%g != %g => Error = %g)\nor\ttsmp != tau (%g != %g => Error = %g)",ssmp,sig,fabs(ssmp-sig),tsmp,tau,fabs(tsmp-tau));
                }
            }
        }

        //Vec3_t n(1.0/x(0), 1.0/x(1), 1.0/x(2));
        //Vec3_t n(-1.0/pow(-x(0),1.0/4.0), -1.0/pow(-x(1),1.0/4.0), -1.0/pow(-x(2),1.0/4.0));
        //Vec3_t n(-1.0/sqrt(-x(0)), -1.0/sqrt(-x(1)), -1.0/sqrt(-x(2)));
        //n /= norm(n);
        //cout << n << endl;
        //n *= sqrt(-I3/I2);

        //double qoct = sqrt(pow(X(0)-X(1),2.0) + pow(X(1)-X(2),2.0) + pow(X(2)-X(0),2.0))/SQ3;
        //double poct = -(X(0)+X(1)+X(2))/SQ3;
        //cout << poct << " " << sig << "  " << qoct << " " << tau << endl;

        //F = I1*I2/I3 - kMN;
        //F = tsmp/ssmp - sqrt(kMN/9.0-1.0);
        F = tau/sig - k;
        //F = tau/sig - kMN;
        //F = qoct/poct - kDP;
    }
    else
    {
        F = 1.0e+15;
        V = 0,0,0;
    }
}

int main(int argc, char **argv) try
{
    // data
    Data dat;

    // number:  nx ny nz
    Array<int> N(10, 10, 30);
    //double scale = 6.0;
    if (argc> 1) N[0]      = atoi(argv[ 1]);
    if (argc> 2) N[1]      = atoi(argv[ 2]);
    if (argc> 3) N[2]      = atoi(argv[ 3]);
    if (argc> 4) dat.Oct   = atoi(argv[ 4]);
    if (argc> 5) dat.Aniso = atoi(argv[ 5]);
    if (argc> 6) dat.Obliq = atoi(argv[ 6]);
    if (argc> 7) dat.c     = atof(argv[ 7]);
    if (argc> 8) dat.k     = atof(argv[ 8]);
    if (argc> 9) dat.a(0)  = atof(argv[ 9]);
    if (argc>10) dat.a(1)  = atof(argv[10]);
    if (argc>11) dat.a(2)  = atof(argv[11]);
    //if (argc>12) scale     = atof(argv[12]);
    if (N[0]<2) throw new Fatal("nx must be greater than 1");
    if (N[1]<2) throw new Fatal("ny must be greater than 1");
    if (N[2]<2) throw new Fatal("nz must be greater than 1");

    cout << "k(mn) = " << sqrt(kMN/9.0-1.0) << endl;
    cout << "k     = " << dat.k << endl;
    cout << "CalcK = " << dat.CalcK() << endl;

    // limits
    double pcam_min = 0.1;
    double pcam_max = 10.0;
    Array<double> L(pcam_min,pcam_max,  0,15,  -PI,PI);

    // grid
    VTK::SGrid gr0(N, L);
    VTK::SGrid gr1(N, L);
    //gr0.ShowPoints ();

    // rotate
    Vec3_t x, l;
    for (int i=0; i<gr0.Size(); ++i)
    {
        gr0.GetPoint (i, x);
        pqth2L       (x(0), x(1), x(2), l, "cam");
        gr0.SetPoint (i, l);
        gr1.SetPoint (i, l);
    }

    // set function
    gr0.SetFunc (&Func,    &dat);
    gr1.SetFunc (&FuncSMP, &dat);
    gr0.FilterV (0.0, 1.0e-2);
    gr1.FilterV (0.0, 1.0e-2);

    // arrows
    //VTK::Arrows arr(gr0);
    //arr.SetScale (scale);
    //arr.SetColor ("red");

    // cylinder
    double poct_max = pcam_max*sqrt(3.0);
    double rad = kDP*poct_max;
    VTK::Cylinder cyl(Vec3_t(0,0,0), Vec3_t(-pcam_max,-pcam_max,-pcam_max), rad, false);
    cyl.SetColor ("black", 0.5);
    cyl.SetWire  ();

    // isosurface
    VTK::IsoSurf is0(gr0);
    VTK::IsoSurf is1(gr1);
    is0.SetColor ("red");
    is1.SetColor ("black",1.0);
    is1.SetWire  ();

    // window and axes
    VTK::Win  win;
    VTK::Axes axe(/*scale*/20, /*hydroline*/true, /*reverse*/true);
    win.SetViewDefault (/*reverse*/true);
    axe.AddTo (win);
    //gr0.AddTo (win);
    //arr.AddTo (win);
    is0.AddTo (win);
    is1.AddTo (win);
    cyl.AddTo (win);
    win.Show  ();
    win.WritePNG ("fcrit.png");

    // end
    return 0;
}
MECHSYS_CATCH
