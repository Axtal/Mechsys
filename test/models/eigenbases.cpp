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

// Std Lib
#include <iostream>

// MechSys
#ifdef USE_VTK
  #include <mechsys/vtk/arrow.h>
  #include <mechsys/vtk/axes.h>
  #include <mechsys/vtk/cube.h>
  #include <mechsys/vtk/cylinder.h>
  #include <mechsys/vtk/plane.h>
  #include <mechsys/vtk/sphere.h>
  #include <mechsys/vtk/spheres.h>
  #include <mechsys/vtk/text.h>
  #include <mechsys/vtk/win.h>
#endif
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using Util::SQ2;
using OrthoSys::O;
using OrthoSys::e0;
using OrthoSys::e1;
using OrthoSys::e2;
using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    /* these give different directions for eigenvectors:
        ./eigenbases 1 1 1 0.00000001
        ./eigenbases 1 1 1 0.0000001
    */

    // stress
    Vec_t sig(4);
    sig = 1., 2., 1., 0.5*SQ2;
    //sig = 1., 2., 1., 1.; // <<<<<< one of the principal values is zero => problems with nsmp

    // aniso vector
    Vec3_t a;
    a = 1., 0., 0.;

    // input
    double alpha = 0.5;
    if (argc>1) sig(0) = atof(argv[1]);
    if (argc>2) sig(1) = atof(argv[2]);
    if (argc>3) sig(2) = atof(argv[3]);
    if (argc>4) sig(3) = atof(argv[4])*SQ2;
    if (argc>5) a(0)   = atof(argv[5]);
    if (argc>6) a(1)   = atof(argv[6]);
    if (argc>7) a(2)   = atof(argv[7]);
    if (argc>8) alpha  = atof(argv[8]);

    // eigenvs
    Mat3_t sigm;
    Vec3_t L, v0, v1, v2;
    Ten2Mat (sig, sigm);
    Eig (sigm, L, v0, v1, v2);
    Mat3_t Q, Qt;
    Q = v0(0), v1(0), v2(0),
        v0(1), v1(1), v2(1),
        v0(2), v1(2), v2(2);
    Trans (Q, Qt);

    // eigen projectors
    Vec3_t La;          // analytic
    Mat3_t P0,P1,P2;
    Vec_t  P0a,P1a,P2a; // analytic
    Dyad (v0,v0, P0);
    Dyad (v1,v1, P1);
    Dyad (v2,v2, P2);
    EigenProjAnalytic (sig, La, P0a, P1a, P2a);

    Mat3_t res, P0i;
    Inv (P0,P0i);
    res = product(P0i,P0);
    cout << "P0i*P0 =\n" << PrintMatrix(res);
    return 0;

    // check
    Mat3_t Lmat, sig_, tmp;
    Lmat = L(0), 0.0,  0.0,
           0.0,  L(1), 0.0,
           0.0,  0.0,  L(2);
    tmp  = product(Q,Lmat);
    sig_ = product(tmp,Qt);
    for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
    {
        if (fabs(sigm(i,j)-sig_(i,j))>1.0e-15) throw new Fatal("sig_ != sigm");
    }

    // normal to SMP
    double I1,I2,I3;
    CharInvs (sig, I1, I2, I3);
    Vec3_t nsmp_123, nsmp;
    nsmp_123 = 1./sqrt(L[0]), 
               1./sqrt(L[1]), 
               1./sqrt(L[2]);
    if (I3<0.0)
    {
        if (fabs(I3)<1.0e-12) I3 = 0.0;
        else throw new Fatal("I3=%g is invalid (negative)", I3);
    }
    nsmp_123 *= sqrt(I3/I2);
    nsmp      = product(Q,nsmp_123);

    // aniso vector w.r.t eigenbasis
    Vec3_t a_123, b;
    a_123 = product(Qt,a);
    b     = alpha*a + nsmp;

    // output
    cout << "sig      ="   << PrintVector(sig);
    cout << "sigm     =\n" << PrintMatrix(sigm);
    cout << "Q*L*Qt   =\n" << PrintMatrix(sig_);
    cout << "Lmat     =\n" << PrintMatrix(Lmat);
    cout << "L        ="   << PrintVector(L);
    cout << "La       ="   << PrintVector(La);
    cout << "Q        =\n" << PrintMatrix(Q);
    cout << "Qt       =\n" << PrintMatrix(Qt);
    cout << "v0       ="   << PrintVector(v0);
    cout << "v1       ="   << PrintVector(v1);
    cout << "v2       ="   << PrintVector(v2);
    cout << "I1,I2,I3 = "  << I1 << ",  " << I2 << ",  " << I3 << endl;
    cout << "nsmp_123 ="   << PrintVector(nsmp_123);
    cout << "nsmp     ="   << PrintVector(nsmp);
    cout << "a        ="   << PrintVector(a);
    cout << "b        ="   << PrintVector(b);
    cout << "a_123    ="   << PrintVector(a_123);
    cout << "P0       =\n" << PrintMatrix(P0);
    cout << "P1       =\n" << PrintMatrix(P1);
    cout << "P2       =\n" << PrintMatrix(P2);
    cout << "P0a      ="   << PrintVector(P0a);
    cout << "P1a      ="   << PrintVector(P1a);
    cout << "P2a      ="   << PrintVector(P2a);

    // visualisation
#ifdef USE_VTK
    VTK::Arrow ar0(O, v0);   ar0.SetColor("red");
    VTK::Arrow ar1(O, v1);   ar1.SetColor("green");
    VTK::Arrow ar2(O, v2);   ar2.SetColor("blue");
    VTK::Arrow ar3(O, nsmp); ar3.SetColor("black");
    VTK::Arrow ar4(O, a);    ar4.SetColor("magenta");
    VTK::Arrow ar5(O, b);    ar5.SetColor("violet");
    VTK::Axes  ax(/*scale*/1, /*hydroline*/true);
    VTK::Win   win;
    ax .AddTo (win);
    ar0.AddTo (win);
    ar1.AddTo (win);
    ar2.AddTo (win);
    ar3.AddTo (win);
    ar4.AddTo (win);
    ar5.AddTo (win);
    win.Show();
#endif

    return 0;
}
MECHSYS_CATCH
