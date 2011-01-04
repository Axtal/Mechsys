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
#include <mechsys/models/elastoplastic.h>
#include <mechsys/models/anisoinvs.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/sgrid.h>
#include <mechsys/vtk/isosurf.h>

using std::cout;
using std::endl;
using Util::PI;
double const TRUE=1.0;

void Func (Vec3_t const & X, double & F, Vec3_t & V, void * UserData)
{
    F = 1.0e+10;
    V = 0.0, 0.0, 0.0;
    if (X(0)<0.0 && X(1)<0.0 && X(2)<0.0)
    {
        ElastoPlastic const * mdl = static_cast<ElastoPlastic*>(UserData);
        EquilibState sta(/*NDim*/3), gra(/*NDim*/3);
        sta.Sig = X(0), X(1), X(2), 0.0, 0.0, 0.0;
        F = mdl->YieldFunc (&sta);
        mdl->Gradients (&gra);
    }
}

void AnisoFunc (Vec3_t const & X, double & F, Vec3_t & V, void * UserData)
{
    F = 1.0e+10;
    V = 0.0, 0.0, 0.0;
    if (X(0)<0.0 && X(1)<0.0 && X(2)<0.0)
    {
        Vec3_t a(0.2, 0.3, 1.0);
        AnisoInvs AI(0.5, 0.3, a, true); // b, alpha, a, obliq
        Vec_t sig(6);
        sig = X(0), X(1), X(2), 0.0, 0.0, 0.0;
        AI.Calc (sig, true); // true => derivs
        double R = 0.2;
        F = AI.sq - AI.sp * R;
        Vec_t dFdSig(AI.dsqdSig - R*AI.dspdSig);
        V = dFdSig(0), dFdSig(1), dFdSig(2);
    }
}

int main (int argc, char **argv) try
{
    // number:  nx ny nz
    Array<int> N(30, 30, 60);
    int tst = 1;
    if (argc>1) tst  = atoi(argv[1]);
    if (argc>2) N[0] = atoi(argv[2]);
    if (argc>3) N[1] = atoi(argv[3]);
    if (argc>4) N[2] = atoi(argv[4]);
    if (N[0]<2) throw new Fatal("nx must be greater than 1");
    if (N[1]<2) throw new Fatal("ny must be greater than 1");
    if (N[2]<2) throw new Fatal("nz must be greater than 1");

    // limits
    Array<double> L(6);
    //     0   1     2   3     4    5
    //   pmi pma   qmi qma  thmi thma
    L =   -5, 10,    0, 15,  -PI,  PI;

    // grid
    VTK::SGrid gri(N, L);

    // rotate
    Vec3_t x, l;
    for (int i=0; i<gri.Size(); ++i)
    {
        gri.GetPoint (i, x);
        pqth2L       (x(0), x(1), x(2), l, "cam");
        gri.SetPoint (i, l);
    }

    // function
    if (tst==1) { gri.SetFunc (&AnisoFunc); }
    else
    {
        SDPair prms;
        //prms.Set("E nu fc c phi", 1000.0, 0.3, FAILCRIT("MC"), 1.0, 30.0);
        prms.Set("E nu MN phi", 1000.0, 0.3, TRUE, 30.0);
        ElastoPlastic mdl(/*NDim*/3, prms);
        gri.SetFunc (&Func, &mdl);
    }

    // filter vectors
    gri.FilterV (0.0, /*Tol*/1.0e-1); 

    // isosurf
    VTK::IsoSurf iso(gri);
    iso.ShowVectors = true;

    // window and axes
    VTK::Win  win;
    VTK::Axes axe(/*scale*/20, /*hydroline*/true, /*reverse*/true);
    win.SetViewDefault (/*reverse*/true);
    axe.AddTo (win);
    //gri.AddTo (win);
    iso.AddTo (win);
    win.Show  ();

    // write file
    gri.WriteVTK ("ysurf");
    cout << "file [1;34m<ysurf.vtk>[0m written" << endl;

    // end
    return 0;
}
MECHSYS_CATCH
