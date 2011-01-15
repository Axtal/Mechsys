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
#include <cmath>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/models/anisoinvs.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>

// MechSys -- VTK
#include <mechsys/vtk/arrow.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/plane.h>
#include <mechsys/vtk/sphere.h>
#include <mechsys/vtk/win.h>

using std::cout;
using std::endl;
using Util::SQ2;
using Util::SQ6;
using Util::PI;

int main(int argc, char **argv) try
{
#ifdef HAS_TENSORS
    VTK::Win win;

    Vec_t sig0(6), sig1(6);
    sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2,  0.0;
    sig1 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.00000001*SQ2;
    //sig0 = -1.5, -2.0,  0.0,  -0.1*SQ2, 0.0, 0.0;
    //sig1 = -1.5, -2.0,  0.0,  -0.1*SQ2, 0.0, -0.00000001*SQ2;

    Vec3_t a(0.2, 0.3, 1.0);
    AnisoInvs AI(0.0, 0.0, a, true); // b, alpha, a, obliq
    cout << "a = " << PrintVector(a);
    cout << endl;

    Vec3_t zero(0., 0., 0.);
    VTK::Arrow vtk_a(zero, a);   vtk_a.SetColor("red");   vtk_a.AddTo(win);

    cout << "sig0 = " << PrintVector(sig0);
    cout << "sig1 = " << PrintVector(sig1);
    cout << endl;

    Vec3_t v0,v1,v2, L, La;
    Vec_t  E0,E1,E2, E0a,E1a,E2a;
    EigenProj         (sig0, L, v0,v1,v2, E0,E1,E2);
    EigenProjAnalytic (sig0, La, E0a,E1a,E2a);
    
    AI.Calc (sig0, /*derivs*/false);
    cout << "sig0: L0    = "  << AI.L(0) << endl;
    cout << "sig0: L1    = "  << AI.L(1) << endl;
    cout << "sig0: L2    = "  << AI.L(2) << endl;
    cout << "sig0: v0    = "  << PrintVector(AI.v0);
    cout << "sig0: v1    = "  << PrintVector(AI.v1);
    cout << "sig0: v2    = "  << PrintVector(AI.v2);
    cout << "sig0: t     = "  << PrintVector(AI.t);
    cout << "sig0: sp    = "  << AI.sp << endl;
    cout << "sig0: sq    = "  << AI.sq << endl;
    cout << "sig0: sq/sp = "  << AI.sq/AI.sp << endl;
    cout << endl;
    VTK::Arrow vtk_v0_0(zero, AI.v0);   vtk_v0_0.SetColor("yellow",0.5);   vtk_v0_0.AddTo(win);
    VTK::Arrow vtk_v1_0(zero, AI.v1);   vtk_v1_0.SetColor("yellow",0.5);   vtk_v1_0.AddTo(win);
    VTK::Arrow vtk_v2_0(zero, AI.v2);   vtk_v2_0.SetColor("yellow",0.5);   vtk_v2_0.AddTo(win);
    //VTK::Arrow vtk_Nu_0(zero, AI.Nu);   vtk_Nu_0.SetColor("yellow",0.5);   vtk_Nu_0.AddTo(win);
    //VTK::Arrow vtk_t_0 (zero, AI.t);    vtk_t_0 .SetColor("yellow",1.0);   vtk_t_0 .AddTo(win);
    //VTK::Arrow vtk_p_0 (zero, AI.p);    vtk_p_0 .SetColor("gold",0.5);     vtk_p_0 .AddTo(win);
    //VTK::Arrow vtk_q_0 (zero, AI.q);    vtk_q_0 .SetColor("gold",0.5);     vtk_q_0 .AddTo(win);
    VTK::Plane vtk_pl0 (zero, AI.Nu);   vtk_pl0 .SetColor("yellow",0.5);   vtk_pl0 .AddTo(win);

    Vec3_t v0_,v1_,v2_, L_, La_;
    Vec_t  E0_,E1_,E2_, E0a_,E1a_,E2a_;
    EigenProj         (sig1, L_, v0_,v1_,v2_, E0_,E1_,E2_);
    EigenProjAnalytic (sig1, La_, E0a_,E1a_,E2a_);

    AI.Calc (sig1, /*derivs*/false);
    cout << "sig1: L0    = " << AI.L(0) << endl;
    cout << "sig1: L1    = " << AI.L(1) << endl;
    cout << "sig1: L2    = " << AI.L(2) << endl;
    cout << "sig1: v0    = " << PrintVector(AI.v0);
    cout << "sig1: v1    = " << PrintVector(AI.v1);
    cout << "sig1: v2    = " << PrintVector(AI.v2);
    cout << "sig1: t     = " << PrintVector(AI.t);
    cout << "sig1: sp    = " << AI.sp << endl;
    cout << "sig1: sq    = " << AI.sq << endl;
    cout << "sig1: sq/sp = " << AI.sq/AI.sp << endl;
    VTK::Arrow vtk_v0_1(zero, AI.v0, 0.1,0.015,0.007);   vtk_v0_1.SetColor("blue");          vtk_v0_1.AddTo(win);
    VTK::Arrow vtk_v1_1(zero, AI.v1, 0.1,0.015,0.007);   vtk_v1_1.SetColor("blue");          vtk_v1_1.AddTo(win);
    VTK::Arrow vtk_v2_1(zero, AI.v2, 0.1,0.015,0.007);   vtk_v2_1.SetColor("blue");          vtk_v2_1.AddTo(win);
    //VTK::Arrow vtk_Nu_1(zero, AI.Nu                 );   vtk_Nu_1.SetColor("blue",  0.5);    vtk_Nu_1.AddTo(win);
    //VTK::Arrow vtk_t_1 (zero, AI.t);                     vtk_t_1 .SetColor("blue",  1.0);    vtk_t_1 .AddTo(win);
    //VTK::Arrow vtk_p_1 (zero, AI.p);                     vtk_p_1 .SetColor("cyan",  0.5);    vtk_p_1 .AddTo(win);
    //VTK::Arrow vtk_q_1 (zero, AI.q);                     vtk_q_1 .SetColor("cyan",  0.5);    vtk_q_1 .AddTo(win);
    VTK::Plane vtk_pl1 (zero, AI.Nu);                    vtk_pl1 .SetColor("blue",  0.5);    vtk_pl1 .AddTo(win);

    cout << endl;
    cout << "sig0: L     = " << PrintVector(L);
    cout << "sig1: L_    = " << PrintVector(L_);
    cout << "sig0: La    = " << PrintVector(La);
    cout << "sig1: La_   = " << PrintVector(La_);
    cout << endl;
    cout << "sig0: E0    = " << PrintVector(E0);
    cout << "sig1: E0    = " << PrintVector(E0_);
    cout << "sig0: E0a   = " << PrintVector(E0a);
    cout << "sig1: E0a   = " << PrintVector(E0a_);
    cout << endl;
    cout << "sig0: E1    = " << PrintVector(E1);
    cout << "sig1: E1    = " << PrintVector(E1_);
    cout << "sig0: E1a   = " << PrintVector(E1a);
    cout << "sig1: E1a   = " << PrintVector(E1a_);
    cout << endl;
    cout << "sig0: E2    = " << PrintVector(E2);
    cout << "sig1: E2    = " << PrintVector(E2_);
    cout << "sig0: E2a   = " << PrintVector(E2a);
    cout << "sig1: E2a   = " << PrintVector(E2a_);
    cout << endl;



    return 0;

    VTK::Axes ax(/*scale*/1, /*hydroline*/true, /*reverse*/false, /*full*/true);
    ax.AddTo (win);
    win.Show();

#else
    cout << "This test needs Tensors" << endl;
#endif
    return 0;
}
MECHSYS_CATCH
