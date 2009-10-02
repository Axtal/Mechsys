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

/*  Bhatti (2005): Example 1.4, p25  *
 *  ===============================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "fem/rod.h"
#include "fem/nlrod.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using Util::_6_3;
using Util::_8s;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    double l = 10.;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (3/*nodes*/, 2/*cells*/);
    mesh.SetVert  (0, -100,  0.0,  2.*l);
    mesh.SetVert  (1, -200,  0.0,     l);
    mesh.SetVert  (2, -300,  0.0,   0.0);
    mesh.SetCell  (0,   -1, /*NVerts*/2, 0, 1);
    mesh.SetCell  (1,   -1, /*NVerts*/2, 1, 2);
    mesh.WriteMPY ("owen_hinton_01",/*OnlyMesh*/false);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob E0 alp A fra", PROB("NLRod"), 200.0, 5.0, 1.0, TRUE);
    //prps.Set(-1, "prob E A fra", PROB("Rod"), 200.0, 1.0, TRUE);

    // domain
    FEM::Domain dom(/*NDim*/2, prps, /*mdls*/Dict(), /*inis*/Dict());
    dom.SetMesh (mesh);
    dom.SetOutNods ("owen_hinton_01",/*NNod*/1,/*IDs*/2);

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy", 0.0,0.0)
       .Set(-200, "ux",    0.0)
       .Set(-300, "ux fy", 0.0, -10.0);
    dom.SetBCs (bcs);

    FEM::Solver sol(dom);
    sol.Initialize ();
    sol.AssembleKA ();

    long neq = sol.NEQ;

    Vec_t F(neq), U(neq), Fext(neq), Fint(neq);
    set_to_zero (F);
    set_to_zero (U);
    set_to_zero (Fext);
    set_to_zero (Fint);

    Vec_t dU(neq);
    set_to_zero(dU);

    Vec_t R(neq);
    Vec_t dUc(neq);

    Array<double> Fincs(2);
    Fincs = -8.0, -2.0;

    // equation for output
    long eqx = dom.Nods[2]->EQ[1];

    // initial output
    std::ofstream of("owen_hinton_01.res", std::ios::out);
    of << _6_3 << "Time" << _8s << "u" << _8s << "fint" << _8s << "fext" << endl;
    of << _6_3 << sol.Time << _8s << U(eqx) << _8s << Fint(eqx) << _8s << Fext(eqx) << endl;

    Vec_t dF(neq);
    size_t NIncs = Fincs.Size();
    for (size_t inc=0; inc<NIncs; ++inc)
    {
        set_to_zero (dF);
        dF(eqx) = Fincs[inc];

        sol.AssembleKA  ();
        UMFPACK::Solve  (sol.A11, dF, dU);
        Sparse::AddMult (sol.K21, dU, dF); // dF2 += K21*dU1
        Sparse::AddMult (sol.K22, dU, dF); // dF2 += K22*dU2

        U += dU;
        F += dF;
        cout << "dU   = " << PrintVector(dU);

        for (size_t i=0; i<dom.Eles.Size(); ++i) dom.Eles[i]->UpdateState (dU, &Fint);
        Fext += dF;
        R     = Fext-Fint;
        double norm_R = Norm(R);

        // output
        sol.Time += 1.0;
        of << _6_3 << sol.Time << _8s << U(eqx) << _8s << Fint(eqx) << _8s << Fext(eqx) << endl;

        size_t it    = 0;
        size_t MaxIt = 15;
        while (it<MaxIt && norm_R>1.0e-8)
        {
            for (size_t i=0; i<dom.Nods.Size();     ++i)
            for (size_t j=0; j<dom.Nods[i]->nDOF(); ++j)
            {
                long eq = dom.Nods[i]->EQ[j];
                if (dom.Nods[i]->pU[j]) R(eq) = 0.0;
            }

            sol.AssembleKA ();
            UMFPACK::Solve (sol.A11, R, dUc);

            U += dUc;
            for (size_t i=0; i<dom.Eles.Size(); ++i) dom.Eles[i]->UpdateState (dUc, &Fint);
            R = Fext-Fint;
            norm_R = Norm(R);

            //cout << "dUc  = " << PrintVector(dUc);
            //cout << "Fext = " << PrintVector(Fext);
            //cout << "Fint = " << PrintVector(Fint);
            //cout << "R    = " << PrintVector(R);
            cout << "[1;32mnorm(R) = " << norm_R << "[0m\n";

            sol.Time += 1.0;
            of << _6_3 << sol.Time << _8s << U(eqx) << _8s << Fint(eqx) << _8s << Fext(eqx) << endl;

            it++;
        }
        if (it>=MaxIt) throw new Fatal("did not converge");
    }

    of.close();

    return 0.0;
}
MECHSYS_CATCH
