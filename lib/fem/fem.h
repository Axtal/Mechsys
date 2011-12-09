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

// MechSys
#ifndef INCLUDE_MODELS_ONLY
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/fem/element.h>
#include <mechsys/fem/rod.h>
#include <mechsys/fem/beam.h>
#include <mechsys/fem/flowelem.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/hydromechelem.h>
#include <mechsys/fem/usigelem.h>
#include <mechsys/fem/usigepselem.h>
#include <mechsys/fem/geomelem.h>
#include <mechsys/fem/elems/lin2.h>
#include <mechsys/fem/elems/tri3.h>
#include <mechsys/fem/elems/tri6.h>
#include <mechsys/fem/elems/tri15.h>
#include <mechsys/fem/elems/quad4.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/elems/hex8.h>
#include <mechsys/fem/elems/hex20.h>
//#include <mechsys/fem/elems/tet10.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solvers/stdsolver.h>
#include <mechsys/fem/solvers/uwpsolver.h>
#include <mechsys/linalg/matvec.h>
#endif
#include <mechsys/models/model.h>
#include <mechsys/models/linflow.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/models/nlelastic.h>
#include <mechsys/models/elastoplastic.h>
#include <mechsys/models/problemep.h>
#include <mechsys/models/camclay.h>
#include <mechsys/models/unsatflow.h>
#include <mechsys/models/unconv04.h>
#include <mechsys/models/bbmx.h>
#ifndef INCLUDE_MODELS_ONLY

#define FEM_ALLOC_DOMAIN(msh, mat, inp,   DOM)                                                \
    FEM::Domain DOM(msh, (*inp.Prps), (*inp.Prms), (*inp.Inis), inp.fnkey.CStr(), inp.OutNods);

#define FEM_ALLOC_SOLVER(inp, dom, PTSOL)                \
    String name = (inp.uwp ? "UWPSolver" : "STDSolver"); \
    SDPair flags;                                        \
    inp.SetSolFlags (flags);                             \
    FEM::Solver * PTSOL = AllocSolver (name, dom, flags);

#define FEM_SOLVE(verbose, inp, ptsol, dom)                                                 \
    for (FEM::Domain::Models_t::const_iterator it=dom.Mdls.begin(); it!=dom.Mdls.end(); ++it) inp.SetSUp (it->second); \
    if (inp.ninc<1) inp.ninc = 1;                                                           \
    {                                                                                       \
        for (size_t i=0; i<inp.Stages.Size(); ++i)                                          \
        {                                                                                   \
            dom.SetBCs ((*inp.Stages[i]));                                                  \
            if (verbose) dom.PrintBCs (std::cout);                                          \
            if (inp.vtufile)                                                                \
            {                                                                               \
                if (inp.dyn) ptsol->DynSolve (inp.tf, inp.dt, inp.dtout, inp.fnkey.CStr()); \
                else         ptsol->Solve    (inp.ninc, inp.fnkey.CStr());                  \
            }                                                                               \
            else                                                                            \
            {                                                                               \
                if (inp.dyn) ptsol->DynSolve (inp.tf, inp.dt, inp.dtout);                   \
                else         ptsol->Solve    (inp.ninc);                                    \
            }                                                                               \
        }                                                                                   \
    }

#endif
