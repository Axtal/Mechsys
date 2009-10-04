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

/*  Owen & Hinton (1980): Example 7.9, p262  *
 *  Finite Elements in Plasticity            *
 *  ======================================== */

// STL
#include <iostream>

// MechSys
#include "mesh/structured.h"
#include "fem/elems/quad8.h"
#include "fem/equilibelem.h"
#include "models/linelastic.h"
#include "models/elastoplastic.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using Util::_4;
using Util::_6_3;
using Util::_8s;

struct DbgDat
{
    size_t        IdxIncOut;
    Array<size_t> IncOut;
    Array<int>    P;
    DbgDat () : IdxIncOut(0) {}
};

void DbgFun (FEM::Solver const & Sol, void * Dat)
{
    // time to output ?
    DbgDat * dat = static_cast<DbgDat*>(Dat);
    size_t idx = dat->IdxIncOut;
    if (idx>=dat->IncOut.Size())   return;
    if (Sol.Inc!=dat->IncOut[idx]) return;

    // header
    String fn;  fn.Printf("owen_hinton_02_P%d.res",dat->P[idx]);
    std::ofstream of(fn.CStr(), std::ios::out);
    of << _4<<"P" << _8s<< "r" << _8s<< "sr" << _8s<< "st" << _8s<< "srt" << "\n";

    // nodes
    /*
    Array<int> nods(7);
    nods = 0, 36, 5, 41, 10, 46, 15;
    for (size_t i=0; i<nods.Size(); ++i)
    {
        FEM::Node const & nod = (*Sol.Dom.Nods[nods[i]]);
        long eqx = nod.EQ[nod.FMap("fx")];
        long eqy = nod.EQ[nod.FMap("fy")];
        double fx = Sol.F(eqx);
        double fy = Sol.F(eqy);
        double fr = sqrt(fx*fx+fy*fy);
        cout << fr << endl;
    }
    */

    // results
    Array<int> eles(4);
    eles = 4, 5, 6, 7;
    for (size_t i=0; i<eles.Size(); ++i)
    {
        FEM::Element const & ele = (*Sol.Dom.Eles[i]);
        Array<SDPair> res;
        ele.GetState (res);
        for (size_t j=0; j<ele.GE->NIP; ++j)
        {
            if (fabs(ele.GE->IPs[j].s)<1.0e-5) // mid point
            {
                // coordinates of IP
                Vec_t X;
                ele.CoordsOfIP (j, X);
                double x  = X(0);
                double y  = X(1);
                double r  = sqrt(x*x+y*y);
                double c  = x/r;
                double s  = y/r;
                double cc = c*c;
                double ss = s*s;
                double cs = c*s;

                // rotation to r-t coordinates
                double sx  = res[j]("sx");
                double sy  = res[j]("sy");
                double sxy = res[j]("sxy");
                double sr  =  cc*sx + ss*sy +  2.0*cs*sxy;
                double st  =  ss*sx + cc*sy -  2.0*cs*sxy;
                double srt = -cs*sx + cs*sy + (cc-ss)*sxy;

                // output
                of << _4<<dat->P[idx] << _8s<< r << _8s<< sr << _8s<< st << _8s<< srt << "\n";
            }
        }
    }
    of.close();
    dat->IdxIncOut++;
}

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    String extra("from pylab import *\nfrom data_handler import *\ndat = read_table('owen_hinton_02_mesh.dat')\nplot(dat['x'],dat['y'],'ro',lw=3)\n");
    Mesh::Structured mesh(/*NDim*/2);
    //mesh.GenQRing (/*O2*/true,/*Nx*/4,/*Ny*/1,/*r*/100.,/*R*/200.,/*Nb*/3,/*Ax*/1.0); // w = 1 + Ax*i
    mesh.GenQRing (/*O2*/true,/*Nx*/0,/*Ny*/1,/*r*/100.,/*R*/200.,/*Nb*/3,/*Ax*/0.0,/*NonLin*/false,/*Wx*/"1.662 2.164 3.034 3.092");
    mesh.WriteMPY ("owen_hinton_02", /*OnlyMesh*/true, extra.CStr());
    mesh.WriteVTU ("owen_hinton_02", /*VolSurfOrBoth*/0);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Quad8"), TRUE);

    // models
    Dict mdls;
    //mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 2.1e+4, 0.3, TRUE);
    mdls.Set(-1, "name E nu sY psa", MODEL("ElastoPlastic"), 2.1e+4, 0.3, 24.0, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(/*NDim*/2, prps, mdls, inis);
    dom.SetMesh    (mesh);
    dom.SetOutNods ("owen_hinton_02", /*NNod*/1, /*IDs*/41);
    dom.SetOutEles ("owen_hinton_02", /*NEle*/1, /*IDs*/4);

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-10, "uy", 0.0)
       .Set(-30, "ux", 0.0)
       .Set(-40, "qn", -20.0);
    dom.SetBCs (bcs);

    // output data
    DbgDat dat;
    dat.IncOut.Resize (4);
    dat.P     .Resize (4);
    dat.IncOut = 3, 5, 6, 8;
    dat.P      = 8, 12, 14, 18;

    // solver
    FEM::Solver sol(dom, &DbgFun, &dat);
    //sol.Scheme = FEM::Solver::NR_t;

    // solve
    sol.Solve (/*NInc*/10);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    // draw elements with IPs
    dom.WriteMPY ("owen_hinton_02_elems");

    return 0;
}
MECHSYS_CATCH
