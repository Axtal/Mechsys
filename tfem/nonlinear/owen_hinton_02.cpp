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

const double DelP = 18.0;
const size_t NInc = 10;

struct DbgDat
{
    Array<int> POut;
    std::ofstream of;
    ~DbgDat () { of.close (); }
     DbgDat ()
    {
         of.open ("owen_hinton_02_n41.res",std::ios::out);
         of<<_6_3<<"Time" << _8s<< "P" <<_8s<<"ur"<< _8s<<"fr_int"<<_8s<<"fr_ext\n"; // radial displacement and forces
         of<<_6_3<<  0    << _8s<<  0  <<_8s<< 0  << _8s<<  0     <<_8s<<  0 << endl;
    }
};

void DbgFun (FEM::Solver const & Sol, void * Dat)
{
    DbgDat * dat = static_cast<DbgDat*>(Dat);

    // current P
    double P = (Sol.Inc+1)*DelP/NInc;

    //////////////////////////////////////////////////////////////////////////////////// Control Node /////
    
    {
        size_t inod = 5;//41;
        FEM::Node const & nod = (*Sol.Dom.Nods[inod]);
        long   eqx    = nod.EQ[nod.FMap("fx")];
        long   eqy    = nod.EQ[nod.FMap("fy")];
        double ux     = Sol.U    (eqx),   uy     = Sol.U    (eqy);
        double fx     = Sol.F    (eqx),   fy     = Sol.F    (eqy);
        double fx_int = Sol.F_int(eqx),   fy_int = Sol.F_int(eqy);
        double ur     = sqrt(ux*ux + uy*uy);
        double fr     = sqrt(fx*fx + fy*fy);
        double fr_int = sqrt(fx_int*fx_int + fy_int*fy_int);
        dat->of << _6_3 << Sol.Time << _8s << P << _8s << ur << _8s << fr_int << _8s << fr << endl;
    }

    /////////////////////////////////////////////////////////////////////////////////////////// Elems /////

    {
        // time to output ?
        if (dat->POut.Find((int)P)<0) return;

        // header
        String fn;  fn.Printf("owen_hinton_02_P%d.res",(int)P);
        std::ofstream of(fn.CStr(), std::ios::out);
        of << _8s<<"P" << _8s<< "r" << _8s<< "sr" << _8s<< "st" << _8s<< "srt" << "\n";

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
                    of << _8s<< P << _8s<< r << _8s<< sr << _8s<< st << _8s<< srt << "\n";
                }
            }
        }
        of.close();
    }
}

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    String extra("from pylab import *\nfrom data_handler import *\ndat = read_table('owen_hinton_02_mesh.dat')\nplot(dat['x'],dat['y'],'ro',lw=3)\n");
    Mesh::Structured mesh(/*NDim*/2);
    //mesh.GenQRing (/*O2*/true,/*Nx*/4,/*Ny*/1,/*r*/100.,/*R*/200.,/*Nb*/3,/*Ax*/1.0); // w = 1 + Ax*i
    mesh.GenQRing (/*O2*/true,/*Nx*/0,/*Ny*/1,/*r*/100.,/*R*/200.,/*Nb*/3,/*Ax*/0.0,/*NonLin*/false,/*Wx*/"1.662 2.164 3.034 3.092");
    mesh.WriteMPY ("owen_hinton_02", /*OnlyMesh*/false);//, extra.CStr());
    mesh.WriteVTU ("owen_hinton_02", /*VolSurfOrBoth*/0);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Quad8"), TRUE);

    // models
    Dict mdls;
    //mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 2.1e+4, 0.3, TRUE);
    mdls.Set(-1, "name E nu fc sY psa", MODEL("ElastoPlastic"), 2.1e+4, 0.3, FAILCRIT("VM"), 24.0, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.SetOutNods ("owen_hinton_02", Array<int>(41, /*JustOne*/true));
    dom.SetOutEles ("owen_hinton_02", Array<int>(4,  /*JustOne*/true));

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-10, "uy", 0.0);
    bcs.Set(-30, "ux", 0.0);
    bcs.Set(-40, "qn", -DelP);
    dom.SetBCs (bcs);

    // output data
    DbgDat dat;
    dat.POut.Resize (4);
    dat.POut = 8, 12, 14, 18;

    // solver
    FEM::Solver sol(dom, &DbgFun, &dat);
    //sol.Scheme = FEM::Solver::NR_t;
    sol.nSS    = 10;
    sol.MaxIt  = 10;

    // solve
    sol.Solve (NInc);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    // draw elements with IPs
    dom.WriteMPY ("owen_hinton_02_elems");

    return 0;
}
MECHSYS_CATCH