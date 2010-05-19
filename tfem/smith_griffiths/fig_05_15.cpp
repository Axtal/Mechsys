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

/*  Smith & Griffiths (2004): Figure 5.15 p181 
 *  ==========================================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>
#include <mechsys/util/fatal.h>
#include <mechsys/draw.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize    (12/*verts*/, 6/*cells*/);
    mesh.SetVert    ( 0,  0,  0.0,   0.0, 0);
    mesh.SetVert    ( 1,  0,  0.0,  -3.0, 0);
    mesh.SetVert    ( 2,  0,  0.0,  -6.0, 0);
    mesh.SetVert    ( 3,  0,  0.0,  -9.0, 0);
    mesh.SetVert    ( 4,  0,  3.0,   0.0, 0);
    mesh.SetVert    ( 5,  0,  3.0,  -3.0, 0);
    mesh.SetVert    ( 6,  0,  3.0,  -6.0, 0);
    mesh.SetVert    ( 7,  0,  3.0,  -9.0, 0);
    mesh.SetVert    ( 8,  0,  6.0,   0.0, 0);
    mesh.SetVert    ( 9,  0,  6.0,  -3.0, 0);
    mesh.SetVert    (10,  0,  6.0,  -6.0, 0);
    mesh.SetVert    (11,  0,  6.0,  -9.0, 0);
    mesh.SetCell    ( 0, -1, Array<int>(0,1, 5, 4));
    mesh.SetCell    ( 1, -1, Array<int>(1,2, 6, 5));
    mesh.SetCell    ( 2, -1, Array<int>(2,3, 7, 6));
    mesh.SetCell    ( 3, -1, Array<int>(4,5, 9, 8));
    mesh.SetCell    ( 4, -1, Array<int>(5,6,10, 9));
    mesh.SetCell    ( 5, -1, Array<int>(6,7,11,10));
    mesh.SetBryTag  ( 0, 0, -10);
    mesh.SetBryTag  ( 0, 3, -30);
    mesh.SetBryTag  ( 1, 0, -10);
    mesh.SetBryTag  ( 2, 0, -10);
    mesh.SetBryTag  ( 2, 1, -20);
    mesh.SetBryTag  ( 3, 2, -10);
    mesh.SetBryTag  ( 4, 2, -10);
    mesh.SetBryTag  ( 5, 1, -20);
    mesh.SetBryTag  ( 5, 2, -10);
    mesh.GenO2Verts ();
    //mesh.WriteMPY   ("fig_05_15");
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa nip", PROB("Equilib"), GEOM("Quad8"), 1.0, 1.0, 4.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"),  1.0e+6, 0.3,  1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    FEM::Solver sol(dom);
    //sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "ux",     0.0);
    bcs.Set( -20, "ux uy",  0.0, 0.0);
    bcs.Set( -30, "qn",    -1.0);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.Solve ();

    //////////////////////////////////////////////////////////////////////////////////////// Check /////
    
    Table nod_sol;
    nod_sol.Set("ux uy", dom.Nods.Size(),
          0.000000000000000E+00, -5.310671749739340E-06,
          0.000000000000000E+00, -3.243122913926113E-06,
          0.000000000000000E+00, -1.378773528167206E-06,
          0.000000000000000E+00,  0.000000000000000E+00,
         -7.221882229919595E-07, -3.342856970059000E-06,
          3.669995660268174E-07, -2.228571333618901E-06,
          1.912334843273510E-07, -1.114285662470972E-06,
          0.000000000000000E+00,  0.000000000000000E+00,
          0.000000000000000E+00, -1.375042190378655E-06,
          0.000000000000000E+00, -1.214019753311685E-06,
          0.000000000000000E+00, -8.497977967747343E-07,
          0.000000000000000E+00,  0.000000000000000E+00,
          0.000000000000000E+00, -4.288483503711751E-06,
         -4.211153193865630E-07, -5.041231308584792E-06,
          0.000000000000000E+00, -2.217378283206796E-06,
          2.708147610339125E-07, -2.873165860694200E-06,
          0.000000000000000E+00, -6.453528854650004E-07,
          1.370485290061350E-07, -1.298706937331023E-06,
          0.000000000000000E+00,  0.000000000000000E+00,
          3.774202249726477E-07, -2.785714138358060E-06,
         -4.211153193865630E-07, -1.644482664256252E-06,
          2.996313570447692E-07, -1.671428489956403E-06,
          2.708147610339123E-07, -1.583976767620262E-06,
          1.121994376475524E-07, -5.571428285393078E-07,
          1.370485290061355E-07, -9.298643811646859E-07,
          0.000000000000000E+00,  0.000000000000000E+00,
          0.000000000000000E+00, -1.282944773004366E-06,
          0.000000000000000E+00, -1.125478696706005E-06,
          0.000000000000000E+00, -4.689327716136146E-07);

    Table ele_sol;
    ele_sol.Set("sx sy sxy", dom.Eles.Size()*4/*nip*/,
        -5.134800470066937E-01, -1.036826042310187E+00, -2.480701376547251E-02,
        -1.512735293091828E-01, -8.970545301297849E-01,  7.274375440183110E-02,
        -3.507449479150401E-01, -7.711415700600247E-01,  1.505498533653201E-01,
        -1.588965142654356E-01, -6.930632798969868E-01,  2.405335636207552E-01,

        -1.371531628337401E-01, -7.668489560025455E-01,  7.382217423331139E-02,
        -1.445917565088939E-01, -6.451832932738965E-01,  3.542425433673822E-02,
        -1.907499825424094E-01, -5.854932290167336E-01,  1.539849745206280E-01,
        -1.912701951212817E-01, -5.522617535799581E-01,  9.349250571924611E-02,

        -1.700148244179110E-01, -5.995548705775071E-01,  1.848021184176790E-02,
        -2.134985721154105E-01, -5.594348209804142E-01,  9.897323581245114E-03,
        -1.970810461196488E-01, -5.337136609835365E-01,  5.416824394600596E-02,
        -2.135697719939453E-01, -5.250366521780223E-01,  3.341624856075757E-02,

        -7.782647027879343E-02, -2.288584057255886E-01,  1.505498417733656E-01,
        -2.696749039283981E-01, -3.069366958886240E-01,  2.405335610829410E-01,
         8.490862881285841E-02,  3.682606652457288E-02, -2.480700217351958E-02,
        -2.772978888846503E-01, -1.029454456558268E-01,  7.274375693964426E-02,

        -2.378214356514237E-01, -4.145067467688754E-01,  1.539849808197870E-01,
        -2.373012230725507E-01, -4.477382222056506E-01,  9.349250327365824E-02,
        -2.914182553600924E-01, -2.331510197830651E-01,  7.382216793415149E-02,
        -2.839796616849384E-01, -3.548166825117125E-01,  3.542425678232621E-02,

        -2.314903720741837E-01, -4.662863148020727E-01,  5.416824545125712E-02,
        -2.150016461998870E-01, -4.749633236075865E-01,  3.341624577098254E-02,
        -2.585565937759212E-01, -4.004451052081013E-01,  1.848021033651687E-02,
        -2.150728460784219E-01, -4.405651548051949E-01,  9.897326371020091E-03);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy",     1.0e-13, 1.0e-12);
    ele_tol.Set("sx sy sxy", 1.0e-7,  1.0e-7,  1.0e-7);

    // return error flag
    bool err_nods = dom.CheckError   (nod_sol, nod_tol);
    bool err_eles = dom.CheckErrorIP (ele_sol, ele_tol);
    return (err_nods || err_eles);
}
MECHSYS_CATCH
