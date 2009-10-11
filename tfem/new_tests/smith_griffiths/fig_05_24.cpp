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

/*  Smith & Griffiths (2004): Figure 5.24 p193 
 *  ==========================================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "fem/elems/hex20.h"
#include "fem/equilibelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "util/maps.h"
#include "util/util.h"
#include "util/fatal.h"
#include "draw.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    Mesh::Generic mesh(/*NDim*/3);
    mesh.SetSize    (70/*verts*/, 6/*cells*/);

    mesh.SetVert    ( 0,-100,  0.0 ,  0.0,  0.0);
    mesh.SetVert    ( 1,-200,  0.25,  0.0,  0.0);
    mesh.SetVert    ( 2,-100,  0.5 ,  0.0,  0.0);
    mesh.SetVert    ( 3,   0,  0.0 ,  0.0, -0.5);
    mesh.SetVert    ( 4,   0,  0.5 ,  0.0, -0.5);
    mesh.SetVert    ( 5,   0,  0.0 ,  0.0, -1.0);
    mesh.SetVert    ( 6,   0,  0.25,  0.0, -1.0);
    mesh.SetVert    ( 7,   0,  0.5 ,  0.0, -1.0);
    mesh.SetVert    ( 8,   0,  0.0 ,  0.0, -1.5);
    mesh.SetVert    ( 9,   0,  0.5 ,  0.0, -1.5);
    mesh.SetVert    (10,   0,  0.0 ,  0.0, -2.0);
    mesh.SetVert    (11,   0,  0.25,  0.0, -2.0);
    mesh.SetVert    (12,   0,  0.5 ,  0.0, -2.0);

    mesh.SetVert    (13,-200,  0.0 ,  0.5,  0.0);
    mesh.SetVert    (14,-200,  0.5 ,  0.5,  0.0);
    mesh.SetVert    (15,   0,  0.0 ,  0.5, -1.0);
    mesh.SetVert    (16,   0,  0.5 ,  0.5, -1.0);
    mesh.SetVert    (17,   0,  0.0 ,  0.5, -2.0);
    mesh.SetVert    (18,   0,  0.5 ,  0.5, -2.0);

    mesh.SetVert    (19,-100,  0.0 ,  1.0,  0.0);
    mesh.SetVert    (20,-200,  0.25,  1.0,  0.0);
    mesh.SetVert    (21,-100,  0.5 ,  1.0,  0.0);
    mesh.SetVert    (22,   0,  0.0 ,  1.0, -0.5);
    mesh.SetVert    (23,   0,  0.5 ,  1.0, -0.5);
    mesh.SetVert    (24,   0,  0.0 ,  1.0, -1.0);
    mesh.SetVert    (25,   0,  0.25,  1.0, -1.0);
    mesh.SetVert    (26,   0,  0.5 ,  1.0, -1.0);
    mesh.SetVert    (27,   0,  0.0 ,  1.0, -1.5);
    mesh.SetVert    (28,   0,  0.5 ,  1.0, -1.5);
    mesh.SetVert    (29,   0,  0.0 ,  1.0, -2.0);
    mesh.SetVert    (30,   0,  0.25,  1.0, -2.0);
    mesh.SetVert    (31,   0,  0.5 ,  1.0, -2.0);

    mesh.SetVert    (32,   0,  0.0 ,  1.5,  0.0);
    mesh.SetVert    (33,   0,  0.5 ,  1.5,  0.0);
    mesh.SetVert    (34,   0,  0.0 ,  1.5, -1.0);
    mesh.SetVert    (35,   0,  0.5 ,  1.5, -1.0);
    mesh.SetVert    (36,   0,  0.0 ,  1.5, -2.0);
    mesh.SetVert    (37,   0,  0.5 ,  1.5, -2.0);

    mesh.SetVert    (38,   0,  0.0 ,  2.0,  0.0);
    mesh.SetVert    (39,   0,  0.25,  2.0,  0.0);
    mesh.SetVert    (40,   0,  0.5 ,  2.0,  0.0);
    mesh.SetVert    (41,   0,  0.0 ,  2.0, -0.5);
    mesh.SetVert    (42,   0,  0.5 ,  2.0, -0.5);
    mesh.SetVert    (43,   0,  0.0 ,  2.0, -1.0);
    mesh.SetVert    (44,   0,  0.25,  2.0, -1.0);
    mesh.SetVert    (45,   0,  0.5 ,  2.0, -1.0);
    mesh.SetVert    (46,   0,  0.0 ,  2.0, -1.5);
    mesh.SetVert    (47,   0,  0.5 ,  2.0, -1.5);
    mesh.SetVert    (48,   0,  0.0 ,  2.0, -2.0);
    mesh.SetVert    (49,   0,  0.25,  2.0, -2.0);
    mesh.SetVert    (50,   0,  0.5 ,  2.0, -2.0);

    mesh.SetVert    (51,   0,  0.0 ,  2.5,  0.0);
    mesh.SetVert    (52,   0,  0.5 ,  2.5,  0.0);
    mesh.SetVert    (53,   0,  0.0 ,  2.5, -1.0);
    mesh.SetVert    (54,   0,  0.5 ,  2.5, -1.0);
    mesh.SetVert    (55,   0,  0.0 ,  2.5, -2.0);
    mesh.SetVert    (56,   0,  0.5 ,  2.5, -2.0);

    mesh.SetVert    (57,   0,  0.0 ,  3.0,  0.0);
    mesh.SetVert    (58,   0,  0.25,  3.0,  0.0);
    mesh.SetVert    (59,   0,  0.5 ,  3.0,  0.0);
    mesh.SetVert    (60,   0,  0.0 ,  3.0, -0.5);
    mesh.SetVert    (61,   0,  0.5 ,  3.0, -0.5);
    mesh.SetVert    (62,   0,  0.0 ,  3.0, -1.0);
    mesh.SetVert    (63,   0,  0.25,  3.0, -1.0);
    mesh.SetVert    (64,   0,  0.5 ,  3.0, -1.0);
    mesh.SetVert    (65,   0,  0.0 ,  3.0, -1.5);
    mesh.SetVert    (66,   0,  0.5 ,  3.0, -1.5);
    mesh.SetVert    (67,   0,  0.0 ,  3.0, -2.0);
    mesh.SetVert    (68,   0,  0.25,  3.0, -2.0);
    mesh.SetVert    (69,   0,  0.5 ,  3.0, -2.0);

    mesh.SetCell    ( 0, -1, /*NVerts*/20, 5,7,26,24,0,2,21,19,6,16,25,15,1,14,20,13,3,4,23,22);
    mesh.SetCell    ( 1, -2, /*NVerts*/20, 10,12,31,29,5,7,26,24,11,18,30,17,6,16,25,15,8,9,28,27);
    mesh.SetCell    ( 2, -1, /*NVerts*/20, 24,26,45,43,19,21,40,38,25,35,44,34,20,33,39,32,22,23,42,41);
    mesh.SetCell    ( 3, -2, /*NVerts*/20, 29,31,50,48,24,26,45,43,30,37,49,36,25,35,44,34,27,28,47,46);
    mesh.SetCell    ( 4, -1, /*NVerts*/20, 43,45,64,62,38,40,59,57,44,54,63,53,39,52,58,51,41,42,61,60);
    mesh.SetCell    ( 5, -2, /*NVerts*/20, 48,50,69,67,43,45,64,62,49,56,68,55,44,54,63,53,46,47,66,65);

    mesh.SetBryTag  (0, 0, -10);
    mesh.SetBryTag  (0, 2, -20);
    //mesh.SetBryTag  (0, 5, -40);
    mesh.SetBryTag  (1, 0, -10);
    mesh.SetBryTag  (1, 2, -20);
    mesh.SetBryTag  (1, 4, -30);

    mesh.SetBryTag  (2, 0, -10);
    mesh.SetBryTag  (3, 0, -10);
    mesh.SetBryTag  (3, 4, -30);

    mesh.SetBryTag  (4, 0, -10);
    mesh.SetBryTag  (5, 0, -10);
    mesh.SetBryTag  (5, 4, -30);

    mesh.WriteVTU   ("fig_05_24");
    //cout << mesh << endl;
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active nip", PROB("Equilib"), GEOM("Hex20"), TRUE, 8.0);
    prps.Set(-2, "prob geom active nip", PROB("Equilib"), GEOM("Hex20"), TRUE, 8.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu", MODEL("LinElastic"),  100.0, 0.3,  TRUE);
    mdls.Set(-2, "name E nu", MODEL("LinElastic"),   50.0, 0.3,  TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);
    inis.Set(-2, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    FEM::Solver sol(dom);
    sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "fz",  0.0417);
    bcs.Set(-200, "fz", -0.1667);
    //bcs.Set( -40, "qy", -1.0);
    bcs.Set( -10, "ux",       0.0);
    bcs.Set( -20, "uy",       0.0);
    bcs.Set( -30, "ux uy uz", 0.0,0.0,0.0);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.Solve ();

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    //dom.PrintResults (cout, Util::_6s);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////
    
    Table nod_sol;
    nod_sol.Set("ux uy uz", dom.Nods.Size(),
                0.000000000000000E+00,    0.000000000000000E+00,   -2.246347084551239E-02,
                1.583533151292792E-03,    0.000000000000000E+00,   -2.255393943959004E-02,
                3.220086429625668E-03,    0.000000000000000E+00,   -2.333080568038476E-02,
                0.000000000000000E+00,    0.000000000000000E+00,   -1.849043628520537E-02,
                1.543646337680335E-03,    0.000000000000000E+00,   -1.883760248896416E-02,
                0.000000000000000E+00,    0.000000000000000E+00,   -1.442659959332960E-02,
                7.581181943686352E-04,    0.000000000000000E+00,   -1.435284871221009E-02,
                1.510770238131118E-03,    0.000000000000000E+00,   -1.410762498329098E-02,
                0.000000000000000E+00,    0.000000000000000E+00,   -6.163718905594008E-03,
                2.791740187814714E-03,    0.000000000000000E+00,   -6.429638686950788E-03,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,   -2.636620496523055E-03,   -2.091497948288540E-02,
                2.681677651464966E-03,   -2.351684531012892E-03,   -2.156679429586825E-02,
                0.000000000000000E+00,    1.920886717070353E-03,   -1.258000900042771E-02,
                1.401580418768850E-03,    2.013340247449022E-03,   -1.228974275154635E-02,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,   -4.404981175789160E-03,   -1.365965203530571E-02,
                6.192483795629938E-04,   -4.364683902429320E-03,   -1.368261677877242E-02,
                1.280307716613955E-03,   -3.937219777613099E-03,   -1.386506269262310E-02,
                0.000000000000000E+00,    1.078170867081458E-03,   -1.150366327069372E-02,
                8.754415356300463E-04,    1.490221892346032E-03,   -1.152005921380261E-02,
                0.000000000000000E+00,    2.958295463657756E-03,   -9.172143775217564E-03,
                5.427998204887557E-04,    2.965766844941001E-03,   -9.093537641138223E-03,
                1.086759767227250E-03,    3.083941056730035E-03,   -8.864376938531468E-03,
                0.000000000000000E+00,    2.375669564035335E-03,   -4.184288525644190E-03,
                2.000002573276872E-03,    2.632279040385318E-03,   -4.380320869217215E-03,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,   -3.580016648413835E-03,   -5.979374273810778E-03,
               -7.434638380963998E-05,   -3.321158639873764E-03,   -5.812858913992063E-03,
                0.000000000000000E+00,    3.150654743473781E-03,   -5.382993118730838E-03,
                7.612000232719073E-04,    3.312568168671743E-03,   -5.104029499792067E-03,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,   -2.120948625021810E-03,   -2.575633224597299E-03,
               -1.816246042762323E-04,   -2.066444640224186E-03,   -2.541068218277476E-03,
               -3.530266971970211E-04,   -2.120321795218551E-03,   -2.337790153492455E-03,
                0.000000000000000E+00,    2.687182058805452E-04,   -2.403797241083586E-03,
                1.000479849273670E-04,    3.146371148261739E-04,   -2.217984973952097E-03,
                0.000000000000000E+00,    2.651365011444304E-03,   -2.133133596422494E-03,
                2.082137603898848E-04,    2.668818739233793E-03,   -2.080360040782984E-03,
                4.197890653658663E-04,    2.802415372347518E-03,   -1.956535480074450E-03,
                0.000000000000000E+00,    2.293415103238367E-03,   -1.130490294385088E-03,
                6.751997500436722E-04,    2.510903712920318E-03,   -1.197902689561233E-03,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,   -1.712291682872705E-03,   -5.376254787559190E-04,
               -9.194692826417251E-05,   -1.719403339987016E-03,   -4.759107524728894E-04,
                0.000000000000000E+00,    2.044723042196544E-03,   -3.706552779963829E-04,
                8.501171983670954E-05,    2.179727451296569E-03,   -2.926182413887153E-04,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,   -1.707806869074250E-03,    1.459455668335570E-03,
               -3.549251809384875E-05,   -1.701950634641145E-03,    1.490896328624764E-03,
               -7.239581275430445E-05,   -1.670036441986721E-03,    1.451953228172754E-03,
                0.000000000000000E+00,    3.922334156895800E-04,    1.136691564328168E-03,
               -1.456178432000260E-04,    4.848115407724173E-04,    1.138450128357737E-03,
                0.000000000000000E+00,    1.667966295760209E-03,    5.260752671179727E-04,
               -9.978009212936460E-05,    1.700229972915349E-03,    5.388904225798224E-04,
               -2.111305284838336E-04,    1.825746337701156E-03,    5.433207525796582E-04,
                0.000000000000000E+00,    1.572057618883754E-03,   -1.028379544719576E-04,
               -7.448353813973922E-05,    1.715763820799789E-03,   -5.845849709250573E-05,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00,
                0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00);

    Table ele_sol;
    ele_sol.Set("sx sy sz sxy syz szx", dom.Eles.Size(),
               -2.671595025408352E-02, -1.646894047557757E-01, -9.088250520062361E-01,  6.145397987621120E-03,  9.597217670407908E-02,  4.352140418857522E-03,
                3.985245032513518E-02, -5.315753362091763E-02, -6.297606614921736E-01, -2.140099351738117E-03,  7.613541237101888E-02,  4.169425308898757E-03,
               -2.481756100845436E-02, -1.260037128751898E-01, -1.051975684912099E-01,  4.840356098789624E-03,  9.398512977721066E-02, -2.813635910514869E-03,
                2.476569870752665E-02, -8.239979525906903E-02, -2.822358640984907E-01, -3.178619845356539E-03,  1.214016176622576E-01,  2.938576876913932E-03,
                3.766727420797759E-03,  8.619087264119318E-03, -1.468896913963496E-02, -9.006024512952835E-04, -8.732613394225753E-03,  8.652315391720999E-04,
                7.406753152614834E-03, -4.389804673637764E-02, -2.831255330484390E-02, -5.639356271726088E-04,  6.083237639099075E-02,  5.077611027338742E-05);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy uz", 1.0e-9, 1.0e-9, 1.0e-8);
    ele_tol.Set("sx sy sz  sxy syz szx", 1.0e-8,1.0e-8,1.0e-8, 1.0e-9,1.0e-8,1.0e-8);

    // return error flag
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH