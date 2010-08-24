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
#include <cstdlib> // for std::rand

// MechSys
#define PARALLEL_DEBUG
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    // to debug: mpirun -np 2 xterm -e gdb --args ./activation01 1

    // input
    bool parallel = false;
    bool nonlin   = false;
    int  nx       = 4;
    int  ny       = 2;
    bool FE       = true;
    bool NR       = false;
    int  nincs    = 1;
    if (argc>1) parallel  = atoi(argv[1]);
    if (argc>2) nonlin    = atoi(argv[2]);
    if (argc>3) nx        = atoi(argv[3]);
    if (argc>4) ny        = atoi(argv[4]);
    if (argc>5) FE        = atoi(argv[5]);
    if (argc>6) NR        = atoi(argv[6]);
    if (argc>7) nincs     = atoi(argv[7]);
    MECHSYS_CATCH_PARALLEL = parallel;

    // mpi
    int  my_id  = -1;
    int  nprocs = 1;
    bool root   = true;
    if (parallel)
    {
#ifdef HAS_MPI
        MECHSYS_MPI_INIT
        my_id  = MPI::COMM_WORLD.Get_rank();
        nprocs = MPI::COMM_WORLD.Get_size();
        if (my_id!=0) root = false;
        if (root) printf("\n%s===================================== Parallel =====================================%s\n",TERM_YELLOW_BLUE,TERM_RST);
#else
        throw new Fatal("main.cpp: this code wasn't compiled with HAS_MPI ==> parallel version is not available");
#endif
    }
    else printf("\n%s====================================== Serial ======================================%s\n",TERM_BLACK_WHITE,TERM_RST);

    // fkey
    String fkey, buf;
    fkey.Printf ("activation01_%d", my_id);

    // mesh
    Array<Mesh::Block> blks(2);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -100.,  0.0, 0.0,
                 -200.,  1.0, 0.0,
                 -300.,  1.0, 0.5,
                 -400.,  0.0, 0.5,  -10.0,-20.0,0.0,-40.0);
    blks[1].Set (/*NDim*/2, /*Tag*/-2, /*NVert*/4,
                 -500.,  0.0, 0.5,
                 -600.,  1.0, 0.5,
                 -700.,  1.0, 1.0,
                 -800.,  0.0, 1.0,  0.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (nx);
    blks[0].SetNy (ny);
    blks[1].SetNx (nx);
    blks[1].SetNy (ny);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.WithInfo = root;
    mesh.Generate (blks,/*O2*/false);
    if (parallel) mesh.PartDomain (nprocs);
    buf = fkey + "_mesh";
    mesh.WriteMPY (buf.CStr());
    mesh.WriteVTU (buf.CStr());

    // initial values
    Dict inis;
    //inis.Set (-1, "geostatic K0 gam y_surf", 1., 0.5, 10.0, 0.5);
    //inis.Set (-2, "sx sy sz", 0., 0., 0.);

    // domain
    FEM::Domain::PARA = parallel;
    Array<int> out_verts(-300,/*justone*/true);
    Dict prps, mdls;
    prps.Set (-1, "prob geom psa active rho", PROB("Equilib"), GEOM("Quad4"), 1., 1., 1.0);
    prps.Set (-2, "prob geom psa active rho", PROB("Equilib"), GEOM("Quad4"), 1., 0., 1.0);
    if (nonlin)
    {
        mdls.Set (-1, "name K0 G0 alp bet psa", MODEL("NLElastic"), 4000.0, 4000.0, 0.4, 0.4, 1.);
        mdls.Set (-2, "name K0 G0 alp bet psa", MODEL("NLElastic"), 4000.0, 4000.0, 0.4, 0.4, 1.);
    }
    else
    {
        mdls.Set (-1, "name E nu psa", MODEL("LinElastic"), 1000.0, 0.2, 1.);
        mdls.Set (-2, "name E nu psa", MODEL("LinElastic"), 1000.0, 0.2, 1.);
    }
    FEM::Domain dom(mesh, prps, mdls, inis, fkey.CStr(), &out_verts);

    // solver
    FEM::Solver sol(dom);
    if (FE) sol.SetScheme("FE");
    if (NR) sol.SetScheme("NR");

    // stage # 1 -------------------------------------
    size_t stg = 1;
    printf("\n%s====================================== Stage # 1 ===================================%s\n",TERM_YELLOW_BLUE,TERM_RST);
    buf.Printf("%s_stage_%d",fkey.CStr(),stg);
    Dict bcs;
    bcs.Set      (-10, "uy", 0.0);
    bcs.Set      (-40, "ux", 0.0);
    bcs.Set      (-1, "gravity", 10.0);
    dom.SetBCs   (bcs);
    sol.Solve    (nincs);
    dom.WriteVTU (buf.CStr());
    dom.PrintResults ("%15.6e", /*onlysummary*/false, /*withelems*/false);

    // check
    Table  nod_sol;
    SDPair nod_tol;
    Array<bool> errors;
    {
        nod_tol.Set ("ux uy Rux Ruy", 1.0e-17, 1.0e-16, 1.0e-14, 1.0e-14);
        nod_sol.Set ("                   ux                       uy                      Rux                     Ruy", /*nrows*/25,
                      0.0,                     0.0,                     3.892689795808790e-02,  6.283246868631101e-01,  //  0
                      1.542029744253440e-04,   0.0,                     0.0,                    1.259944524157890e+00,  //  1
                      3.187008560576240e-04,   0.0,                     0.0,                    1.265961771803640e+00,  //  2
                      5.060750411671041e-04,   0.0,                     0.0,                    1.253272398942170e+00,  //  3
                      7.234217992983710e-04,   0.0,                     0.0,                    5.924966182331980e-01,  //  4
                      0.0,                    -8.854485051446230e-04,   4.744665593691200e-04,  0.0,                    //  5
                      1.502269998821290e-04,  -8.889267220985640e-04,   0.0,                    0.0,                    //  6
                      3.029642044612330e-04,  -8.974323056218250e-04,   0.0,                    0.0,                    //  7
                      4.606176361447270e-04,  -8.980214840474780e-04,   0.0,                    0.0,                    //  8
                      6.108135565641700e-04,  -8.293493102852790e-04,   0.0,                    0.0,                    //  9
                      0.0,                    -1.205592244608900e-03,  -3.940136451745700e-02,  0.0,                    // 10
                      1.491185603563630e-04,  -1.209509893574970e-03,   0.0,                    0.0,                    // 11
                      2.861668803803680e-04,  -1.218572549683000e-03,   0.0,                    0.0,                    // 12
                      3.905365902813930e-04,  -1.209630483826970e-03,   0.0,                    0.0,                    // 13
                      4.549510875732880e-04,  -1.118981901221230e-03,   0.0,                    0.0,                    // 14
                      0.0,                     0.0,                     0.0,                    0.0,                    // 15
                      0.0,                     0.0,                     0.0,                    0.0,                    // 16
                      0.0,                     0.0,                     0.0,                    0.0,                    // 17
                      0.0,                     0.0,                     0.0,                    0.0,                    // 18
                      0.0,                     0.0,                     0.0,                    0.0,                    // 19
                      0.0,                     0.0,                     0.0,                    0.0,                    // 20
                      0.0,                     0.0,                     0.0,                    0.0,                    // 21
                      0.0,                     0.0,                     0.0,                    0.0,                    // 22
                      0.0,                     0.0,                     0.0,                    0.0,                    // 23
                      0.0,                     0.0,                     0.0,                    0.0);                   // 24
    }
    errors.Push (dom.CheckErrorNods (nod_sol, nod_tol));

    for (size_t i=0; i<0; ++i)
    {
        // stage # stg -------------------------------------
        stg++;
        printf("\n%s====================================== Stage # %zd === (activate) ====================%s\n",TERM_YELLOW_BLUE,stg,TERM_RST);
        buf.Printf("%s_stage_%d",fkey.CStr(),stg);
        bcs.clear();
        bcs.Set      (-10, "uy", 0.0);
        bcs.Set      (-40, "ux", 0.0);
        bcs.Set      (-2, "activate gravity", 1., 10.0);
        dom.SetBCs   (bcs);
        sol.Solve    (nincs);
        dom.WriteVTU (buf.CStr());
        //dom.PrintResults ("%15.6e", /*onlysummary*/false, /*withelems*/false);

        // check
        {
            nod_sol.Set ("                   ux                       uy                      Rux                     Ruy", /*nrows*/25,
                          0.0,                     0.0,                     5.173603378432450e-02,  1.274489475357090e+00,  //  0
                          4.447289583704260e-04,   0.0,                     0.0,                    2.546829644806320e+00,  //  1
                          9.047155851342100e-04,   0.0,                     0.0,                    2.531736894158710e+00,  //  2
                          1.389689368991100e-03,   0.0,                     0.0,                    2.471338432409340e+00,  //  3
                          1.891875770854880e-03,   0.0,                     0.0,                    1.175605553268530e+00,  //  4
                          0.0,                    -2.118324813332090e-03,   3.891375844225420e-02,  0.0,                    //  5
                          4.268132233651130e-04,  -2.118520772708360e-03,   0.0,                    0.0,                    //  6
                          8.612182040420320e-04,  -2.112908251925140e-03,   0.0,                    0.0,                    //  7
                          1.307033080993060e-03,  -2.074172644917210e-03,   0.0,                    0.0,                    //  8
                          1.742396791772380e-03,  -1.929039988223290e-03,   0.0,                    0.0,                    //  9
                          0.0,                    -3.656932334878470e-03,  -5.409225759231010e-03,  0.0,                    // 10
                          3.872419712641690e-04,  -3.655372920066210e-03,   0.0,                    0.0,                    // 11
                          7.654335842111140e-04,  -3.640381240641550e-03,   0.0,                    0.0,                    // 12
                          1.116568909128100e-03,  -3.566403128802770e-03,   0.0,                    0.0,                    // 13
                          1.434161521855390e-03,  -3.321460805164210e-03,   0.0,                    0.0,                    // 14
                          0.0,                    -3.364334671438990e-03,  -3.960790237148980e-02,  0.0,                    // 15
                          1.873621738370680e-04,  -3.357632583057290e-03,   0.0,                    0.0,                    // 16
                          3.696140223637040e-04,  -3.327844815131080e-03,   0.0,                    0.0,                    // 17
                          5.415092298132580e-04,  -3.244417792054800e-03,   0.0,                    0.0,                    // 18
                          6.934921977662340e-04,  -3.046758325150510e-03,   0.0,                    0.0,                    // 19
                          0.0,                    -3.692443879127750e-03,  -4.563266409585790e-02,  0.0,                    // 20
                          1.493990589845870e-04,  -3.683830758413770e-03,   0.0,                    0.0,                    // 21
                          2.767392524126690e-04,  -3.647988336203710e-03,   0.0,                    0.0,                    // 22
                          3.650010669630910e-04,  -3.550035458852070e-03,   0.0,                    0.0,                    // 23
                          4.229742939303970e-04,  -3.343847013933140e-03,   0.0,                    0.0);                   // 24
        }
        errors.Push (dom.CheckErrorNods (nod_sol, nod_tol));

        // stage # stg -------------------------------------
        stg++;
        printf("\n%s====================================== Stage # %zd === (deactivate) ==================%s\n",TERM_YELLOW_BLUE,stg,TERM_RST);
        buf.Printf("%s_stage_%d",fkey.CStr(),stg);
        bcs.clear();
        bcs.Set      (-10, "uy", 0.0);
        bcs.Set      (-40, "ux", 0.0);
        bcs.Set      (-2, "deactivate gravity", 1., 10.0);
        dom.SetBCs   (bcs);
        sol.Solve    (nincs);
        dom.WriteVTU (buf.CStr());
        //dom.PrintResults ("%15.6e", /*onlysummary*/false, /*withelems*/false);

        // check
        {
            nod_sol.Set ("                   ux                       uy                      Rux                     Ruy", /*nrows*/25,
                          0.0,                     0.0,                     3.892689795808790e-02,  6.283246868631101e-01,  //  0
                          1.542029744253440e-04,   0.0,                     0.0,                    1.259944524157890e+00,  //  1
                          3.187008560576240e-04,   0.0,                     0.0,                    1.265961771803640e+00,  //  2
                          5.060750411671041e-04,   0.0,                     0.0,                    1.253272398942170e+00,  //  3
                          7.234217992983710e-04,   0.0,                     0.0,                    5.924966182331980e-01,  //  4
                          0.0,                    -8.854485051446230e-04,   4.744665593691200e-04,  0.0,                    //  5
                          1.502269998821290e-04,  -8.889267220985640e-04,   0.0,                    0.0,                    //  6
                          3.029642044612330e-04,  -8.974323056218250e-04,   0.0,                    0.0,                    //  7
                          4.606176361447270e-04,  -8.980214840474780e-04,   0.0,                    0.0,                    //  8
                          6.108135565641700e-04,  -8.293493102852790e-04,   0.0,                    0.0,                    //  9
                          0.0,                    -1.205592244608900e-03,  -3.940136451745700e-02,  0.0,                    // 10
                          1.491185603563630e-04,  -1.209509893574970e-03,   0.0,                    0.0,                    // 11
                          2.861668803803680e-04,  -1.218572549683000e-03,   0.0,                    0.0,                    // 12
                          3.905365902813930e-04,  -1.209630483826970e-03,   0.0,                    0.0,                    // 13
                          4.549510875732880e-04,  -1.118981901221230e-03,   0.0,                    0.0,                    // 14
                          0.0,                     0.0,                     0.0,                    0.0,                    // 15
                          0.0,                     0.0,                     0.0,                    0.0,                    // 16
                          0.0,                     0.0,                     0.0,                    0.0,                    // 17
                          0.0,                     0.0,                     0.0,                    0.0,                    // 18
                          0.0,                     0.0,                     0.0,                    0.0,                    // 19
                          0.0,                     0.0,                     0.0,                    0.0,                    // 20
                          0.0,                     0.0,                     0.0,                    0.0,                    // 21
                          0.0,                     0.0,                     0.0,                    0.0,                    // 22
                          0.0,                     0.0,                     0.0,                    0.0,                    // 23
                          0.0,                     0.0,                     0.0,                    0.0);                   // 24
        }
        errors.Push (dom.CheckErrorNods (nod_sol, nod_tol));
    }

    // end
#ifdef HAS_MPI
    if (parallel) MPI::Finalize();
#endif
    for (size_t i=0; i<errors.Size(); ++i) if (errors[i]) return 1;
    return 0;
}
MECHSYS_CATCH
