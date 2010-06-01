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

/*  Abbo & Sloan (1996): Figure 5 p1751 
 *  ===================================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    // input
    bool   quad8      = false;
    int    ninc       = 10;
    bool   calc_error = false;
    double STOL       = 1.0e-7;
    bool   NR         = false;
    bool   FE         = false;
    bool   CorR       = true;
    if (argc>1) quad8      = atoi(argv[1]);
    if (argc>2) ninc       = atoi(argv[2]);
    if (argc>3) calc_error = atoi(argv[3]);
    if (argc>4) STOL       = atof(argv[4]);
    if (argc>5) NR         = atoi(argv[5]);
    if (argc>6) FE         = atoi(argv[6]);
    if (argc>7) CorR       = atoi(argv[7]);

    // file key
    String fkey("abbo_sloan_02");
    if (quad8) 
    {
        cout << "\n[1;35m=================================== Quad8 ======================================[0m\n";
        fkey.append("_quad8");
    }
    else
    {
        cout << "\n[1;35m=================================== Tri15 ======================================[0m\n";
        fkey.append("_tri15");
    }

    // constants
    double B = 2.0;   // footing width
    double b = B/2.0; // footing half width

    // mesh
    Mesh::Generic * mesh;
    if (quad8)
    {
        Array<Mesh::Block> blks(2);
        blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                     0.,  0.0, 0.0,
                     0.,    b, 0.0,
                     0.,    b, 5*B,
                     0.,  0.0, 5*B,  -20., 0., -30., -10.);
        blks[1].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                     0.,    b, 0.0,
                     0.,  5*B, 0.0,
                     0.,  5*B, 5*B,
                     0.,    b, 5*B,  -20., -10., 0., 0.);
        blks[0].SetNx ( 2);
        blks[0].SetNy ( 4, /*Ax*/-1.0, /*NonLin*/true);
        blks[1].SetNx ( 4, /*Ax*/ 1.0, /*NonLin*/true);
        blks[1].SetNy ( 4, /*Ax*/-1.0, /*NonLin*/true);
        Mesh::Structured * msh = new Mesh::Structured (2); // 2D
        msh->Generate (blks,/*O2*/true);
        mesh = msh;
    }
    else
    {
        Table tab;
        tab.Read ("abbo_sloan_fig5.points");
        Mesh::Unstructured * msh = new Mesh::Unstructured (2); // 2D
        msh->Delaunay    (tab("x"), tab("y"), /*Tag*/-1);
        msh->TagHSeg     (/*y*/10.0, /*xMin*/0.0, /*xMax*/b, -30);
        msh->GroundTags  (/*L*/-10, /*R*/-10, /*B*/-20);
        msh->Tri3ToTri6  ();
        msh->Tri6ToTri15 ();
        mesh = msh;
    }
    mesh->Check    ();
    mesh->WriteMPY (fkey.CStr(), /*tags*/true, /*ids*/true, /*shares*/false);

    // parameters
    double c   = 10.0;
    double G   = 4000.0;
    double nu  = 0.3;
    double E   = 2.0*G*(1.0+nu);
    double phi = 30.0;
    double del = 0.16*B;

    // props and domain
    double geom = (quad8 ? GEOM("Quad8") : GEOM("Tri15"));
    double nip  = (quad8 ? 4.0           : 16.0);
    Dict prps, mdls;
    prps.Set(-1, "prob geom psa nip", PROB("Equilib"), geom, 1.0, nip);
    //mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), E, nu, 1.0);
    mdls.Set(-1, "name E nu fc c phi psa", MODEL("ElastoPlastic"), E, nu, FAILCRIT("MC"), c, phi, 1.0);
    FEM::Domain dom((*mesh), prps, mdls, /*inis*/Dict());

    // solver
    FEM::Solver sol(dom);
    if (NR) sol.SetScheme("NR");
    if (FE) sol.SetScheme("FE");
    sol.STOL  = STOL;
    sol.CorR  = CorR;
    sol.MaxIt = 30;

    // output nodes
    Array<int> nods = (quad8 ? Array<int>(12,43,13,44,14) : Array<int>(0,207,60,206,13,240,69,239,12));
    dom.SetOutNods (fkey.CStr(), nods);

    // solve
    Dict bcs;
    bcs.Set      (-10, "ux", 0.0);
    bcs.Set      (-20, "uy", 0.0);
    bcs.Set      (-30, "uy", -del);
    dom.SetBCs   (bcs);
    sol.Solve    (ninc);
    dom.WriteVTU (fkey.CStr());

    // results
    std::ostringstream disp, stress;
    size_t idx_sx=-1, idx_sy=-1;
    for (size_t i=0; i<dom.EleKeys.Size(); ++i)
    {
        if (dom.EleKeys[i]=="sx") idx_sx = i;
        if (dom.EleKeys[i]=="sy") idx_sy = i;
    }
    for (size_t i=0; i<dom.Nods.Size(); ++i)
    {
        disp   << Util::_20_15 << dom.Nods[i]->U[dom.Nods[i]->UMap("ux")] << " ";
        disp   << Util::_20_15 << dom.Nods[i]->U[dom.Nods[i]->UMap("uy")] << std::endl;
        stress << Util::_20_15 << dom.NodResults(i,idx_sx)                << " ";
        stress << Util::_20_15 << dom.NodResults(i,idx_sy)                << std::endl;
    }
    String buf;
    buf.Printf("%s.disp",fkey.CStr());
    std::ofstream of1(buf.CStr(), std::ios::out);
    of1 << disp.str();
    of1.close();
    buf.Printf("%s.stress",fkey.CStr());
    std::ofstream of2(buf.CStr(), std::ios::out);
    of2 << stress.str();
    of2.close();

    // error
    if (calc_error)
    {
        size_t nn = dom.Nods.Size();
        Vec_t u(nn*2), u_ref(nn*2);
        Vec_t s(nn*2), s_ref(nn*2);
        std::ifstream fu((quad8 ? "abbo_sloan_01_quad8_1000000.disp"   : "abbo_sloan_01_tri15_1000000.disp"),   std::ios::in);
        std::ifstream fs((quad8 ? "abbo_sloan_01_quad8_1000000.stress" : "abbo_sloan_01_tri15_1000000.stress"), std::ios::in);
        for (size_t i=0; i<dom.Nods.Size(); ++i)
        {
            u(0+i*2) = dom.Nods[i]->U[dom.Nods[i]->UMap("ux")];
            u(1+i*2) = dom.Nods[i]->U[dom.Nods[i]->UMap("uy")];
            s(0+i*2) = dom.NodResults(i,idx_sx);
            s(1+i*2) = dom.NodResults(i,idx_sy);
            fu >> u_ref(0+i*2) >> u_ref(1+i*2);
            fs >> s_ref(0+i*2) >> s_ref(1+i*2);
        }
        Vec_t u_dif(u-u_ref);
        Vec_t s_dif(s-s_ref);
        double error_u = Norm(u_dif)/Norm(u_ref);
        double error_s = Norm(s_dif)/Norm(s_ref);
        std::cout << "Error (disp)   = " << Util::_8s << error_u << std::endl;
        std::cout << "Error (stress) = " << Util::_8s << error_s << std::endl;
        std::cout << "Error (funb)   = " << Util::_8s << sol.NormR/sol.MaxNormF << std::endl;
    }

    // pressure
    //double L  = b/2.0;
    //size_t nn = (nods.Size()+1)/2;
    //Vec_t f1(nn), f2(nn);
    //for (size_t i=0; i<nn; ++i) f1(i) = dom.Nods[nods[i]]->F[dom.Nods[nods[i]]->FMap("fx")];
    //double q1 = (quad8 ? sqrt(2.0)*Norm(f1)/L : 90.0*Norm(f1)/(sqrt(2290.0)*L));
    //std::cout << "Error (q/c)    = " << Util::_8s << fabs(q1/c-30.1396) << "      q/c = " << Util::_10_6 << (q1/c) << "    (30.1396)" << std::endl;

    // end
    delete mesh;
    return 0;
}
MECHSYS_CATCH
