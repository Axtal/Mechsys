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

/*  Abbo & Sloan (1996): Figure 1 p1747 
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
    bool   calc_error = true;
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
    String fkey("abbo_sloan_01");
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
    String buf;
    if (!NR && !FE) buf.Printf("_%d_%g",ninc,STOL);
    else            buf.Printf("_%d",   ninc);
    fkey.append(buf);

    // constants
    double a  = 1;
    double b  = a*2;
    int    nx = 5;
    double t  = (b-a)/5;

    // mesh
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 0.,  a, 0.0,
                 0.,  b, 0.0,
                 0.,  b,   t,
                 0.,  a,   t,  -10., 0., -20., -30.);
    blks[0].SetNx (nx);
    blks[0].SetNy (1);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/true);
    if (!quad8)
    {
        mesh.Quad8ToTri6 ();
        mesh.Tri6ToTri15 ();
    }
    mesh.Check    ();
    mesh.WriteMPY (fkey.CStr(), /*tags*/true, /*ids*/true, /*shares*/true);

    // parameters
    double E   = 100000.;
    double nu  = 0.3;
    double sY  = 100.;
    double Hp  = 0.0;
    double c   = 10.0;
    double phi = 30.0;
    double del = 0.3*a*1.0e-3;

    // output vertices
    Array<int> out_verts = (quad8 ? Array<int>(0,6,22) : Array<int>(0,6,22,44,45));

    // props and domain
    double geom = (quad8 ? GEOM("Quad8") : GEOM("Tri15"));
    double nip  = (quad8 ? 4.0           : 16.0);
    Dict prps, mdls;
    prps.Set(-1, "prob geom axs nip", PROB("Equilib"), geom, 1.0, nip);
    //mdls.Set(-1, "name E nu fc sY Hp axs", MODEL("ElastoPlastic"), E, nu, FAILCRIT("VM"), sY, Hp, 1.0);
    mdls.Set(-1, "name E nu fc c phi axs", MODEL("ElastoPlastic"), E, nu, FAILCRIT("MC"), c, phi, 1.0);
    //mdls.Set(-1, "name E nu axs", MODEL("LinElastic"), E, nu, 1.0);
    FEM::Domain dom(mesh, prps, mdls, /*inis*/Dict(), fkey.CStr(), &out_verts);

    // solver
    FEM::Solver sol(dom);
    if (NR) sol.SetScheme("NR");
    if (FE) sol.SetScheme("FE");
    sol.STOL = STOL;
    sol.CorR = CorR;

    // solve
    Dict bcs;
    bcs.Set      (-10, "uy", 0.0);
    bcs.Set      (-20, "uy", 0.0);
    bcs.Set      (-30, "ux", del);
    dom.SetBCs   (bcs);
    sol.Solve    (ninc);
    dom.WriteVTU (fkey.CStr());

    // results
    std::ostringstream disp, stress;
    size_t idx_sx=-1, idx_sy=-1;
    for (size_t i=0; i<dom.AllEKeys.Size(); ++i)
    {
        if (dom.AllEKeys[i]=="sx") idx_sx = i;
        if (dom.AllEKeys[i]=="sy") idx_sy = i;
    }
    for (size_t i=0; i<dom.Nods.Size(); ++i)
    {
        disp   << Util::_20_15 << dom.Nods[i]->U("ux") << " ";
        disp   << Util::_20_15 << dom.Nods[i]->U("uy") << std::endl;
        //stress << Util::_20_15 << dom.NodResults(i,idx_sx)                << " ";
        //stress << Util::_20_15 << dom.NodResults(i,idx_sy)                << std::endl;
    }
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
            u(0+i*2) = dom.Nods[i]->U("ux");
            u(1+i*2) = dom.Nods[i]->U("uy");
            //s(0+i*2) = dom.NodResults(i,idx_sx);
            //s(1+i*2) = dom.NodResults(i,idx_sy);
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
    Vec_t f(out_verts.Size());
    for (size_t i=0; i<out_verts.Size(); ++i) f(i) = dom.Nods[out_verts[i]]->F("fx");
    double q = (quad8 ? sqrt(2.0)*Norm(f)/t : 90.0*Norm(f)/(sqrt(2290.0)*t));
    std::cout << "Error (q/c)    = " << Util::_8s << fabs(q/c-1.0174) << "      q/c = " << Util::_10_6 << (q/c) << "    (1.0174)" << std::endl;

    // end
    return 0;
}
MECHSYS_CATCH
