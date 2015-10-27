/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
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

// Std lib
#include <math.h>

// MechSys
#include <mechsys/dem/graph.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using namespace DEM;

int main(int argc, char **argv) try
{
    double error = 0.0;
    double tol   = 1.0e-15;

    // test vertex-edge distance
    {
        // vertex and edge
        Vec3_t V(1,1,0),V1(1,1,1),V2(0,0,0);
        Vec3_t *a = &V1,*b = &V2;
        Edge   E(a,b);

        // distance
        Vec3_t Xi,Xf;
        Distance (V,E,Xi,Xf);

        // create edge for drawing
        a = &Xi;
        b = &Xf;
        Edge D(a,b);

        // draw
        std::ofstream of("test_distances_VE.pov",std::ios::out);
        POVHeader (of);
        POVSetCam (of, Vec3_t(4,2,2), Vec3_t(0,0,0));
        POVDrawVert  (V, of, /*R*/0.2, "Green");
        E.Draw    (of,    /*R*/0.1, "Red");
        D.Draw    (of,    /*R*/0.1, "Blue");
        of.close  ();

        std::ofstream of2("test_distances_VE.bpy",std::ios::out);
        BPYHeader   (of2);
        BPYDrawVert (V, of2, /*R*/0.2);
        E.Draw      (of2,    /*R*/0.1, "", true);
        D.Draw      (of2,    /*R*/0.1, "", true);
        of2.close   ();

        // check
        double d = Distance(V,E);
        error += fabs(d - sqrt(2.0/3.0));
    }

    // test edge-edge distance
    {
        // two edges
        Vec3_t a(-1.0,-1.0,1.0),b(1.0,1.0,1.0),c(-1.0,1.0,0.0),d(1.0,-1.0,0.0);
        Edge E1(a,b);
        Edge E2(c,d);

        // distance
        Vec3_t Xi,Xf;
        Distance (E1,E2,Xi,Xf);

        // create edge for drawing
        Edge D(Xi,Xf);

        std::ofstream of2("test_distances_EE.bpy",std::ios::out);
        BPYHeader   (of2);
        E1.Draw     (of2,    /*R*/0.2, "", true);
        E2.Draw     (of2,    /*R*/0.2, "", true);
        D.Draw      (of2,    /*R*/0.1, "", true);
        of2.close   ();

        error += fabs(Distance(E1,E2)-1.0);
    }

    // test vertex-face distance
    {
        // vertex
        Vec3_t V(0.0,0.0,1.0);

        // face
        Array<Vec3_t> C(4); // connectivity
        C[0] = -1.0,-1.0,0.0;
        C[1] = -1.0,1.0,0.0;
        C[2] =  1.0,1.0,0.0;
        C[3] = 1.0,-1.0,0.0;

        Face F(C);

        // distance
        Vec3_t Xi,Xf;
        Distance (V,F,Xi,Xf);

        // create edge for drawing
        Edge D(Xi,Xf);

        // draw
        std::ofstream of("test_distances_VF.bpy",std::ios::out);
        BPYHeader   (of);
        BPYDrawVert (V, of, /*R*/0.2);
        F.Draw      (of,    /*R*/0.1, "",true);
        D.Draw      (of,    /*R*/0.1, "",true);
        of.close    ();

        error += fabs(Distance(V,F)-1.0);
    }

    // test torus-vertex distance
    {
        // vertex
        Vec3_t V(10.0,  10.0, 0.0);

        // torus
        Vec3_t P0( 0.0, 0.0, 0.0);
        Vec3_t P1( 1.0, 0.0, 0.0);
        Vec3_t P2( 0.0, 0.0, 1.0);

        Torus T( P0, P1, P2);

        // distance
        Vec3_t Xi,Xf;
        Distance (V,T,Xi,Xf);

        error += fabs(Distance(V,T)-sqrt(100+81));
    }

    // test cylinder-vertex distance
    {
        // vertex
        Vec3_t V(-10.0, 0.5, -10.0);

        // torus
        Vec3_t P00( 0.0, 0.0, 0.0);
        Vec3_t P01( 1.0, 0.0, 0.0);
        Vec3_t P02( 0.0, 0.0, 1.0);
        Vec3_t P03( 0.0, 0.0,-1.0);

        Torus T0( P00, P01, P02);

        //Vec3_t P10( 0.0, 1.0, 0.0);
        //Vec3_t P11( 0.5, 1.0, 0.0);
        //Vec3_t P12( 0.0, 1.0, 0.5);
        //Vec3_t P13( 0.0, 1.0,-0.5);
        Vec3_t P10( 0.0, 1.0, 0.0);
        Vec3_t P11( 1.0, 1.0, 0.0);
        Vec3_t P12( 0.0, 1.0, 1.0);
        Vec3_t P13( 0.0, 1.0,-1.0);

        Torus T1( P10, P11, P12);

        Cylinder C(T0,T1,P03,P13);

        // distance
        Vec3_t Xi,Xf;
        Distance(V,C,Xi,Xf);

        error += fabs(Distance(V,C)-sqrt(200)+1);
    }

    cout << "error = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    return (error>tol ? 1 : 0);
}
MECHSYS_CATCH
