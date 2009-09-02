/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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
#include "dem/graph.h"
#include "dem/distance.h"
#include "util/fatal.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    double error = 0.0;
    double tol   = 1.0e-15;

    // test vertex-edge distance
    {
        // vertex and edge
        Vec3_t V(1,1,0);
        Edge   E(Vec3_t(1,1,1),Vec3_t(0,0,0));

        // distance
        Vec3_t Xi,Xf;
        Distance (V,E,Xi,Xf);

        // create edge for drawing
        Edge D(Xi,Xf);

        // draw
        std::ofstream of("test_distances_VE.pov",std::ios::out);
        PovHeader (of);
        PovSetCam (of, Vec3_t(4,2,2), Vec3_t(0,0,0));
        PovDrawVert  (V, of, /*R*/0.2, "Green");
        E.Draw    (of,    /*R*/0.1, "Red");
        D.Draw    (of,    /*R*/0.1, "Blue");
        of.close  ();

        std::ofstream of2("test_distances_VE.py",std::ios::out);
        BlenderHeader   (of2);
        BlenderDrawVert (V, of2, /*R*/0.2);
        E.Draw          (of2,    /*R*/0.1, "", true);
        D.Draw          (of2,    /*R*/0.1, "", true);
        of2.close       ();

        // check
        double d = Distance(V,E);
        error += fabs(d - sqrt(2.0/3.0));
    }

    // test edge-edge distance
    {
        // two edges
        Edge E1(Vec3_t(1.0,1.0,1.0), Vec3_t(0,0,0));
        Edge E2(Vec3_t(0.5,1.0,0.0), Vec3_t(2,0,0));

        // distance
        Vec3_t Xi,Xf;
        Distance (E1,E2,Xi,Xf);

        // create edge for drawing
        Edge D(Xi,Xf);

        /*
        Graph g("test_distances_EE",false);
        Vec3_t p(0,0,10),v(0.5,0.5,0.5);
        g.SetCamera(p,v);
        g.DrawEdge(E1,0.2,"Blue");
        g.DrawEdge(E2,0.2,"Blue");
        g.DrawEdge(D,0.1,"Blue");
        g.Close();
        */

        // check
        //error += adsfasdfad
    }

    // test vertex-face distance
    {
        // vertex
        Vec3_t V(1,1,1);

        // face
        Array<Vec3_t> C(3); // connectivity
        C[0] = 1.0,0.0,0.0;
        C[1] = 0.0,0.0,0.0;
        C[2] = 0.0,1.0,0.0;
        Face F(C);

        // distance
        Vec3_t Xi,Xf;
        Distance (V,F,Xi,Xf);

        // create edge for drawing
        Edge D(Xi,Xf);

        // draw
        std::ofstream of("test_distances_VF.pov",std::ios::out);
        PovHeader   (of);
        PovSetCam   (of, Vec3_t(4,2,2), Vec3_t(0,0,0));
        PovDrawVert (V, of, /*R*/0.2, "Green");
        F.Draw      (of,    /*R*/0.1, "Red");
        D.Draw      (of,    /*R*/0.1, "Blue");
        of.close    ();
    }

    cout << "error = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    return (error>tol ? 1 : 0);
}
MECHSYS_CATCH
