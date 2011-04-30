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


#ifndef MECHSYS_LBM_DEM_H
#define MECHSYS_LBM_DEM_H

// Mechsys
#include <mechsys/lbm/Lattice.h>

class Disk
{
public:
    //Constructor
    Disk(int Tag, Vec3_t const & X, Vec3_t const & V, Vec3_t const & W, double rho, double R, double dt);

    //Methods
    void ImprintDisk(Lattice & Lat);
    void Translate  (double dt);


    // Data
    int  Tag;              ///< Id of the particle
    Vec3_t X;              ///< position of the particle
    Vec3_t Xb;             ///< position of the particle before
    Vec3_t V;              ///< velocity of the CM
    Vec3_t W;              ///< angular velocity
    Vec3_t F;              ///< Force
    Vec3_t T;              ///< Torque
    double R;              ///< Disk radius
    double M;              ///< mass of the disk
    double I;              ///< inertia moment of the particle
};

inline Disk::Disk(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt)
{
    Tag = TheTag;
    X   = TheX;
    V   = TheV;
    W   = TheW;
    R   = TheR;
    M   = M_PI*R*R*Therho;
    I   = 0.5*M*R*R;
    Xb  = X - dt*V;
    F   = 0.0,0.0,0.0;
    T   = 0.0,0.0,0.0;
}

inline void Disk::ImprintDisk(Lattice & Lat)
{
    for (size_t i=0;i<Lat.Cells.Size();i++)
    {
        Cell  * cell = Lat.Cells[i];
        double x     = Lat.dx*cell->Index(0);
        double y     = Lat.dx*cell->Index(1);
        double z     = Lat.dx*cell->Index(2);
        Vec3_t  C(x,y,z);

        Array<Vec3_t> P(4);

        P[0] = C - 0.5*Lat.dx*OrthoSys::e0 - 0.5*Lat.dx*OrthoSys::e1;
        P[1] = C + 0.5*Lat.dx*OrthoSys::e0 - 0.5*Lat.dx*OrthoSys::e1;
        P[2] = C + 0.5*Lat.dx*OrthoSys::e0 + 0.5*Lat.dx*OrthoSys::e1;
        P[3] = C - 0.5*Lat.dx*OrthoSys::e0 + 0.5*Lat.dx*OrthoSys::e1;

        double len = 0.0;

        for (size_t j=0;j<4;j++)
        {
            Vec3_t D = P[(j+1)%4] - P[j];
            double a = dot(D,D);
            double b = 2*dot(P[j]-X,D);
            double c = dot(P[j]-X,P[j]-X) - R*R;
            if (b*b-4*a*c>0.0)
            {
                double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
                double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
                if (ta>1.0&&tb>1.0) continue;
                if (ta<0.0&&tb<0.0) continue;
                if (ta<0.0) ta = 0.0;
                if (tb>1.0) tb = 1.0;
                len += norm((tb-ta)*D);
            }
        }

        if (len>0.0)
        {
            cell->Gamma   = len/(4*Lat.dx);
            Vec3_t B      = C - X;
            cell->VelP    = V + cross(W,B);
            Vec3_t V;
            double rho = cell->VelDen(V);
            double Bn  = (cell->Gamma*(cell->Tau-0.5))/((1.0-cell->Gamma)+(cell->Tau-0.5));
            for (size_t j=0;j<cell->Nneigh;j++)
            {
                double Feqn    = cell->Feq(j,                   V,rho);
                double Fvpp    = cell->Feq(cell->Op[j],cell->VelP,rho);
                double Fvp     = cell->Feq(j          ,cell->VelP,rho);
                //cell->Omeis[j] = Fvp - cell->F[j] + (1.0 - 1.0/cell->Tau)*(cell->F[j] - Feqn);
                cell->Omeis[j] = cell->F[cell->Op[j]] - Fvpp - (cell->F[j] - Fvp);
                F -= Bn*cell->Omeis[j]*cell->C[j];
            }
        }

    }
}

inline void Disk::Translate(double dt)
{
    //std::cout << F(0) << " " << M << " " << V(0) << std::endl;
    Vec3_t Xa = 2*X - Xb + F*(dt*dt/M);
    Vec3_t tp = Xa - X;
    V         = 0.5*(Xa - Xb)/dt;
    Xb        = X;
    X         = Xa;
}

#endif
