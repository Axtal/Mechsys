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

// STD
#include <algorithm>

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
    void FixVelocity() {vf=true,true,true; wf=true,true,true;};
    bool IsFree () {return !vf(0)&&!vf(1)&&!vf(2)&&!wf(0)&&!wf(1)&&!wf(2);}; ///< Ask if the particle has any constrain in its movement

    


    // Data
    int  Tag;              ///< Id of the particle
    Vec3_t X0;             ///< initial position of the particle
    Vec3_t X;              ///< position of the particle
    Vec3_t Xb;             ///< position of the particle before
    Vec3_t V;              ///< velocity of the CM
    Vec3_t W;              ///< angular velocity
    Vec3_t Wb;             ///< angular velocity
    Vec3_t F;              ///< Force
    Vec3_t Ff;             ///< fixed Force
    Vec3_t T;              ///< Torque
    Vec3_t Tf;             ///< fixed Torque
    double R;              ///< Disk radius
    double M;              ///< mass of the disk
    double I;              ///< inertia moment of the particle
    double GT;             ///< dissipation constant for the torque
    double Kn;             ///< Stiffness constant
    bVec3_t vf;            ///< prescribed velocities
    bVec3_t wf;            ///< prescribed angular velocities
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
    Wb  = W;
    Ff  = 0.0,0.0,0.0;
    Tf  = 0.0,0.0,0.0;
    vf  = false,false,false;
    wf  = false,false,false;
    GT  = 1.0e5;
    Kn  = 1.0e3;
}

inline void Disk::ImprintDisk(Lattice & Lat)
{
    for (size_t n=std::max(0.0,double(X(0)-R-Lat.dx)/Lat.dx);n<=std::min(double(Lat.Ndim(0)-1),double(X(0)+R+Lat.dx)/Lat.dx);n++)
    for (size_t m=std::max(0.0,double(X(1)-R-Lat.dx)/Lat.dx);m<=std::min(double(Lat.Ndim(1)-1),double(X(1)+R+Lat.dx)/Lat.dx);m++)
    //for (size_t i=0;i<Lat.Cells.Size();i++)
    {
        Cell  * cell = Lat.GetCell(iVec3_t(n,m,0));
        //Cell  * cell = Lat.Cells[i];
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
            cell->Gamma   = std::max(len/(4*Lat.dx),cell->Gamma);
            if (fabs(cell->Gamma-1.0)<1.0e-12&&fabs(Lat.G)>1.0e-12) continue;
            Vec3_t B      = C - X;
            Vec3_t VelP   = V + cross(W,B);
            Vec3_t V;
            double rho = cell->VelDen(V);
            double Bn  = (cell->Gamma*(cell->Tau-0.5))/((1.0-cell->Gamma)+(cell->Tau-0.5));
            for (size_t j=0;j<cell->Nneigh;j++)
            {
                //double Feqn    = cell->Feq(j,                   V,rho);
                double Fvpp    = cell->Feq(cell->Op[j],VelP,rho);
                double Fvp     = cell->Feq(j          ,VelP,rho);
                //cell->Omeis[j] = Fvp - cell->F[j] + (1.0 - 1.0/cell->Tau)*(cell->F[j] - Feqn);
                cell->Omeis[j] = cell->F[cell->Op[j]] - Fvpp - (cell->F[j] - Fvp);
                Vec3_t Flbm    = -Bn*cell->Omeis[j]*cell->C[j];
                F             += Flbm;
                T             += cross(B,Flbm);
            }
        }
    }
}

inline void Disk::Translate(double dt)
{
    //std::cout << F(0) << " " << M << " " << V(0) << std::endl;
    if (vf(0)) F(0) = 0.0;
    if (vf(1)) F(1) = 0.0;
    if (vf(2)) F(2) = 0.0;
    Vec3_t Xa = 2*X - Xb + F*(dt*dt/M);
    Vec3_t tp = Xa - X;
    V         = 0.5*(Xa - Xb)/dt;
    Xb        = X;
    X         = Xa;

    //if (wf(0)) T(0) = 0.0;
    //if (wf(1)) T(1) = 0.0;
    //if (wf(2)) T(2) = 0.0;
    //Vec3_t Wa = Wb + 2*dt*(T-GT*W)/I;
    //Wb        = W;
    //W         = Wa;

    //std::cout << T(2) << " " << W(2) << " " << Wb(2) << std::endl;
}

#endif
