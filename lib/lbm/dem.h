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

// MechSys
#include <mechsys/lbm/mixture.h>
#include <mechsys/dem/distance.h>

namespace LBM
{
class Disk
{
public:
    //Constructor
    Disk() {};      ///< Default constructor
    Disk(Vec3_t const & x0, Vec3_t const & v0, double R0, double M0, double K0, double dt);


    // Methods;
    double InsideFluid (size_t i, size_t j, double h);
    void StartForce(Vec3_t f = Vec3_t(0.0,0.0,0.0)) {F = f;}
    void DrawDisk(LBM::Lattice & l, double dt);
    void Move(double dt);


    // Data

    Vec3_t x; ///< Position of disk in space
    Vec3_t v; ///< Velocity of the particle
    double R; ///< Radius of the disk
    double M; ///< Mass of the disk
    double K; ///< Disk stiffness
    Vec3_t F; ///< Force over the particle
    Vec3_t xb;///< Position at the previous tme step given by the 


        
};

// Implentation

inline double Disk::InsideFluid(size_t i, size_t j, double h)
{
    Vec3_t pbl(i*h    ,j*h    ,0.0);
    Vec3_t pbr((i+1)*h,j*h    ,0.0);
    Vec3_t ptl(i*h    ,(j+1)*h,0.0);
    Vec3_t ptr((i+1)*h,(j+1)*h,0.0);
    Edge Xmin(pbl,ptl);
    Edge Xmax(pbr,ptr);
    Edge Ymin(pbl,pbr);
    Edge Ymax(ptl,ptr);

    Vec3_t xi;

    Vec3_t xfxmin;
    Distance(x,Xmin,xi,xfxmin);
    double dxmin = norm(xfxmin-xi);
    
    Vec3_t xfxmax;
    Distance(x,Xmax,xi,xfxmax);
    double dxmax = norm(xfxmax-xi);

    Vec3_t xfymin;
    Distance(x,Ymin,xi,xfymin);
    double dymin = norm(xfymin-xi);

    Vec3_t xfymax;
    Distance(x,Ymax,xi,xfymax);
    double dymax = norm(xfymax-xi);

    if ((dxmin>R)&&(dxmax>R)&&(dymin>R)&&(dymax>R))     return -1.0;
    if ((dxmin<=R)&&(dxmax<=R)&&(dymin<=R)&&(dymax<=R)) return  1.0;
    if ((dxmin<=R)&&(dxmax>R)&&(dymin<=R)&&(dymax>R))
    {
        double dx = i*h - x(0);
        double dy = j*h - x(1);
        double lx = sqrt(R*R-dy*dy);
        return (0.5*(lx*sqrt(R*R - lx*lx) + R*R*atan(lx/sqrt(R*R - lx*lx)))-0.5*(dx*sqrt(R*R - dx*dx) + R*R*atan(dx/sqrt(R*R - dx*dx)))-(lx-dx)*dy)/(h*h);
    }
    if ((dxmin<=R)&&(dxmax>R)&&(dymin>R)&&(dymax<=R))
    {
        double dx = i*h - x(0);
        double dy = x(1) - (j+1)*h;
        double lx = sqrt(R*R-dy*dy);
        return (0.5*(lx*sqrt(R*R - lx*lx) + R*R*atan(lx/sqrt(R*R - lx*lx)))-0.5*(dx*sqrt(R*R - dx*dx) + R*R*atan(dx/sqrt(R*R - dx*dx)))-(lx-dx)*dy)/(h*h);
    }
    if ((dxmin>R)&&(dxmax<=R)&&(dymin>R)&&(dymax<=R))
    {
        double dx = x(0) - (i+1)*h;
        double dy = x(1) - (j+1)*h;
        double lx = sqrt(R*R-dy*dy);
        return (0.5*(lx*sqrt(R*R - lx*lx) + R*R*atan(lx/sqrt(R*R - lx*lx)))-0.5*(dx*sqrt(R*R - dx*dx) + R*R*atan(dx/sqrt(R*R - dx*dx)))-(lx-dx)*dy)/(h*h);
    }
    if ((dxmin>R)&&(dxmax<=R)&&(dymin<=R)&&(dymax>R))
    {
        double dx = x(0) - (i+1)*h;
        double dy = j*h - x(1);
        double lx = sqrt(R*R-dy*dy);
        return (0.5*(lx*sqrt(R*R - lx*lx) + R*R*atan(lx/sqrt(R*R - lx*lx)))-0.5*(dx*sqrt(R*R - dx*dx) + R*R*atan(dx/sqrt(R*R - dx*dx)))-(lx-dx)*dy)/(h*h);
    }
    if ((dxmin<=R)&&(dxmax>R)&&(dymin>R)&&(dymax>R))
    {
        double d = dxmin;
        double theta = acos(d/R);
        return (R*(R*theta-R*sin(theta)))/(h*h);
    }
    if ((dxmin>R)&&(dxmax<=R)&&(dymin>R)&&(dymax>R))
    {
        double d = dxmax;
        double theta = acos(d/R);
        return (R*(R*theta-R*sin(theta)))/(h*h);
    }
    if ((dxmin>R)&&(dxmax>R)&&(dymin<=R)&&(dymax>R))
    {
        double d = dymin;
        double theta = acos(d/R);
        return (R*(R*theta-R*sin(theta)))/(h*h);
    }
    if ((dxmin>R)&&(dxmax>R)&&(dymin>R)&&(dymax<=R))
    {
        double d = dymax;
        double theta = acos(d/R);
        return (R*(R*theta-R*sin(theta)))/(h*h);
    }
    if ((dxmin<=R)&&(dxmax>R)&&(dymin<=R)&&(dymax<=R))
    {
        double ls = (j+1)*h;
        double li = j*h;
        return (0.5*(ls*sqrt(R*R - ls*ls) + R*R*atan(ls/sqrt(R*R - ls*ls)))-0.5*(li*sqrt(R*R - li*li) + R*R*atan(li/sqrt(R*R - li*li)))-(i*h-x(0))*h)/(h*h);
    }
    if ((dxmin<=R)&&(dxmax<=R)&&(dymin>R)&&(dymax<=R))
    {
        double ls = (i+1)*h;
        double li = i*h;
        return (0.5*(ls*sqrt(R*R - ls*ls) + R*R*atan(ls/sqrt(R*R - ls*ls)))-0.5*(li*sqrt(R*R - li*li) + R*R*atan(li/sqrt(R*R - li*li)))-(x(1)-(j+1)*h)*h)/(h*h);
    }
    if ((dxmin>R)&&(dxmax<=R)&&(dymin<=R)&&(dymax<=R))
    {
        double ls = (j+1)*h;
        double li = j*h;
        return (0.5*(ls*sqrt(R*R - ls*ls) + R*R*atan(ls/sqrt(R*R - ls*ls)))-0.5*(li*sqrt(R*R - li*li) + R*R*atan(li/sqrt(R*R - li*li)))-(x(0)-(i+1)*h)*h)/(h*h);
    }
    if ((dxmin<=R)&&(dxmax<=R)&&(dymin<=R)&&(dymax>R))
    {
        double ls = (i+1)*h;
        double li = i*h;
        return (0.5*(ls*sqrt(R*R - ls*ls) + R*R*atan(ls/sqrt(R*R - ls*ls)))-0.5*(li*sqrt(R*R - li*li) + R*R*atan(li/sqrt(R*R - li*li)))-((j+1)*h-x(1))*h)/(h*h);
    }
}

inline void CalcForce(Disk * const A,Disk * const B)
{
    double K = 2*(A->K*B->K/(A->K+B->K)); //Equivalent stiffness
    double G = 16;                        //Dissipation constant
    double dist = norm(B->x - A->x);      //Distance between the two spheres
    Vec3_t n = (B->x - A->x)/dist;        //Normal Vector
    double delta = A->R + B->R - dist;    //Overlapping distance
    if (delta>0.0)
    {
        Vec3_t F = K*delta*n;
        Vec3_t vrel = B->v - A->v;
        F -= G*dot(vrel,n)*n;
        A->F -= F;
        B->F += F;
    }
}

Disk::Disk(Vec3_t const & x0, Vec3_t const & v0, double R0, double M0, double K0, double dt)
    : x(x0), v(v0), R(R0), M(M0), K(K0)
{
    xb = x - v*dt;
} 

inline void Disk::DrawDisk(LBM::Lattice & l, double dt)
{
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		// Set solids and set the magnus velocity
		if (pow(x(0)-double(i+0.5),2.0) + pow(x(1)-double(j+0.5),2.0) <= pow(R,2.0)) l.GetCell(i,j)->SetSolid(v(0), v(1));
        //double vr = InsideFluid(i,j,l._h);
        //if (vr > 0.0) l.GetCell(i,j)->SetSolid(v(0), v(1),vr);
	}

	 //Calculate de force in the solid
	for (size_t i=0; i<l.Nx(); ++i)
	for (size_t j=0; j<l.Ny(); ++j)
	{
		LBM::Cell * c = l.GetCell(i,j);
		if (c->IsSolid())
		{
			double rho = 0.0;
			for (size_t k=0; k<l.NNeigh(); k++) rho += c->F(k);
            
            double gamma = c->_gamma;
            double tau = l._tau;
            double beta = gamma*(tau-0.5)/((1-gamma)+(tau-0.5));
            Vec3_t v;
            c->Velocity(v,l._Cs);
			
			for (size_t k=0; k<l.NNeigh(); ++k)
			{
				if (l.GetCell(c->Neigh(k))->IsSolid()) continue;
				size_t op = c->Opp(k);
				double alpha = 6.0*c->W(op)*rho;
				F(0) += (2.0*c->F(op) - alpha*(c->C(op,0)*v(0)+c->C(op,1)*v(1)))/dt*c->C(op,0);
				F(1) += (2.0*c->F(op) - alpha*(c->C(op,0)*v(0)+c->C(op,1)*v(1)))/dt*c->C(op,1);

                //double feqnop = c->EqFun (Cell::OPPOSITE2D[k],v,rho,l._Cs);
                //double feqnb = c->EqFun (k,c->_vel_bc,rho,l._Cs);
                //double gamma = c->_gamma;
                //double beta = gamma*(l._tau-0.5)/((1-gamma)+(tau-0.5));
                //double fm = c->F(Cell::OPPOSITE2D[k])-c->F(k)+feqnb-feqnop;
                //F(0) += l._Cs*l._h*l._h*beta*fm*c->C(k,0);
                //F(1) += l._Cs*l._h*l._h*beta*fm*c->C(k,1);
			}
		}
	}

}

inline void Disk::Move(double dt)
{
    Vec3_t temp,xa;
    xa    = 2*x - xb + F*(dt*dt/M);
    temp  = xa - x;
    v    = 0.5*(xa - xb)/dt;
    xb   = x;
    x    = xa;
}

}
#endif // MECHSYS_LBM_DEM_H
