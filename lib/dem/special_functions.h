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

#ifndef MECHSYS_DEM_SPECIAL_H
#define MECHSYS_DEM_SPECIAL_H

// MechSys
#include <mechsys/dem/distance.h>

// Stl
#include <algorithm>

inline void Erosion(Array<Vec3_t> & V, Array<Array<int> > & E, Array<Array <int> > & F, double R) // Mathematical morphology erosion
{
    if (V.Size()<=3) throw new Fatal("There are no enough vertices to work with");
    Array<Face*> Faces;
    for (size_t i=0; i<F.Size(); i++)
    {
        Array<Vec3_t*> verts(F[i].Size());
        for (size_t j=0; j<F[i].Size(); ++j) verts[j] = (&V[F[i][j]]);
        Faces.Push (new Face(verts));
    }
    Array<Vec3_t> Normal(3);
    Array<Array <int> > VFlist;
    Array<Vec3_t> Vtemp;
    for (size_t i=0; i<F.Size()-2; i++)
    {
        Normal[0] = cross(Faces[i]->Edges[0]->dL,Faces[i]->Edges[1]->dL);
        Normal[0] = Normal[0]/norm(Normal[0]);
        Normal[0] = *Faces[i]->Edges[0]->X0 - R*Normal[0];
        for (size_t j=i+1; j<F.Size()-1; j++) 
        {
            Normal[1] = cross(Faces[j]->Edges[0]->dL,Faces[j]->Edges[1]->dL);
            Normal[1] = Normal[1]/norm(Normal[1]);
            Normal[1] = *Faces[j]->Edges[0]->X0 - R*Normal[1];
            for (size_t k=j+1; k<F.Size(); k++)
            {
                Normal[2] = cross(Faces[k]->Edges[0]->dL,Faces[k]->Edges[1]->dL);
                Normal[2] = Normal[2]/norm(Normal[2]);
                Normal[2] = *Faces[k]->Edges[0]->X0 - R*Normal[2];
                Mat_t M(9,9);
                Vec_t X(9),B(9);
                for (size_t l = 0;l<9;l++)
                {
                    for (size_t m = 0;m<9;m++)
                    {
                        M(l,m) = 0;
                    }
                    X(l) = 0;
                    B(l) = 0;
                }
                M(0,0) = Faces[i]->Edges[0]->dL(0);
                M(0,1) = Faces[i]->Edges[1]->dL(0);
                M(0,6) = -1;
                M(1,0) = Faces[i]->Edges[0]->dL(1);
                M(1,1) = Faces[i]->Edges[1]->dL(1);
                M(1,7) = -1;
                M(2,0) = Faces[i]->Edges[0]->dL(2);
                M(2,1) = Faces[i]->Edges[1]->dL(2);
                M(2,8) = -1;
                M(3,2) = Faces[j]->Edges[0]->dL(0);
                M(3,3) = Faces[j]->Edges[1]->dL(0);
                M(3,6) = -1;
                M(4,2) = Faces[j]->Edges[0]->dL(1);
                M(4,3) = Faces[j]->Edges[1]->dL(1);
                M(4,7) = -1;
                M(5,2) = Faces[j]->Edges[0]->dL(2);
                M(5,3) = Faces[j]->Edges[1]->dL(2);
                M(5,8) = -1;
                M(6,4) = Faces[k]->Edges[0]->dL(0);
                M(6,5) = Faces[k]->Edges[1]->dL(0);
                M(6,6) = -1;
                M(7,4) = Faces[k]->Edges[0]->dL(1);
                M(7,5) = Faces[k]->Edges[1]->dL(1);
                M(7,7) = -1;
                M(8,4) = Faces[k]->Edges[0]->dL(2);
                M(8,5) = Faces[k]->Edges[1]->dL(2);
                M(8,8) = -1;
                B(0) = -Normal[0](0);
                B(1) = -Normal[0](1);
                B(2) = -Normal[0](2);
                B(3) = -Normal[1](0);
                B(4) = -Normal[1](1);
                B(5) = -Normal[1](2);
                B(6) = -Normal[2](0);
                B(7) = -Normal[2](1);
                B(8) = -Normal[2](2);
                try { Sol(M,B,X); }
                catch (Fatal * fatal) { continue; }
                Vec3_t Inter;
                Inter(0) = X(6);
                Inter(1) = X(7);
                Inter(2) = X(8);
                bool inside = true;
                for (size_t l = 0;l<F.Size();l++)
                {
                    Vec3_t ct(0,0,0);
                    for (size_t m = 0;m<Faces[l]->Edges.Size();m++)
                    {
                        ct += *Faces[l]->Edges[m]->X0;
                    }
                    ct /= Faces[l]->Edges.Size();
                    Vec3_t temp = Inter - ct;
                    Vec3_t N = cross(Faces[l]->Edges[0]->dL,Faces[l]->Edges[1]->dL);
                    if (dot(temp,N)>0) inside = false;
                }
                if (inside) 
                {
                    bool belong = true;
                    for (size_t l = 0;l<F.Size();l++)
                    {
                        if ((Distance(Inter,*Faces[l])<R)&&(l!=i)&&(l!=j)&&(l!=k))
                        {
                            belong = false;
                        }
                    }
                    if (belong)
                    {
                        Vtemp.Push(Inter);
                        Array<int> VFlistaux(0);
                        VFlistaux.Push(i);
                        VFlistaux.Push(j);
                        VFlistaux.Push(k);
                        VFlist.Push(VFlistaux);
                    }
                }
            }
        }
    }
    V = Vtemp;
    if (V.Size()<=3) throw new Fatal("The erosion gave too few vertices to build a convex hull, try a smaller erosion parameter");
    Vec3_t cm(0,0,0);
    for (size_t i = 0; i < V.Size(); i++)
    {
        cm += V[i];
    }
    cm /=V.Size();
    E.Resize(0);
    for (size_t i = 0; i < V.Size()-1; i++)
    {
        for (size_t j = i+1; j < V.Size(); j++)
        {
            bool first = true;
            for (size_t k = 0; k < 3; k++)
            {
                if (VFlist[j].Find(VFlist[i][k])!=-1)
                {
                    if (!first)
                    {
                        Array<int> Eaux(2);
                        Eaux[0] = i;
                        Eaux[1] = j;
                        E.Push(Eaux);
                    }
                    first = false;
                }
            }
        }
    }
    Array<Array <int> > Ftemp;
    Array <int> Faux;
    for (size_t i = 0;i < F.Size();i++)
    {
        Faux.Resize(0);
        for (size_t j = 0;j<V.Size();j++)
        {
            if (VFlist[j].Find(i)!=-1)
            {
                Faux.Push(j);
            }
        }
        if(Faux.Size()!=0) Ftemp.Push(Faux);
    }

    for (size_t i=0;i<Ftemp.Size();i++)
    {
        Array<double> Angles(Ftemp[i].Size());
        Vec3_t ct(0,0,0);
        for (size_t j = 0;j<Ftemp[i].Size();j++)
        {
            ct += V[Ftemp[i][j]];
        }
        ct /= Ftemp[i].Size();
        Vec3_t axis = V[Ftemp[i][0]] - ct;
        Vec3_t inward = cm - ct;
        Angles[0]=0.0;
        for (size_t j = 1;j<Ftemp[i].Size();j++)
        {
            Vec3_t t1 = V[Ftemp[i][j]] - ct;
            Angles[j] = acos(dot(axis,t1)/(norm(axis)*norm(t1)));
            Vec3_t t2 = cross(axis,t1);
            if (dot(t2,inward)>0) Angles[j] = 2*M_PI - Angles[j];
        }
        bool swap = false;
        while (!swap)
        {
            swap = true;
            for (size_t j=0;j<Ftemp[i].Size()-1;j++)
            {
                if (Angles[j]>Angles[j+1])
                {
                    double temp1 = Angles[j];
                    Angles[j] = Angles[j+1];
                    Angles[j+1] = temp1;
                    int temp2 = Ftemp[i][j];
                    Ftemp[i][j] = Ftemp[i][j+1];
                    Ftemp[i][j+1] = temp2;
                    swap = false;
                }
            }
        }
    }
    F = Ftemp;
}

inline void PolyhedraMP(Array<Vec3_t> & V, Array<Array <int> > & F, double & vol, Vec3_t & CM, Mat3_t & It) // Calculate mass properties of general polyhedra
{
    vol = 0.0;
    CM = 0.0,0.0,0.0;    
    Array<Face*> Faces;
    for (size_t i=0; i<F.Size(); i++)
    {
        Array<Vec3_t*> verts(F[i].Size());
        for (size_t j=0; j<F[i].Size(); ++j) verts[j] = (&V[F[i][j]]);
        Faces.Push (new Face(verts));
    }
    for (size_t i=0;i<F.Size();i++)
    {
        Vec3_t C;
        Faces[i]->Centroid(C);
        double x0 = C(0);
        double y0 = C(1);
        double z0 = C(2);
        for (size_t j=0;j<Faces[i]->Edges.Size();j++)
        {
            double x1 = (*Faces[i]->Edges[j]->X0)(0);
            double x2 = (*Faces[i]->Edges[(j+1)%(Faces[i]->Edges.Size())]->X0)(0);
            double y1 = (*Faces[i]->Edges[j]->X0)(1);
            double y2 = (*Faces[i]->Edges[(j+1)%(Faces[i]->Edges.Size())]->X0)(1);
            double z1 = (*Faces[i]->Edges[j]->X0)(2);
            double z2 = (*Faces[i]->Edges[(j+1)%(Faces[i]->Edges.Size())]->X0)(2);
            double d0 = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
            double d1 = (z1-z0)*(x2-x0)-(z2-z0)*(x1-x0);
            double d2 = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
            vol += d0*(x0+x1+x2)/6;
            double ixx = d0*(x0*x0+x0*x1+x1*x1+x2*(x0+x1+x2))/12.0;
            double iyy = d1*(y0*y0+y0*y1+y1*y1+y2*(y0+y1+y2))/12.0;
            double izz = d2*(z0*z0+z0*z1+z1*z1+z2*(z0+z1+z2))/12.0;
            CM += 0.5*Vec3_t(ixx,iyy,izz);
        }
    }
    CM/=vol;
    It(0,0) = 0.0;
    It(1,1) = 0.0;
    It(2,2) = 0.0;
    It(1,0) = 0.0;
    It(2,0) = 0.0;
    It(2,1) = 0.0;
    for (size_t i=0;i<Faces.Size();i++)
    {
        Vec3_t C;
        Faces[i]->Centroid(C);
        double x0 = C(0)-CM(0);
        double y0 = C(1)-CM(1);
        double z0 = C(2)-CM(2);
        for (size_t j=0;j<Faces[i]->Edges.Size();j++)
        {
            double x1 = (*Faces[i]->Edges[j]->X0)(0)-CM(0);
            double x2 = (*Faces[i]->Edges[(j+1)%(Faces[i]->Edges.Size())]->X0)(0)-CM(0);
            double y1 = (*Faces[i]->Edges[j]->X0)(1)-CM(1);
            double y2 = (*Faces[i]->Edges[(j+1)%(Faces[i]->Edges.Size())]->X0)(1)-CM(1);
            double z1 = (*Faces[i]->Edges[j]->X0)(2)-CM(2);
            double z2 = (*Faces[i]->Edges[(j+1)%(Faces[i]->Edges.Size())]->X0)(2)-CM(2);
            double d0 = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
            double d1 = (z1-z0)*(x2-x0)-(z2-z0)*(x1-x0);
            double d2 = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
            double ixxx = d0*(x0*x0*x0+x0*x0*x1+x0*x1*x1+x1*x1*x1+x2*(x0*x0+x0*x1+x1*x1+x2*(x0+x1+x2)))/20.0;
            double iyyy = d1*(y0*y0*y0+y0*y0*y1+y0*y1*y1+y1*y1*y1+y2*(y0*y0+y0*y1+y1*y1+y2*(y0+y1+y2)))/20.0;
            double izzz = d2*(z0*z0*z0+z0*z0*z1+z0*z1*z1+z1*z1*z1+z2*(z0*z0+z0*z1+z1*z1+z2*(z0+z1+z2)))/20.0;
            double ixxy = d0*(y0*((x0*x0+x0*x1+x1*x1+x2*(x0+x1+x2))+x0*(x0+x1+x2+x0))+y1*((x0*x0+x0*x1+x1*x1+x2*(x0+x1+x2))+x1*(x0+x1+x2+x1))+y2*(x0*x0+x0*x1+x1*x1+x2*(x0+x1+x2))+x2*(x0+x1+x2+x2))/60.0;
            double iyyz = d1*(z0*((y0*y0+y0*y1+y1*y1+y2*(y0+y1+y2))+y0*(y0+y1+y2+y0))+z1*((y0*y0+y0*y1+y1*y1+y2*(y0+y1+y2))+y1*(y0+y1+y2+y1))+z2*(y0*y0+y0*y1+y1*y1+y2*(y0+y1+y2))+y2*(y0+y1+y2+y2))/60.0;
            double izzx = d2*(x0*((z0*z0+z0*z1+z1*z1+z2*(z0+z1+z2))+z0*(z0+z1+z2+z0))+x1*((z0*z0+z0*z1+z1*z1+z2*(z0+z1+z2))+z1*(z0+z1+z2+z1))+x2*(z0*z0+z0*z1+z1*z1+z2*(z0+z1+z2))+z2*(z0+z1+z2+z2))/60.0;

            It(0,0) += (1.0/3.0)*(iyyy+izzz);
            It(1,1) += (1.0/3.0)*(ixxx+izzz);
            It(2,2) += (1.0/3.0)*(iyyy+izzz);
            It(1,0) -= 0.5*ixxy;
            It(2,0) -= 0.5*izzx;
            It(2,1) -= 0.5*iyyz;
        }
    }
    It(0,1) = It(1,0);
    It(0,2) = It(2,0);
    It(1,2) = It(2,1);
}

inline double ReducedValue(double A, double B) //Reduced Value for two quantities
{
    double result = A*B;
    if (result>0.0) result/=(A+B);
    return result;
}

#endif //MECHSYS_DEM_SPECIAL_H

