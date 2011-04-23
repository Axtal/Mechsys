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

#ifndef MECHSYS_DEM_INTERACTON_H
#define MECHSYS_DEM_INTERACTON_H

// Std Lib
#include <math.h>
#include <map>
#include <vector>
#include <utility>

// MechSys
#include <mechsys/dem/particle.h>

// typedefs
typedef std::map<std::pair<int,int>,Vec3_t> FrictionMap_t;
typedef Array<pair<int,int> > ListContacts_t;

class Interacton   //General class for interactons
{
public:

    // Constructor and destructor
    Interacton () {}           ///< Default constructor
    virtual ~Interacton () {}  ///< Destructor

    // Methods
    virtual bool UpdateContacts   (double alpha) =0;    ///< Update contacts by verlet algorithm
    virtual void CalcForce        (double dt = 0.0) =0; ///< Calculates the contact force between particles
    virtual void UpdateParameters () =0;                ///< Update the parameters

    // Data
    Particle     * P1;        ///< First particle
    Particle     * P2;        ///< Second particle
    size_t         I1;        ///< Index of the first particle
    size_t         I2;        ///< Index of the second particle
};

class CInteracton: public Interacton // Interacton for collision
{
public:
    // Constructor
    CInteracton (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles
    CInteracton () {};

    // Methods
    virtual bool UpdateContacts (double alpha);    ///< Update contacts by verlet algorithm
    virtual void CalcForce      (double dt = 0.0); ///< Calculates the contact force between particles
    virtual void UpdateParameters ();              ///< Update the parameters

    // Data
    double         Kn;        ///< Normal stiffness
    double         Kt;        ///< Tengential stiffness
    double         Gn;        ///< Normal viscous coefficient
    double         Gt;        ///< Tangential viscous coefficient
    double         Mu;        ///< Microscpic coefficient of friction
    double         Epot;      ///< Potential elastic energy
    double         dEvis;     ///< Energy dissipated in viscosity at time step
    double         dEfric;    ///< Energy dissipated by friction at time step
    size_t         Nc;        ///< Number of contacts
    size_t         Nsc;       ///< Number of sliding contacts
    size_t         Nr;        ///< Number of rolling contacts (only for spheres)
    Vec3_t         Fn;        ///< Normal force between elements
    Vec3_t         Fnet;      ///< Net normal force
    Vec3_t         Ftnet;     ///< Net tangential force
    Vec3_t         Xc;        ///< Net Position of the contact
    Mat3_t         B;         ///< Branch tensor for the study of isotropy
    ListContacts_t Lee;       ///< List of edge-edge contacts 
    ListContacts_t Lvf;       ///< List of vertex-face contacts 
    ListContacts_t Lfv;       ///< List of face-vertex contacts
    ListContacts_t Lvt;       ///< List of vertex-torus contacts
    ListContacts_t Ltv;       ///< List of torus-vertex contacts
    FrictionMap_t  Fdee;      ///< Static friction displacement for pair of edges
    FrictionMap_t  Fdvf;      ///< Static friction displacement for pair of vertex-face
    FrictionMap_t  Fdfv;      ///< Static friction displacement for pair of face-vertex
protected:
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, ListContacts_t & L, double dt);
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_contacts        (FeatureA_T & A, FeatureB_T & B, ListContacts_t & L, double alpha);
};

class CInteractonSphere: public CInteracton // Collision interacton for spheres
{
public:
    // Methods 
    CInteractonSphere (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles
    bool UpdateContacts (double alpha);                 ///< Update contacts by verlet algorithm
    void CalcForce (double dt = 0.0);                   ///< Calculates the contact force between particles
    void UpdateParameters ();                           ///< Update the parameters

    // Data
    FrictionMap_t  Fdvv;                                ///< Static Friction displacement for the vertex vertex pair
    ListContacts_t Lvv;                                 ///< List of vertices
    Vec3_t         Fdr;                                 ///< Rolling displacement 
    double         beta;                                ///< Rolling stiffness coefficient
    double         eta;                                 ///< Plastic moment coefficient

protected:
    void _update_rolling_resistance(double dt);               ///< Calculates the rolling resistance torque
};

class BInteracton: public Interacton // Interacton for cohesion
{
public:
    BInteracton (Particle * Pt1, Particle * Pt2, size_t Fi1, size_t Fi2); ///< Constructor requires pointers to both particles and alos the indixes of the size they share

    // Methods
    bool UpdateContacts (double alpha);    ///< Update contacts by verlet algorithm
    void CalcForce      (double dt = 0.0); ///< Calculates the contact force between particles
    void UpdateParameters ();              ///< Update the parameters in case they change

    // Data
    size_t F1;                             ///< Index of the shared face for particle 1
    size_t F2;                             ///< Index of the shared face for particle 2
    double Area;                           ///< Area of the shared side
    double Bn;                             ///< Elastic normal constant for the cohesion
    double Bt;                             ///< Elastic tangential constant for the cohesion
    double Bm;                             ///< Elastic for the cohesion torque
    double Gn;                             ///< Dissipative normal constant
    double Gt;                             ///< Dissipative tangential constant
    double L0;                             ///< Equilibrium distance
    double An;                             ///< Angular displacement
    double eps;                            ///< Maximun strain before fracture
    bool valid;                            ///< Check if the bound has not been broken
    double s1,t1;                          ///< Planar coordinates for face F1
    double s2,t2;                          ///< Planar coordinates for face F2
    Vec3_t Fnet;                           ///< Net force excerted by the interacton

};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// Collision interacton

inline CInteracton::CInteracton (Particle * Pt1, Particle * Pt2)
{
    P1              = Pt1;
    P2              = Pt2;
    I1              = P1->Index;
    I2              = P2->Index;
    double r1       = pow(P1->Props.V,1.0/3.0);
    double r2       = pow(P2->Props.V,1.0/3.0);
    Kn              = (r1+r2)*ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
    Kt              = (r1+r2)*ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
    Gn              = 2*ReducedValue(Pt1->Props.Gn,Pt2->Props.Gn)*ReducedValue(Pt1->Props.m,Pt2->Props.m);
    Gt              = 2*ReducedValue(Pt1->Props.Gt,Pt2->Props.Gt)*ReducedValue(Pt1->Props.m,Pt2->Props.m);
    Mu              = 2*ReducedValue(Pt1->Props.Mu,Pt2->Props.Mu);
    CalcForce(0.0);
}

inline bool CInteracton::UpdateContacts (double alpha)
{
    Lee.Resize(0);
    Lvf.Resize(0);
    Lfv.Resize(0);
    Lvt.Resize(0);
    Ltv.Resize(0);
    if (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*alpha)
    {
        _update_contacts (P1->Edges,P2->Edges,Lee,alpha);
        _update_contacts (P1->Verts,P2->Faces,Lvf,alpha);
        _update_contacts (P1->Faces,P2->Verts,Lfv,alpha);
        _update_contacts (P1->Verts,P2->Tori ,Lvt,alpha);
        _update_contacts (P1->Tori ,P2->Verts,Ltv,alpha);
        if (Lee.Size()>0||Lvf.Size()>0||Lfv.Size()>0||Ltv.Size()>0||Lvt.Size()>0) return true;
        else return false;
    }
    else return false;
}

inline void CInteracton::CalcForce (double dt)
{
    Epot   = 0.0;
    dEvis  = 0.0;
    dEfric = 0.0;
    Nc     = 0.0;
    Nsc    = 0.0;
    Nr     = 0.0;
    Fnet   = 0.0;
    Ftnet  = 0.0;
    Xc     = 0.0;
    _update_disp_calc_force (P1->Edges,P2->Edges,Fdee,Lee,dt);
    _update_disp_calc_force (P1->Verts,P2->Faces,Fdvf,Lvf,dt);
    _update_disp_calc_force (P1->Faces,P2->Verts,Fdfv,Lfv,dt);
    _update_disp_calc_force (P1->Verts,P2->Tori ,Fdvf,Lvt,dt);
    _update_disp_calc_force (P1->Tori ,P2->Verts,Fdfv,Ltv,dt);

    //If there is at least a contact, increase the coordination number of the particles
    if (Nc>0) 
    {
        P1->Cn++;
        P2->Cn++;
    }
}

inline void CInteracton::UpdateParameters ()
{
    Kn              = 2*ReducedValue(P1->Props.Kn,P2->Props.Kn);
    Kt              = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn)*ReducedValue(P1->Props.m,P2->Props.m);
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt)*ReducedValue(P1->Props.m,P2->Props.m);
    Mu              = 2*ReducedValue(P1->Props.Mu,P2->Props.Mu);
}

template<typename FeatureA_T, typename FeatureB_T>
inline void CInteracton::_update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, ListContacts_t & L, double dt)
{
    // update
    for (size_t k=0; k<L.Size(); ++k)
    {
        size_t i = L[k].first;
        size_t j = L[k].second;
        Vec3_t xi, xf;
        Distance ((*A[i]), (*B[j]), xi, xf);
        double dist  = norm(xf-xi);
        double delta = P1->Props.R + P2->Props.R - dist;
        if (delta>0)
        {
            if (delta > 0.8*(P1->Props.R+P2->Props.R)) throw new Fatal("Interacton::_update_disp_calc_force: Maximun overlap detected between particles %d(%d) and %d(%d)",P1->Index,P1->Tag,P2->Index,P2->Tag);

            // Count a contact
            Nc++;

            // update force
            Vec3_t n = (xf-xi)/dist;
            Vec3_t x = xi+n*((P1->Props.R*P1->Props.R-P2->Props.R*P2->Props.R+dist*dist)/(2*dist));
            Xc += x;
            Vec3_t t1,t2,x1,x2;
            Rotation(P1->w,P1->Q,t1);
            Rotation(P2->w,P2->Q,t2);
            x1 = x - P1->x;
            x2 = x - P2->x;
            Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
            Vec3_t vt = vrel - dot(n,vrel)*n;
            Fn = Kn*delta*n;
            Fnet += Fn;
            pair<int,int> p;
            p = make_pair(i,j);
            if (FMap.count(p)==0) FMap[p] = 0.0,0.0,0.0;
            FMap[p] += vt*dt;
            FMap[p] -= dot(FMap[p],n)*n;
            Vec3_t tan = FMap[p];
            if (norm(tan)>0.0) tan/=norm(tan);
            if (norm(FMap[p])>Mu*norm(Fn)/Kt)
            {
                // Count a sliding contact
                Nsc++;
                FMap[p] = Mu*norm(Fn)/Kt*tan;
                dEfric += Kt*dot(FMap[p],vt)*dt;
            }
            Ftnet += Kt*FMap[p];
            Vec3_t F = Fn + Kt*FMap[p] + Gn*dot(n,vrel)*n + Gt*vt;
            if (dot(F,n)<0) F-=dot(F,n)*n;
#ifdef USE_THREAD
            std::lock_guard<std::mutex> lk1(P1->mtex);
            std::lock_guard<std::mutex> lk2(P2->mtex);
#endif
            P1->F += -F;
            P2->F +=  F;
            dEvis += (Gn*dot(vrel-vt,vrel-vt)+Gt*dot(vt,vt))*dt;
            // torque
            Vec3_t T, Tt;
            Tt = cross (x1,F);
            Quaternion_t q;
            Conjugate (P1->Q,q);
            Rotation  (Tt,q,T);
            P1->T -= T;
            Tt = cross (x2,F);
            Conjugate (P2->Q,q);
            Rotation  (Tt,q,T);
            P2->T += T;
            //Transfering the branch vector information
            Vec3_t nor = n;
            for (size_t m=0;m<3;m++)
            {
                for (size_t n=0;n<3;n++)
                {
                    P1->M(m,n)  -= F(m)*x1(n);
                    P2->M(m,n)  += F(m)*x2(n);
                    P1->B(m,n)  += nor(m)*nor(n);
                    P2->B(m,n)  += nor(m)*nor(n);
                    this->B(m,n) = nor(m)*nor(n);
                }
            }

            // potential energy
            Epot += 0.5*Kn*delta*delta+0.5*Kt*dot(FMap[p],FMap[p]);
        }
    }
}

template<typename FeatureA_T, typename FeatureB_T>
inline void CInteracton::_update_contacts (FeatureA_T & A, FeatureB_T & B, ListContacts_t & L, double alpha)
{
    for (size_t i=0; i<A.Size(); ++i)
    for (size_t j=0; j<B.Size(); ++j)
    {
        if (Distance ((*A[i]), (*B[j]))<=P1->Props.R+P2->Props.R+2*alpha)
        {
            pair<int,int> p;
            p = make_pair(i,j);
            L.Push(p);
        }
    }

}

//Collision interacton for spheres

inline CInteractonSphere::CInteractonSphere (Particle * Pt1, Particle * Pt2)
{
    P1   = Pt1;
    P2   = Pt2;
    I1 = P1->Index;
    I2 = P2->Index;
    Kn   = 2*ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
    Kt   = 2*ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
    Gn   = 2*ReducedValue(Pt1->Props.Gn,Pt2->Props.Gn)*ReducedValue(Pt1->Props.m,Pt2->Props.m);
    Gt   = 2*ReducedValue(Pt1->Props.Gt,Pt2->Props.Gt)*ReducedValue(Pt1->Props.m,Pt2->Props.m);
    Mu   = 2*ReducedValue(Pt1->Props.Mu,Pt2->Props.Mu);
    beta = 2*ReducedValue(Pt1->Props.Beta,Pt2->Props.Beta);
    eta  = 2*ReducedValue(Pt1->Props.Eta,Pt2->Props.Eta);
    Nc = 0;
    Nsc = 0;

    Epot = 0.0;
    Fdr  = 0.0, 0.0, 0.0;
    Lvv.Push(make_pair(0,0));

    CalcForce(0.0);
}

inline void CInteractonSphere::_update_rolling_resistance(double dt)
{
    Vec3_t t1,t2;
    Rotation(P1->w,P1->Q,t1);
    Rotation(P2->w,P2->Q,t2);
    Vec3_t Normal = Fn/norm(Fn);
    Vec3_t Vr = P1->Props.R*P2->Props.R*cross(Vec3_t(t1 - t2),Normal)/(P1->Props.R+P2->Props.R);
    Fdr += Vr*dt;
    Fdr -= dot(Fdr,Normal)*Normal;
    Vec3_t tan = Fdr;
    if (norm(tan)>0.0) tan/=norm(tan);
    double Kr = beta*Kt;
    if (norm(Fdr)>eta*Mu*norm(Fn)/Kr)
    {
        Fdr = eta*Mu*norm(Fn)/Kr*tan;
        Nr++;
    }
    Vec3_t Ft = -Kr*Fdr;

    Vec3_t Tt = P1->Props.R*cross(Normal,Ft);
    Vec3_t T;
#ifdef USE_THREAD
    std::lock_guard<std::mutex> lk1(P1->mtex);
    std::lock_guard<std::mutex> lk2(P2->mtex);
#endif
    Quaternion_t q;
    Conjugate (P1->Q,q);
    Rotation  (Tt,q,T);
    P1->T += T;

    Tt = P2->Props.R*cross(Normal,Ft);
    Conjugate (P2->Q,q);
    Rotation  (Tt,q,T);
    P2->T -= T;
}

inline void CInteractonSphere::CalcForce(double dt)
{
    Epot   = 0.0;
    dEvis  = 0.0;
    dEfric = 0.0;
    Nc     = 0.0;
    Nsc    = 0.0;
    Nr     = 0.0;
    Fnet   = 0.0;
    Ftnet  = 0.0;
    Xc     = 0.0;
    _update_disp_calc_force (P1->Verts,P2->Verts,Fdvv,Lvv,dt);
    if (Epot>0.0) _update_rolling_resistance(dt);


    //If there is at least a contact, increase the coordination number of the particles
    if (Nc>0) 
    {
        P1->Cn++;
        P2->Cn++;
    }
}

inline bool CInteractonSphere::UpdateContacts (double alpha)
{
    if (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*alpha) return true;
    else return false;
}

inline void CInteractonSphere::UpdateParameters ()
{
    Kn   = 2*ReducedValue(P1->Props.Kn,P2->Props.Kn);
    Kt   = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
    Gn   = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn)*ReducedValue(P1->Props.m,P2->Props.m);
    Gt   = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt)*ReducedValue(P1->Props.m,P2->Props.m);
    Mu   = 2*ReducedValue(P1->Props.Mu,P2->Props.Mu);
    beta = 2*ReducedValue(P1->Props.Beta,P2->Props.Beta);
    eta  = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);
}

//Cohesion interacton

inline BInteracton::BInteracton (Particle * Pt1, Particle * Pt2, size_t Fi1, size_t Fi2)
{
    P1              = Pt1;
    P2              = Pt2;
    F1              = Fi1;
    if (Fi2>=P2->Faces.Size())
    {
        F2 = P2->Faces.Size()-1;
        Area            = P1->Faces[F1]->Area();
    }
    else
    {
        F2 = Fi2;
        Area            = 0.5*(P1->Faces[F1]->Area()+P2->Faces[F2]->Area());
    }
    I1              = P1->Index;
    I2              = P2->Index;
    Bn              = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn)*Area;
    Bt              = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt)*Area;
    Bm              = 2*ReducedValue(P1->Props.Bm,P2->Props.Bm)*Area;
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn)*ReducedValue(P1->Props.m,P2->Props.m);
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt)*ReducedValue(P1->Props.m,P2->Props.m);
    eps             = 2*ReducedValue(P1->Props.eps,P2->Props.eps);

    Vec3_t n1,n2;
    P1->Faces[F1]->Normal(n1);
    P2->Faces[F2]->Normal(n2);
    Vec3_t c1,c2;
    An         = 0.0;
    valid      = true;
    P1->Faces[F1]->Centroid(c1);
    P2->Faces[F2]->Centroid(c2);
    L0         = dot(n1,c2-c1);
    if (Fi2>=P2->Faces.Size()) c2 = c1;
    Vec3_t V   = 0.5*(c1+c2);

    Face * F   = P1->Faces[F1];
    double a   = dot(*F->Edges[0]->X0-V , F->Edges[0]->dL);
    double b   = dot(F->Edges[0]->dL    , F->Edges[0]->dL);
    double c   = dot(F->Edges[0]->dL    , F->Edges[1]->dL);
    double d   = dot(*F->Edges[0]->X0-V , F->Edges[1]->dL);
    double f   = dot(F->Edges[1]->dL    , F->Edges[1]->dL);
    s1         = (c*d-a*f)/(b*f-c*c);
    t1         = (a*c-b*d)/(b*f-c*c);
    
    Vec3_t p1  = -s1*F->Edges[0]->dL - t1*F->Edges[1]->dL;

    F          = P2->Faces[F2];
    a          = dot(*F->Edges[0]->X0-V , F->Edges[0]->dL);
    b          = dot(F->Edges[0]->dL    , F->Edges[0]->dL);
    c          = dot(F->Edges[0]->dL    , F->Edges[1]->dL);
    d          = dot(*F->Edges[0]->X0-V , F->Edges[1]->dL);
    f          = dot(F->Edges[1]->dL    , F->Edges[1]->dL);
    s2         = (c*d-a*f)/(b*f-c*c);
    t2         = (a*c-b*d)/(b*f-c*c);

    Vec3_t p2  = -s2*F->Edges[0]->dL - t2*F->Edges[1]->dL;
    
    double arg = dot(p1,p2)/(norm(p1)*norm(p2));
    if (arg> 1.0) arg= 1.0;
    if (arg<-1.0) arg=-1.0;
    An         = acos(arg);
    if (dot(cross(p1,p2),n1)<0) An = -An;
}

inline bool BInteracton::UpdateContacts (double alpha)
{
    return valid;
}

inline void BInteracton::CalcForce(double dt)
{
    if (valid)
    {
        // Calculate the normal vector and centroid of the contact face
        Vec3_t n1,n2,n;
        P1->Faces[F1]->Normal(n1);
        P2->Faces[F2]->Normal(n2);
        n = 0.5*(n1-n2);
        n/=norm(n);

        Face * F    = P1->Faces[F1];
        Vec3_t pro1 = *F->Edges[0]->X0 + s1*F->Edges[0]->dL + t1*F->Edges[1]->dL;
        Vec3_t p1   = -s1*F->Edges[0]->dL - t1*F->Edges[1]->dL;
        F           = P2->Faces[F2];
        Vec3_t pro2 = *F->Edges[0]->X0 + s2*F->Edges[0]->dL + t2*F->Edges[1]->dL;
        Vec3_t p2   = -s2*F->Edges[0]->dL - t2*F->Edges[1]->dL;

        // Normal force
        double delta = (dot(pro2-pro1,n)-L0)/L0;
        if (delta<0.0) delta = 0.0;
        Vec3_t Fn = -Bn*(delta)*n;

        // Tangential Force
        Vec3_t x = 0.5*(pro1+pro2);
        Vec3_t x1,x2;
        x1 = x - P1->x;
        x2 = x - P2->x;
        Vec3_t td   = pro2-pro1-dot(pro2-pro1,n)*n;
        Vec3_t Ft   = -Bt*td/L0;

#ifdef USE_THREAD
        std::lock_guard<std::mutex> lk1(P1->mtex);
        std::lock_guard<std::mutex> lk2(P2->mtex);
#endif
        //Adding forces and torques
        Fnet   = Fn+Ft;
        P1->F -= Fnet;
        P2->F += Fnet;
        Vec3_t T, Tt;
        Tt = cross (x1,Fnet);
        Quaternion_t q1;
        Conjugate (P1->Q,q1);
        Rotation  (Tt,q1,T);
        P1->T -= T;
        Tt = cross (x2,Fnet);
        Quaternion_t q2;
        Conjugate (P2->Q,q2);
        Rotation  (Tt,q2,T);
        P2->T += T;

        //Torque
        double arg = dot(p1,p2)/(norm(p1)*norm(p2));
        if (arg> 1.0) arg= 1.0;
        if (arg<-1.0) arg=-1.0;
        double Ant = acos(arg);
        if (dot(cross(p1,p2),n)<0) Ant = -Ant;
        Tt         = -Bm*(Ant-An)*n/L0;
        Rotation  (Tt,q1,T);
        P1->T -= T;
        Rotation  (Tt,q2,T);
        P2->T += T;

        //Breaking point
        if ((L0*fabs(delta)+norm(td)+0.0*fabs(Ant-An)*L0>L0*eps)&&(eps>=0.0))
        {
            valid = false;
        }
    }
}

inline void BInteracton::UpdateParameters ()
{
    Bn              = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn)*Area;
    Bt              = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt)*Area;
    Bm              = 2*ReducedValue(P1->Props.Bm,P2->Props.Bm)*Area;
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn)*ReducedValue(P1->Props.m,P2->Props.m);
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt)*ReducedValue(P1->Props.m,P2->Props.m);
    eps             = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
}
#endif //  MECHSYS_DEM_INTERACTON_H
