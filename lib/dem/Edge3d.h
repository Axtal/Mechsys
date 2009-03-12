#ifndef __EDGE3D_H
#define __EDGE3D_H


#include "Quaternion.h"
#include <iostream>

using namespace std;

/** Class defining a single line segment (edge) in 3D */
class CEdge3d {
    protected:
        Vec3 m_xi,m_xf,m_dx;
        double m_lenght;
    public:
        CEdge3d(void);
        CEdge3d(const Vec3 &,const Vec3 &);
        Vec3 I(void) {return m_xi;};
        Vec3 F(void) {return m_xf;};
        Vec3 rc(void) {return (m_xi+m_xf)/2;};
        Vec3 DX(void) {return m_dx;};
        double L(void) {return m_lenght;};
        void translate (const Vec3 & a) {m_xi+=a; m_xf+=a;};
        void rotation (Quaternion & Q,Vec3 & v);
        void invert(void) {m_dx=m_xi; m_xi=m_xf; m_xf=m_dx; m_dx=m_xf-m_xi;};
        friend CEdge3d distance_point_Edge3d(const Vec3 &,const CEdge3d &);
        friend CEdge3d distance_Edge3d_Edge3d(const CEdge3d &,const CEdge3d &);
        friend CEdge3d rotation_edge(CEdge3d &,Quaternion &,Vec3 &);
        friend Vec3 cross(const Vec3&,const Vec3&);
        friend double dot(const Vec3&,const Vec3&);
        friend class CFace3d;
        friend class CGraph;
        friend class Vec3;
};




#endif
