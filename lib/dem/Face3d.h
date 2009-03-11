#ifndef __FACE3D_H
#define __FACE3D_H

#include <vector>
#include <map>
#include "Edge3d.h"

class CFace3d {
    protected:
        int m_Ns;
        vector<CEdge3d> m_sides;
        //CEdge3d m_sides[4];
    public:
        CFace3d() {};
        CFace3d(Vec3 *,int N);
        CFace3d(int N) {m_Ns=N;};
        void translation(Vec3 & a);
        void rotation(Quaternion &,Vec3 &);
        friend CFace3d rotation_face(const CFace3d &,Quaternion & Q,Vec3 & v);
        friend CEdge3d distance_point_Face3d(Vec3 &,CFace3d &);
        friend class CGraph;
        friend class CParticle;
        friend class CEdge3d;
};


#endif
