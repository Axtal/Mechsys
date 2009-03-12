#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "Face3d.h"

class CParticle {
    protected:
        int m_Nf,m_Ne,m_Nv;
        //CFace3d m_Face[8],m_Faceold[8];
        //CEdge3d m_Edge[12],m_Edgeold[12];
        //Vec3 m_ve[8],m_veold[8];
        vector<CFace3d> m_Face,m_Faceold;
        vector<CEdge3d> m_Edge,m_Edgeold;
        vector<Vec3> m_ve,m_veold;
        double m_m,m_R,m_Ixx,m_Iyy,m_Izz,m_Ekin,m_Erot,m_dmax;
        Vec3 m_r,m_v,m_f,m_T,m_Td,m_ome,m_fome,m_r1,m_r2;
        Quaternion m_Q;
        bool m_large;
    public:
        CParticle (void) {};
        CParticle(double l0, double m0, double R0,Vec3 &r0, Vec3 &v0,  Vec3 &O,Quaternion &Q0,string p);//General particle constructor
        CParticle(double a,double b,double R0,Vec3 &r0,Quaternion &Q0,const Vec3 & v0=Vec3(0,0,0));//plane cosntructor
        CParticle(double a,double b,double c,double R0,Vec3 &r0,Quaternion &Q0);//hopper constructor
        CParticle(double a,double b,int N,double R0,double l0,Vec3 &r0,Vec3 &v0,Quaternion &Q0);//rough surface cosntructor;
        void rotation_quaternion(Quaternion &);
        void rotation(double dt);
        void start(double);
        void translation(double);
        void translation(Vec3 &);
        void startforce(const Vec3 & F=Vec3(0,0,0)) {m_f=F; m_T=Vec3(0,0,0);};
        void setT(const Vec3 & F=Vec3(0,0,0)) {m_T=F;};
        void setm(double m) {m_m=m;};
        void setv(const Vec3 & F=Vec3(0,0,0)) {m_v=F;};
        double getXmax(void) {
            double x=getVer(0).X();
            for (int i=1;i<m_Nv;i++) {
                if (getVer(i).X()>x) x=getVer(i).X();
            }
            return x+m_R;
        }
        double getXmin(void) {
            double x=getVer(0).X();
            for (int i=1;i<m_Nv;i++) {
                if (getVer(i).X()<x) x=getVer(i).X();
            }
            return x-m_R;
        }
        double getEkin(void) {return m_Ekin;};
        double getErot(void) {return m_Erot;};
        double getm(void) {return m_m;};
        double getIxx(void) {return m_Ixx;};
        double getIyy(void) {return m_Iyy;};
        double getIzz(void) {return m_Izz;};
        double getdmax(void) {return m_dmax;};
        double getPower(void) {return dot(m_f,m_v)+dot(m_T,m_ome);};
        Vec3 getVer(int i) {return m_ve[i];};
        Vec3 getr(void) {return m_r;};
        Vec3 getome(void) {return m_ome;};
        Vec3 getfome(void) {return m_fome;};
        Vec3 getv(void) {return m_v;};
        Vec3 getT(void) {return m_T;};
        Vec3 getTd(void) {return m_Td;};
        Vec3 getF(void) {return m_f;};
        Vec3 getp(void) {return m_m*m_v;};
        Vec3 getL(void) {return Vec3(m_Ixx*m_ome.X(),m_Iyy*m_ome.Y(),m_Izz*m_ome.Z());};
        Vec3 getLs(void) {Vec3 v=getL(); v=m_m*cross(m_r,m_v)+rotate(m_Q,v); return v;};
        Vec3 getvertex(int i) {return m_ve[i];};
        Quaternion getQ(void) {return m_Q;};
        friend class CGraph;
        friend class CInteracton;
            

};


#endif
