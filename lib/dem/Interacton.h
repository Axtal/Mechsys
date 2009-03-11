#ifndef __INTERACTON_H
#define __INTERACTON_H

#include "Particle.h"

class CInteracton {
    protected:
        CParticle *m_P1,*m_P2;
        double m_kn,m_kt,m_Epot,m_Edis,m_gn,m_gt,m_mu;
        map<pair<int,int>,Vec3> m_free,m_frvf,m_frfv;
        Vec3 m_fn;
        bool m_c;
    public:
        CInteracton (void) {};
        CInteracton(CParticle &p1,CParticle &p2);
        void calcForce(double);
        double getEpot(void) {return m_Epot;};
        Vec3 getFn(void) {return m_fn;};
        bool getContact(void) {return m_c;};
        double getmu(void) {return m_mu;};
        void setmu(double mu) {m_mu=mu;};
	    void setkn(double k) {m_kn=k;};
	    void setkt(double k) {m_kt=k;};
        friend class CGraph;
};



#endif
