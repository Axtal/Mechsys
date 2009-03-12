/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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
