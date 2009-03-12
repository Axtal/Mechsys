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

#ifndef __GRAPH_H
#define __GRAPH_H

#include <stdlib.h>
#include <fstream>
#include <string>
#include "Interacton.h"

class CGraph {
    protected:
        ofstream m_file;
        string m_filename;
    public:
        CGraph (char *a);
        void Graph_Face3d(CFace3d &,double,char *);
        void Graph_Edge3d(CEdge3d &,double,char *);
        void Graph_Vec3(Vec3 &,double,char *);
        void Graph_Polygon(Vec3 *,int,char *);
        void Graph_Particle(CParticle &,char *);
        void Graph_Contact(CInteracton &,char *);
        void Animation(const int & t);
        void close(void) {m_file.close();};
        
};
#endif
