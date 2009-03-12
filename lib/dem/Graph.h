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
