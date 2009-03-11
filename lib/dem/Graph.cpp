#include "Graph.h"


CGraph::CGraph(char *a) {
    m_filename="";
    m_filename.append(a);
    m_file.open(a);
    m_file << "#include \"colors.inc\" \n";
    m_file << "#include \"glass.inc\" \n";
    m_file << "background {color White} \n";
    m_file << "light_source{<3,0,0> color White shadowless}  \n";
    m_file << "light_source{<-3,0,0> color White shadowless}  \n";
    m_file << "light_source{<0,3,0> color White shadowless}  \n";
    m_file << "light_source{<0,-3,0> color White shadowless}  \n";
    //file << "global_settings {ambient_light White} \n";
    m_file << "camera { location <0,17,3> sky <0,0,1> look_at <0,0,0> }\n";
}


void CGraph::Graph_Vec3(Vec3 & E,double R,char *s) {
    m_file << "sphere  { <"<<E.X()<<","<<E.Y()<<","<<E.Z()<<">,"<<R<<"\n pigment { color "<<s<<" } }\n";   
}

void CGraph::Graph_Edge3d(CEdge3d & E,double R,char *s) {
    Vec3 IN=E.I(),FI=E.F();
    m_file << "cylinder { <"<<IN.X()<<","<<IN.Y()<<","<<IN.Z()<<">,<"<<FI.X()<<","<<FI.Y()<<","<<FI.Z()<<">,"<<R<<"\n pigment { color "<<s<<" } }\n";
}

void CGraph::Graph_Face3d(CFace3d & E,double R,char *s) {
    int i;
    Vec3 vi[10],vs[10],n;
    n=cross(E.m_sides[0].m_dx,E.m_sides[1].m_dx);
    n=n.unit();
    //n.show_on_screen(); cout << endl;
    for(i=0;i<=E.m_Ns;i++) {
        //Graph_Edge3d(E.m_sides[i],R,s);
        //Graph_Vec3(E.m_sides[i].m_xi,R,s);
        vi[i]=E.m_sides[i%E.m_Ns].I()-R*n;
        vs[i]=E.m_sides[i%E.m_Ns].I()+R*n;
    }
    Graph_Polygon(vi,E.m_Ns,s);
    Graph_Polygon(vs,E.m_Ns,s);
    
}

void CGraph::Graph_Polygon(Vec3 *v,int N,char *s) {
    int i;
    m_file << "polygon {"<<N<<", \n";
    i=0;
    m_file << "<"<<v[i].X()<<","<<v[i].Y()<<","<<v[i].Z()<<">";
    for(i=1;i<N;i++) {
        m_file << ",<"<<v[i].X()<<","<<v[i].Y()<<","<<v[i].Z()<<">";
    }
    m_file <<"\n pigment { color "<<s<<" } }\n";
}

void CGraph::Graph_Particle(CParticle & P,char *s) {
    int i;
    //cout<<P.m_R<< endl;
    for (i=0;i<P.m_Nf;i++) {
        if(P.m_large) Graph_Face3d(P.m_Face[i],P.m_R,"Col_Glass_Bluish");
        else Graph_Face3d(P.m_Face[i],P.m_R,s);
    }
    for (i=0;i<P.m_Ne;i++) {
        Graph_Edge3d(P.m_Edge[i],P.m_R,s);
    }
    for (i=0;i<P.m_Nv;i++) {
        Graph_Vec3(P.m_ve[i],P.m_R,s);
    }
}

void CGraph::Graph_Contact(CInteracton & I,char *s) {
    if (I.m_c) {
        Vec3 xi,xf;
        xi=I.m_P1->m_r;
        xf=I.m_P2->m_r;
        CEdge3d E(xi,xf);
        double R=0.001*I.m_fn.norm();
        Graph_Edge3d(E,R,s);
    }
}


void CGraph::Animation(const int & t) {
       string p="",s="";
       char b[10];
       p=m_filename;
       sprintf(b,"%d",t);
       p.append(b);
       p.append(".png");
       const char *a=p.c_str();
       s="povray +W480 +H480 +FN +I";
       s=s+m_filename;
       s=s+" +O";
       s.append(a);
       const char *c=s.c_str();
       system(c);
}





