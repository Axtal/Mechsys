#include "Particle.h"
#include <stdlib.h>

CParticle::CParticle(double l0, double m0, double R0,Vec3 &r0, Vec3 &v0,  Vec3 &O,Quaternion &Q0,string p) {
    double wx,wy,wz;
    int i,j,n;
    //cout << m_Face<< " " <<m_Faceold<< endl;
    m_large=false;
    //if(p=="tetra") {
    vector<CFace3d>::iterator itf;
    vector<CEdge3d>::iterator ite;
    vector<Vec3>::iterator itv;
    Vec3 tempv;
    if(p=="yermis") {
        m_dmax=l0/2+R0;
        m_Nf=0;
        m_Ne=3;
        m_Nv=6;
        m_Q=Q0.conj();
        m_r=r0;
        //m_r.show_on_screen(); cout <<endl;
        m_v=v0;
        m_ome=m_fome=rotate(m_Q,O);
        m_Q=Q0;
        m_T=Vec3(0,0,0);
        m_f=Vec3(0,0,0);
        m_R=R0;
        m_m=3*((4./3.)*M_PI*m0*R0*R0*R0+m0*M_PI*R0*R0*l0);
        
        tempv=Vec3(0,0,-l0/2);
        itv=m_veold.end();
        m_veold.insert(itv,tempv);
        
        tempv=Vec3(0,0,l0/2);
        itv=m_veold.end();
        m_veold.insert(itv,tempv);
        
        tempv=Vec3(0,-l0/2,0);
        itv=m_veold.end();
        m_veold.insert(itv,tempv);
        
        tempv=Vec3(0,l0/2,0);
        itv=m_veold.end();
        m_veold.insert(itv,tempv);
        
        tempv=Vec3(-l0/2,0,0);
        itv=m_veold.end();
        m_veold.insert(itv,tempv);
        
        tempv=Vec3(l0/2,0,0);
        itv=m_veold.end();
        m_veold.insert(itv,tempv);
        
        for(i=0;i<m_Nv;i++) {
            itv=m_ve.end();
            m_ve.insert(itv,rotate(m_Q,m_veold[i])+m_r);
        }
        
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[0],m_veold[1]));
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[2],m_veold[3]));
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[4],m_veold[5]));
        
        for(i=0;i<m_Ne;i++) {
            m_Edgeold[i].translate(m_r);
            ite=m_Edge.end();
            m_Edge.insert(ite,rotation_edge(m_Edgeold[i],m_Q,m_r));
        }
        m_Izz=m_Ixx=m_Iyy=(2./5.)*M_PI*m_R*m_R*l0*m0*m_R*m_R+0.5*(4*M_PI/3)*m_R*m_R*m_R*m0*m_R*m_R+(1./12.)*M_PI*m_R*m_R*l0*m0*(3*m_R*m_R+l0*l0)+(2./5.)*(4*M_PI/3)*m_R*m_R*m_R*m0*m_R*m_R+(4*M_PI/3)*m_R*m_R*m_R*m0*0.25*l0*l0;
       
        wx=m_ome.X();
        wy=m_ome.Y();
        wz=m_ome.Z();

        m_Erot=0.5*(m_Ixx*wx*wx+m_Iyy*wy*wy+m_Izz*wz*wz);
        m_Ekin=0.5*m_m*dot(m_v,m_v);
    }
    
        
    if(p=="rice") {
        m_dmax=l0/2+R0;
        m_Nf=0;
        m_Ne=1;
        m_Nv=2;
        m_Q=Q0.conj();
        m_r=r0;
        //m_r.show_on_screen(); cout <<endl;
        m_v=v0;
        m_ome=m_fome=rotate(m_Q,O);
        m_Q=Q0;
        m_T=Vec3(0,0,0);
        m_f=Vec3(0,0,0);
        m_R=R0;
        m_m=(4./3.)*M_PI*m0*R0*R0*R0+m0*M_PI*R0*R0*l0;

        tempv=Vec3(0,0,-l0/2);
        itv=m_veold.end();
        m_veold.insert(itv,tempv);

        tempv=Vec3(0,0,l0/2);
        itv=m_veold.end();
        m_veold.insert(itv,tempv);

        for(i=0;i<m_Nv;i++) {
            itv=m_ve.end();
            m_ve.insert(itv,rotate(m_Q,m_veold[i])+m_r);
        }

        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[0],m_veold[1]));

        for(i=0;i<m_Ne;i++) {
            m_Edgeold[i].translate(m_r);
            ite=m_Edge.end();
            m_Edge.insert(ite,rotation_edge(m_Edgeold[i],m_Q,m_r));
        }
        m_Izz=(2./5.+0.5)*m_m*m_R*m_R;
        m_Ixx=m_Iyy=(2./5.)*m_m*m_R*m_R+0.25*m_m*l0*l0+(1./12.)*m_m*(3*m_R*m_R+l0*l0);

        wx=m_ome.X();
        wy=m_ome.Y();
        wz=m_ome.Z();

        m_Erot=0.5*(m_Ixx*wx*wx+m_Iyy*wy*wy+m_Izz*wz*wz);
        m_Ekin=0.5*m_m*dot(m_v,m_v);
    }
    
    if(p=="tetra") {
        m_dmax=sqrt(3./8.)*l0+R0;
        m_Nf=4;
        m_Ne=6;
        m_Nv=4;
        Vec3 v[3];
        //m_Face=new CFace3d[m_Nf];
        //m_Faceold=new CFace3d[m_Nf];
        //cout << m_Face<< " " <<m_Faceold<< endl;
        m_Q=Q0.conj();
        m_r=r0;
        //m_r.show_on_screen(); cout <<endl;
        m_v=v0;
        m_ome=m_fome=rotate(m_Q,O);
        m_Q=Q0;
        m_T=Vec3(0,0,0);
        m_f=Vec3(0,0,0);
        m_R=R0;
        m_m=l0*l0*l0*sqrt(2)/12.;
        
        v[0]=Vec3(l0/2,0,l0/sqrt(8));
        v[1]=Vec3(-l0/2,0,l0/sqrt(8));
        v[2]=Vec3(0,-l0/2,-l0/sqrt(8));
        itf=m_Faceold.end();
        m_Faceold.insert(itf,CFace3d(v,3));
       
        v[0]=Vec3(-l0/2,0,l0/sqrt(8));
        v[1]=Vec3(l0/2,0,l0/sqrt(8));
        v[2]=Vec3(0,l0/2,-l0/sqrt(8));
        itf=m_Faceold.end();
        m_Faceold.insert(itf,CFace3d(v,3));

        v[0]=Vec3(-l0/2,0,l0/sqrt(8));
        v[1]=Vec3(0,l0/2,-l0/sqrt(8));
        v[2]=Vec3(0,-l0/2,-l0/sqrt(8));
        itf=m_Faceold.end();
        m_Faceold.insert(itf,CFace3d(v,3));

        v[0]=Vec3(l0/2,0,l0/sqrt(8));
        v[1]=Vec3(0,-l0/2,-l0/sqrt(8));
        v[2]=Vec3(0,l0/2,-l0/sqrt(8));
        itf=m_Faceold.end();
        m_Faceold.insert(itf,CFace3d(v,3));
        
        
        for(i=0;i<m_Nf;i++) {
            m_Faceold[i].translation(m_r);
            itf=m_Face.end();
            m_Face.insert(itf,rotation_face(m_Faceold[i],m_Q,m_r));
            //cout << "1"<<endl;
        }
        
        itv=m_veold.end();
        m_veold.insert(itv,Vec3(l0/2,0,l0/sqrt(8)));
        
        itv=m_veold.end();
        m_veold.insert(itv,Vec3(-l0/2,0,l0/sqrt(8)));
        
        itv=m_veold.end();
        m_veold.insert(itv,Vec3(0,-l0/2,-l0/sqrt(8)));
        
        itv=m_veold.end();
        m_veold.insert(itv,Vec3(0,l0/2,-l0/sqrt(8)));
        
        for(i=0;i<m_Nv;i++) {
            itv=m_ve.end();
            m_ve.insert(itv,rotate(m_Q,m_veold[i])+m_r);
        }

        //cout<<"1"<<endl;
        for(i=0;i<m_Nv-1;i++) {
            for(j=i+1;j<m_Nv;j++) {
                ite=m_Edgeold.end();
                m_Edgeold.insert(ite,CEdge3d(m_veold[i],m_veold[j]));
                //cout << i << " " << j <<endl;
            }
        }
        for(i=0;i<m_Ne;i++) {
            m_Edgeold[i].translate(m_r);
            ite=m_Edge.end();
            m_Edge.insert(ite,rotation_edge(m_Edgeold[i],m_Q,m_r));
        }


        m_Ixx=m_Iyy=m_Izz=m_m*l0*l0/20;

        wx=m_ome.X();
        wy=m_ome.Y();
        wz=m_ome.Z();

        m_Erot=0.5*(m_Ixx*wx*wx+m_Iyy*wy*wy+m_Izz*wz*wz);
        m_Ekin=0.5*m_m*dot(m_v,m_v);


    }
    
    if (p=="tippetop") {
        m_dmax=2*R0;
        m_dmax=sqrt(3./8.)*l0+R0;
        m_Nf=0;
        m_Ne=0;
        m_Nv=1;
        m_Q=Q0.conj();
        m_r=r0;
        m_v=v0;
        m_ome=m_fome=rotate(m_Q,O);
        m_Q=Q0;
        m_T=Vec3(0,0,0);
        m_f=Vec3(0,0,0);
        m_R=R0;
        m_m=m0;
        
        itv=m_veold.end();
        m_veold.insert(itv,Vec3(0,0,l0));
        
        itv=m_ve.end();
        m_ve.insert(itv,rotate(m_Q,m_veold[0])+m_r);
        
    
        m_Ixx=m_Iyy=m_Izz=(2./5.)*m_m*m_R*m_R;

        wx=m_ome.X();
        wy=m_ome.Y();
        wz=m_ome.Z();

        m_Erot=0.5*(m_Ixx*wx*wx+m_Iyy*wy*wy+m_Izz*wz*wz);
        m_Ekin=0.5*m_m*dot(m_v,m_v);
        
        
    }
    
    if (p=="pirinola") {
        m_dmax=2*l0;
        Vec3 v[6];
        m_Nf=13;
        m_Ne=24;
        m_Nv=13;
        m_Q=Q0.conj();
        m_r=r0;
        m_v=v0;
        m_ome=m_fome=O;
        m_Q=Q0;
        m_T=Vec3(0,0,0);
        m_f=Vec3(0,0,0);
        m_R=R0;
        m_m=m0;

        itv=m_veold.end();
        m_veold.insert(itv,Vec3(0,0,-l0));
        for (i=0;i<6;i++) {
            itv=m_veold.end();
            m_veold.insert(itv,Vec3(l0*cos(M_PI*i/3)/2,l0*sin(M_PI*i/3)/2,-l0/2));
            itv=m_veold.end();
            m_veold.insert(itv,Vec3(l0*cos(M_PI*i/3)/2,l0*sin(M_PI*i/3)/2,l0/2));
            v[0]=Vec3(l0*cos(M_PI*i/3)/2,l0*sin(M_PI*i/3)/2,-l0/2);
            v[1]=Vec3(l0*cos(M_PI*i/3)/2,l0*sin(M_PI*i/3)/2,l0/2);
            v[2]=Vec3(l0*cos(M_PI*(i+1)/3)/2,l0*sin(M_PI*(i+1)/3)/2,l0/2);
            v[3]=Vec3(l0*cos(M_PI*(i+1)/3)/2,l0*sin(M_PI*(i+1)/3)/2,-l0/2);
            ite=m_Edgeold.end();
            m_Edgeold.insert(ite,CEdge3d(v[0],v[1]));
            ite=m_Edgeold.end();
            m_Edgeold.insert(ite,CEdge3d(v[0],v[3]));
            ite=m_Edgeold.end();
            m_Edgeold.insert(ite,CEdge3d(v[1],v[2]));
            itf=m_Faceold.end();
            m_Faceold.insert(itf,CFace3d(v,4));
            v[0]=Vec3(l0*cos(M_PI*i/3)/2,l0*sin(M_PI*i/3)/2,-l0/2);
            v[1]=Vec3(l0*cos(M_PI*(i+1)/3)/2,l0*sin(M_PI*(i+1)/3)/2,-l0/2);
            v[2]=Vec3(0,0,-l0);
            ite=m_Edgeold.end();
            m_Edgeold.insert(ite,CEdge3d(v[0],v[2]));
            itf=m_Faceold.end();
            m_Faceold.insert(itf,CFace3d(v,3));
        }
        for (i=0;i<6;i++) {
            v[i]=Vec3(l0*cos(M_PI*i/3)/2,l0*sin(M_PI*i/3)/2,l0/2);
        }
        itf=m_Faceold.end();
        m_Faceold.insert(itf,CFace3d(v,6));
        
        for(i=0;i<m_Ne;i++) {
            m_Edgeold[i].translate(m_r);
            ite=m_Edge.end();
            m_Edge.insert(ite,rotation_edge(m_Edgeold[i],m_Q,m_r));
        }
        
        for(i=0;i<m_Nf;i++) {
            m_Faceold[i].translation(m_r);
            itf=m_Face.end();
            m_Face.insert(itf,rotation_face(m_Faceold[i],m_Q,m_r));
        }
        
        for(i=0;i<m_Nv;i++) {
            itv=m_ve.end();
            m_ve.insert(itv,rotate(m_Q,m_veold[i])+m_r);
        }
        
        
        

        m_Ixx=m_Iyy=(1./12.)*m_m*(3*(l0/2)*(l0/2)+l0*l0);
        m_Izz=0.125*m_m*l0*l0;
        wx=m_ome.X();
        wy=m_ome.Y();
        wz=m_ome.Z();

        m_Erot=0.5*(m_Ixx*wx*wx+m_Iyy*wy*wy+m_Izz*wz*wz);
        m_Ekin=0.5*m_m*dot(m_v,m_v);


    }


}

CParticle::CParticle(double a,double b,double R0,Vec3 &r0,Quaternion &Q0, const Vec3 & v0) {
    double wx,wy,wz;
    vector<CFace3d>::iterator itf;
    vector<CEdge3d>::iterator ite;
    vector<Vec3>::iterator itv;
    m_Nf=1;
    m_Ne=4;
    m_Nv=4;
    Vec3 v[4];
    int i;
    m_r=r0;
    m_R=R0;
    m_Q=Q0;
    m_large=true;
    m_m=1;
    m_dmax=sqrt((a/2)*(a/2)+(b/2)*(b/2))+R0;
    m_v=v0;
    m_T=Vec3(0,0,0);
    m_f=Vec3(0,0,0);
    m_ome=m_fome=Vec3(0,0,0);    
    v[0]=Vec3(a/2,b/2,0);
    v[1]=Vec3(-a/2,b/2,0);
    v[2]=Vec3(-a/2,-b/2,0);
    v[3]=Vec3(a/2,-b/2,0);
    itf=m_Faceold.end();
    m_Faceold.insert(itf,CFace3d(v,4));
    
    i=0;
    m_Faceold[i].translation(m_r);
    itf=m_Face.end();
    m_Face.insert(itf,rotation_face(m_Faceold[i],m_Q,m_r));
    
    for(i=0;i<m_Nv;i++) {
        itv=m_veold.end();
        m_veold.insert(itv,v[i]);
    }
    for(i=0;i<m_Nv;i++) {
        itv=m_ve.end();
        m_ve.insert(itv,rotate(m_Q,m_veold[i])+m_r);
    }
   
    
    for(i=0;i<m_Nv;i++) {
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[i],m_veold[(i+1)%m_Nv]));
        m_Edgeold[i].translate(m_r);
        ite=m_Edge.end();
        m_Edge.insert(ite,rotation_edge(m_Edgeold[i],m_Q,m_r));
    }
    
    wx=m_ome.X();
    wy=m_ome.Y();
    wz=m_ome.Z();
    
    m_Erot=0.5*(m_Ixx*wx*wx+m_Iyy*wy*wy+m_Izz*wz*wz);
    m_Ekin=0.5*m_m*dot(m_v,m_v);
}

CParticle::CParticle(double a,double b,double c,double R0,Vec3 &r0,Quaternion &Q0) {
    double wx,wy,wz;
    vector<CFace3d>::iterator itf;
    vector<CEdge3d>::iterator ite;
    vector<Vec3>::iterator itv;
    Vec3 v[4];
    int i,j;
    m_Nf=4;
    m_Ne=12;
    m_Nv=8;
    m_r=r0;
    m_R=R0;
    //cout << m_R<< endl;
    m_Q=Q0;
    m_large=true;
    m_dmax=sqrt((a/2)*(a/2)+(c/2)*(c/2))+R0;
    m_v=Vec3(0,0,0);
    m_ome=m_fome=Vec3(0,0,0);
    
    itv=m_veold.end();
    m_veold.insert(itv,Vec3(a/2,a/2,c/2));
    
    itv=m_veold.end();
    m_veold.insert(itv,Vec3(a/2,-a/2,c/2));
    
    itv=m_veold.end();
    m_veold.insert(itv,Vec3(-a/2,-a/2,c/2));
    
    itv=m_veold.end();
    m_veold.insert(itv,Vec3(-a/2,a/2,c/2));
    
    itv=m_veold.end();
    m_veold.insert(itv,Vec3(b/2,b/2,-c/2));
    
    itv=m_veold.end();
    m_veold.insert(itv,Vec3(b/2,-b/2,-c/2));
    
    itv=m_veold.end();
    m_veold.insert(itv,Vec3(-b/2,-b/2,-c/2));
    
    itv=m_veold.end();
    m_veold.insert(itv,Vec3(-b/2,b/2,-c/2));
    
    for(i=0;i<m_Nv;i++) {
        itv=m_ve.end();
        m_ve.insert(itv,rotate(m_Q,m_veold[i])+m_r);
    }
    
    v[0]=Vec3(a/2,a/2,c/2);
    v[1]=Vec3(b/2,b/2,-c/2);
    v[2]=Vec3(b/2,-b/2,-c/2);
    v[3]=Vec3(a/2,-a/2,c/2);
    itf=m_Faceold.end();
    m_Faceold.insert(itf,CFace3d(v,4));

    v[0]=Vec3(-a/2,-a/2,c/2);
    v[1]=Vec3(-b/2,-b/2,-c/2);
    v[2]=Vec3(b/2,-b/2,-c/2);
    v[3]=Vec3(a/2,-a/2,c/2);
    itf=m_Faceold.end();
    m_Faceold.insert(itf,CFace3d(v,4));

    v[0]=Vec3(-a/2,a/2,c/2);
    v[1]=Vec3(-b/2,b/2,-c/2);
    v[2]=Vec3(-b/2,-b/2,-c/2);
    v[3]=Vec3(-a/2,-a/2,c/2);
    itf=m_Faceold.end();
    m_Faceold.insert(itf,CFace3d(v,4));

    v[0]=Vec3(-a/2,a/2,c/2);
    v[1]=Vec3(-b/2,b/2,-c/2);
    v[2]=Vec3(b/2,b/2,-c/2);
    v[3]=Vec3(a/2,a/2,c/2);
    itf=m_Faceold.end();
    m_Faceold.insert(itf,CFace3d(v,4));

    for(i=0;i<4;i++) {
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[i],m_veold[(i+1)%4]));
    }

    for(i=0;i<4;i++) {
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[i+4],m_veold[(i+1)%4+4]));
    }

    for(i=0;i<4;i++) {
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[i],m_veold[i+4]));
    }

    for(i=0;i<m_Nf;i++) {
        m_Faceold[i].translation(m_r);
    }
    for(i=0;i<m_Nf;i++) {
        itf=m_Face.end();
        m_Face.insert(itf,rotation_face(m_Faceold[i],m_Q,m_r));
    }
    for(i=0;i<m_Nv;i++) {
        itv=m_ve.end();
        m_ve.insert(itv,rotate(m_Q,m_veold[i])+m_r);
    }

    for(i=0;i<m_Ne;i++) {
        ite=m_Edge.end();
        m_Edgeold[i].translate(m_r);
        m_Edge.insert(ite,rotation_edge(m_Edgeold[i],m_Q,m_r));
    }
    
    wx=m_ome.X();
    wy=m_ome.Y();
    wz=m_ome.Z();
    
    m_Erot=0.5*(m_Ixx*wx*wx+m_Iyy*wy*wy+m_Izz*wz*wz);
    m_Ekin=0.5*m_m*dot(m_v,m_v);

}

CParticle::CParticle(double a,double b,int N,double R0,double l0,Vec3 &r0,Vec3 &v0,Quaternion &Q0) {
    srand(100);
    Vec3 rab,rar,r,d,I(1,0,0),J(0,1,0),K(0,0,1),Z(0,0,0);
    vector<CFace3d>::iterator itf;
    vector<CEdge3d>::iterator ite;
    vector<Vec3>::iterator itv;
    Vec3 v[4];
    int i,j;
    m_Nf=1;
    m_Ne=N+4;
    m_Nv=N+4;
    m_r=r0;
    m_R=R0;
    m_Q=Q0;
    m_large=false;
    m_dmax=sqrt((a/2)*(a/2)+(b/2)*(b/2))+R0;
    m_v=v0;
    m_ome=m_fome=Vec3(0,0,0);
    m_m=1;
    v[0]=Vec3(a/2,b/2,0);
    v[1]=Vec3(-a/2,b/2,0);
    v[2]=Vec3(-a/2,-b/2,0);
    v[3]=Vec3(a/2,-b/2,0);
    itf=m_Faceold.end();
    m_Faceold.insert(itf,CFace3d(v,4));
    i=0;
    m_Faceold[i].translation(m_r);
    itf=m_Face.end();
    m_Face.insert(itf,rotation_face(m_Faceold[i],m_Q,m_r));
    for(i=0;i<4;i++) {
        itv=m_veold.end();
        m_veold.insert(itv,v[i]);
    }
    for(i=0;i<4;i++) {
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(m_veold[i],m_veold[(i+1)%4]));
    }
    
    
    for(i=0;i<N;i++) {
        r=(a*rand()/RAND_MAX-a/2)*I+(b*rand()/RAND_MAX-b/2)*J;
        rab=r;
        rar=r+0.01*rand()/RAND_MAX*I+0.01*rand()/RAND_MAX*J+l0*K;
        ite=m_Edgeold.end();
        m_Edgeold.insert(ite,CEdge3d(rab,rar));
        itv=m_veold.end();
        m_veold.insert(itv,rar);
    }
    
    for(i=0;i<m_Nv;i++) {
        itv=m_ve.end();
        m_ve.insert(itv,rotate(m_Q,m_veold[i])+m_r);
    }
    
    for(i=0;i<m_Ne;i++) {
        m_Edgeold[i].translate(m_r);
        ite=m_Edge.end();
        m_Edge.insert(ite,rotation_edge(m_Edgeold[i],m_Q,m_r));
    }
    
    
    
    
    
}


void CParticle::rotation_quaternion(Quaternion & Q) {
    int i;
    //m_Faceold[0].m_sides[0].I().show_on_screen(); cout << "1"<<endl;
    for (i=0;i<m_Nf;i++) {
        m_Face[i]=rotation_face(m_Faceold[i],m_Q,m_r);
    }
    for (i=0;i<m_Ne;i++) {
        m_Edge[i]=rotation_edge(m_Edgeold[i],m_Q,m_r);
    }
    for (i=0;i<m_Nv;i++) {
        m_ve[i]=rotate(m_Q,m_veold[i])+m_r;
    }
       //m_Faceold[0].m_sides[0].I().show_on_screen(); cout << "2"<<endl;
    
}


void CParticle::rotation(double dt) {
    double q0,q1,q2,q3,wx,wy,wz;
    //T.show_on_screen(); cout << endl;
    q0=0.5*m_Q.getscalar();
    q1=0.5*m_Q.getvector().X();
    q2=0.5*m_Q.getvector().Y();
    q3=0.5*m_Q.getvector().Z();
    m_Td=Vec3((m_T.X()+(m_Iyy-m_Izz)*m_fome.Y()*m_fome.Z())/m_Ixx,(m_T.Y()+(m_Izz-m_Ixx)*m_fome.X()*m_fome.Z())/m_Iyy,(m_T.Z()+(m_Ixx-m_Iyy)*m_fome.X()*m_fome.Y())/m_Izz);
    //m_Td.show_on_screen(); cout << endl;
    m_ome=m_fome+0.5*dt*m_Td;

    //m_fome.show_on_screen(); cout << endl;
    //m_ome.show_on_screen(); cout << endl;

    //m_Q.show_on_screen();
    wx=m_ome.X();
    wy=m_ome.Y();
    wz=m_ome.Z();
    Quaternion dq(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz),qm;

    m_fome=m_fome+m_Td*dt;
    qm=m_Q+dq*(0.5*dt);

    q0=0.5*qm.getscalar();
    q1=0.5*qm.getvector().X();
    q2=0.5*qm.getvector().Y();
    q3=0.5*qm.getvector().Z();
    wx=m_fome.X();
    wy=m_fome.Y();
    wz=m_fome.Z();
    dq=Quaternion(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz);

    m_Q=(qm+dq*0.5*dt);
    m_Q.normalize();
    rotation_quaternion(m_Q);
    //cout << 1 << " "; m_xi0.show_on_screen(); cout<< endl;
    //cout << 2 << " "; m_xf0.show_on_screen(); cout<< endl;
    m_Erot=0.5*(m_Ixx*wx*wx+m_Iyy*wy*wy+m_Izz*wz*wz);

}

void CParticle::start(double dt) {
    m_r1=m_r-m_v*dt;
}

void CParticle::translation(double dt) {
    int i;
    Vec3 t;
    m_r2=2*m_r-m_r1+m_f*(dt*dt/m_m);
    m_v=0.5*(m_r2-m_r1)/dt;
    t=m_r2-m_r;
    m_r1=m_r;
    m_r=m_r2;
    m_Ekin=0.5*m_m*dot(m_v,m_v);
    for (i=0;i<m_Nf;i++) {
        m_Faceold[i].translation(t);
    }
    for (i=0;i<m_Ne;i++) {
        m_Edgeold[i].translate(t);
    }
    if (m_large) {
        for (i=0;i<m_Nf;i++) {
            m_Face[i].translation(t);
        }
        for (i=0;i<m_Ne;i++) {
            m_Edge[i].translate(t);
        }
        for (i=0;i<m_Nv;i++) {
            m_ve[i]+=t;
        }
    }
}





void CParticle::translation(Vec3 &t) {
    int i;
    m_r+=t;
    m_r1+=t;
    for (i=0;i<m_Nf;i++) {
        m_Faceold[i].translation(t);
    }
    for (i=0;i<m_Ne;i++) {
        m_Edgeold[i].translate(t);
    }
    rotation_quaternion(m_Q);
}

