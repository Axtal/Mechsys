#include "Interacton.h"

CInteracton::CInteracton(CParticle & p1,CParticle & p2)
{
    m_P1=&p1;
    m_P2=&p2;
    m_kn=10000.;
    m_kt=5000.;
    m_gn=16;
    m_gt=8;
    m_mu=0.4;
    //m_gn=0;
    //m_gt=0;
    m_Epot=0;
    m_c=false;
}


void CInteracton::calcForce(double dt) {
    double delta,R1,R2,dist,alpha;
    CEdge3d F;
    Vec3 f,ft,n,r,T,vrel,vt,t;
    Quaternion q;
    R1=(m_P1->m_R);
    R2=(m_P2->m_R);
    int i,j;
    m_Epot=0;
    m_c=false;
    pair<int,int> p;
    
    
    m_fn=Vec3(0,0,0);
    r=m_P1->m_r-m_P2->m_r;
    vrel=m_P2->m_v-m_P1->m_v;
    //vrel.show_on_screen(); cout << endl;
    dist=r.norm();
    delta=m_P1->m_dmax+m_P2->m_dmax;
    
    if(dist<=delta) {
        for (i=0;i<m_P1->m_Ne;i++) {
            for (j=0;j<m_P2->m_Ne;j++) {
                F=distance_Edge3d_Edge3d(m_P1->m_Edge[i],m_P2->m_Edge[j]);
                dist=F.DX().norm();
                delta=R1+R2-dist;
                if (delta>0) {
                    
                    m_c=true;
                    n=F.DX().unit();
                    f=n*(m_kn*delta);
                    m_fn+=f;
                    //f.show_on_screen(); cout << i << " " <<j<<" "<< delta<<endl;
                    r=F.I()+n*((R1*R1-R2*R2+dist*dist)/(2*dist));
                    alpha=1;
                    //alpha=(R1+R2)/dist;
                    
                    vrel=-(alpha*(m_P2->m_v-m_P1->m_v)+cross(rotate(m_P2->m_Q,m_P2->m_ome),r-m_P2->m_r)-cross(rotate(m_P1->m_Q,m_P1->m_ome),r-m_P1->m_r));
                    vt=vrel-dot(n,vrel)*n;
                    p=make_pair(i,j);
                    
                    m_free[p]+=vt*dt;
                    m_free[p]-=dot(m_free[p],n)*n;
                    t=m_free[p].unit();
                    if(m_free[p].norm()>m_mu*f.norm()/m_kt) {
                        m_free[p]=t*m_mu*f.norm()/m_kt;
                    }
                    ft=m_kt*m_free[p];
                    //cout<< m_mu<<" "<<m_free[p]<<endl;
                    m_Epot+=0.5*m_kn*delta*delta+0.5*m_kt*m_free[p]*m_free[p];
                    
                    
                    m_P1->m_f-=f+m_gn*dot(n,vrel)*n+m_gt*vt+ft;
                    m_P2->m_f+=f+m_gn*dot(n,vrel)*n+m_gt*vt+ft;
                
                    T=cross(r-m_P1->m_r,f+m_gn*dot(n,vrel)*n+m_gt*vt+ft);
                    q=m_P1->m_Q.conj();
                    m_P1->m_T-=rotate(q,T);
                
                    T=cross(r-m_P2->m_r,f+m_gn*dot(n,vrel)*n+m_gt*vt+ft);
                    q=m_P2->m_Q.conj();
                    m_P2->m_T+=rotate(q,T);
                
                }
            }
        }
    
        for (i=0;i<m_P1->m_Nv;i++) {
            for (j=0;j<m_P2->m_Nf;j++) {
                F=distance_point_Face3d(m_P1->m_ve[i],m_P2->m_Face[j]);
                dist=F.DX().norm();
                delta=R1+R2-dist;
                if (delta>0) {
                    
                    m_c=true;
                    n=F.DX().unit();
                    f=n*(m_kn*delta);
                    m_fn+=f;
                    //f.show_on_screen(); cout << i << " " <<j<<" "<< delta<<endl;
                    r=F.I()+n*((R1*R1-R2*R2+dist*dist)/(2*dist));
                    alpha=1;
                    //alpha=(R1+R2)/dist;

                    vrel=-(alpha*(m_P2->m_v-m_P1->m_v)+cross(rotate(m_P2->m_Q,m_P2->m_ome),r-m_P2->m_r)-cross(rotate(m_P1->m_Q,m_P1->m_ome),r-m_P1->m_r));
                    vt=vrel-dot(n,vrel)*n;
                    p=make_pair(i,j);
                    
                    m_frvf[p]+=vt*dt;
                    m_frvf[p]-=dot(m_frvf[p],n)*n;
                    t=m_frvf[p].unit();
                    if(m_frvf[p].norm()>m_mu*f.norm()/m_kt) {
                        m_frvf[p]=t*m_mu*f.norm()/m_kt;
                    }
                    ft=m_kt*m_frvf[p];
                    
                    m_Epot+=0.5*m_kn*delta*delta+0.5*m_kt*m_frvf[p]*m_frvf[p];
                    
                    m_P1->m_f-=f+m_gn*dot(n,vrel)*n+m_gt*vt+ft;
                    m_P2->m_f+=f+m_gn*dot(n,vrel)*n+m_gt*vt+ft;

                    T=cross(r-m_P1->m_r,f+m_gn*dot(n,vrel)*n+m_gt*vt+ft);
                    q=m_P1->m_Q.conj();
                    m_P1->m_T-=rotate(q,T);

                    T=cross(r-m_P2->m_r,f+m_gn*dot(n,vrel)*n+m_gt*vt+ft);
                    q=m_P2->m_Q.conj();
                    m_P2->m_T+=rotate(q,T);

                }
            }
        }
    
        for (i=0;i<m_P1->m_Nf;i++) {
            for (j=0;j<m_P2->m_Nv;j++) {
                F=distance_point_Face3d(m_P2->m_ve[j],m_P1->m_Face[i]);
                F.invert();
                dist=F.DX().norm();
                delta=R1+R2-dist;
                //f.show_on_screen(); cout << i << " " <<j<<" "<< delta<<endl;
                if (delta>0) {
                    m_c=true;
                    n=F.DX().unit();
                    f=n*(m_kn*delta);
                    m_fn+=f;
                    //if(delta>R1) f.show_on_screen(); cout << i << " " <<j<<" "<< delta<<endl;
                    //m_P2->m_ve[i].show_on_screen(); cout << i << " " <<j<<" "<< delta<<endl;
                    r=F.I()+n*((R1*R1-R2*R2+dist*dist)/(2*dist));
                    alpha=1;
                    //alpha=(R1+R2)/dist;

                    vrel=-(alpha*(m_P2->m_v-m_P1->m_v)+cross(rotate(m_P2->m_Q,m_P2->m_ome),r-m_P2->m_r)-cross(rotate(m_P1->m_Q,m_P1->m_ome),r-m_P1->m_r));
                    vt=vrel-dot(n,vrel)*n;
                    p=make_pair(i,j);
                    
                    m_frfv[p]+=vt*dt;
                    m_frfv[p]-=dot(m_frfv[p],n)*n;
                    t=m_frfv[p].unit();
                    if(m_frfv[p].norm()>m_mu*f.norm()/m_kt) {
                        m_frfv[p]=t*m_mu*f.norm()/m_kt;
                    }
                    ft=m_kt*m_frfv[p];
                    
                    m_Epot+=0.5*m_kn*delta*delta+0.5*m_kt*m_frfv[p]*m_frfv[p];

                    m_P1->m_f-=f+m_gn*dot(n,vrel)*n+m_gt*vt+ft;
                    m_P2->m_f+=f+m_gn*dot(n,vrel)*n+m_gt*vt+ft;

                    T=cross(r-m_P1->m_r,f+m_gn*dot(n,vrel)*n+m_gt*vt+ft);
                    q=m_P1->m_Q.conj();
                    m_P1->m_T-=rotate(q,T);

                    T=cross(r-m_P2->m_r,f+m_gn*dot(n,vrel)*n+m_gt*vt+ft);
                    q=m_P2->m_Q.conj();
                    m_P2->m_T+=rotate(q,T);


                }
            }
        }
    }    
}
