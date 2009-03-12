#include "Edge3d.h"

CEdge3d::CEdge3d (const Vec3 & a,const Vec3 & b) {
    m_xi=a;
    m_xf=b;
    m_dx=m_xf-m_xi;
    m_lenght=m_dx.norm();   
}

CEdge3d::CEdge3d (void) {
    
}

CEdge3d distance_Edge3d_Edge3d(const CEdge3d & E0,const CEdge3d & E1) {
    bool ex=false;
    double t,s,a,b,c,d,e;
    Vec3 I,F,os;
    CEdge3d temp,l1,l2,l3,l4;
    c=dot(E0.m_dx,E0.m_dx);
    e=dot(E0.m_dx,E1.m_dx);
    d=dot(E1.m_dx,E1.m_dx);
    a=dot(E0.m_dx,E0.m_xi-E1.m_xi);
    b=dot(E1.m_dx,E0.m_xi-E1.m_xi);
    t=(c*b-e*a)/(c*d-e*e);
    s=(e*b-a*d)/(c*d-e*e);
    //cout <<s<<" "<<t<< endl;
    
    if ((s>0)&&(s<1)&&(t>0)&&(t<1)) {
        I=E0.m_xi+E0.m_dx*s;
        F=E1.m_xi+E1.m_dx*t;
        //I.show_on_screen(); cout<< endl;
        //F.show_on_screen(); cout<< endl;
        //os.show_on_screen(); cout<< endl;
        
        temp=CEdge3d(I,F);
        return temp;
    }
    else {
        l1=distance_point_Edge3d(E0.m_xi,E1);
        l2=distance_point_Edge3d(E0.m_xf,E1);
        l3=distance_point_Edge3d(E1.m_xi,E0);
        l4=distance_point_Edge3d(E1.m_xf,E0);
        l3.invert();
        l4.invert();
        if((l1.L()<=l2.L())&&(l1.L()<=l3.L())&&(l1.L()<=l4.L())) return l1;
        if((l2.L()<=l1.L())&&(l2.L()<=l3.L())&&(l2.L()<=l4.L())) return l2;
        if((l3.L()<=l1.L())&&(l3.L()<=l2.L())&&(l3.L()<=l4.L())) return l3;
        if((l4.L()<=l1.L())&&(l4.L()<=l2.L())&&(l4.L()<=l3.L())) return l4;
    }
    
    
}

void CEdge3d::rotation (Quaternion & Q,Vec3 & v) {
    Vec3 temp;
    temp=m_xi-v;
    m_xi=rotate(Q,temp)+v;
    temp=m_xf-v;
    m_xf=rotate(Q,temp)+v;
    m_dx=m_xf-m_xi;
    
}

CEdge3d distance_point_Edge3d(const Vec3 &V,const CEdge3d &E) {
    double t;
    Vec3 I,F;
    t=(dot(V,E.m_dx)-dot(E.m_xi,E.m_dx))/(dot(E.m_dx,E.m_dx));
    I=V;
    if (t<0) F=E.m_xi;
    else if (t>1) F=E.m_xf;
    else F=E.m_xi+E.m_dx*t;
    CEdge3d temp(I,F);
    return temp;
}

CEdge3d rotation_edge(CEdge3d & E,Quaternion & Q,Vec3 & v) {
    int i;
    CEdge3d temp;
    temp=E;
    temp.rotation(Q,v);
    return temp;
}
