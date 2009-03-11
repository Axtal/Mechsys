#include "Face3d.h"
CFace3d::CFace3d (Vec3 *v,int N) {
    int i;
    
    CEdge3d temp;
    m_Ns=N;
    vector<CEdge3d>::iterator it;
    it=m_sides.end();
    //m_sides=new CEdge3d[m_Ns];
    //cout<< "1" <<endl;
    for(i=0;i<N;i++) {
        temp=CEdge3d(v[i],v[(i+1)%N]);
        //v[i].show_on_screen(); cout << endl;
        //m_sides[i]=temp;
        m_sides.insert(it,temp);
        it=m_sides.end();
        //m_sides.push_back(temp);
        //it++;
    }
    //cout << (int) m_sides.size() <<endl;
}



CEdge3d distance_point_Face3d(Vec3 & v,CFace3d & F) {
    Vec3 pro,nor;
    bool inside=true;
    double s,t,a,b,c,d,f,lt,ld;
    int i;
    nor=cross(F.m_sides[0].DX(),F.m_sides[1].DX());
    nor=nor.unit();
    a=dot(F.m_sides[0].I()-v,F.m_sides[0].DX());
    b=dot(F.m_sides[0].DX(),F.m_sides[0].DX());
    c=dot(F.m_sides[0].DX(),F.m_sides[1].DX());
    d=dot(F.m_sides[0].I()-v,F.m_sides[1].DX());
    f=dot(F.m_sides[1].DX(),F.m_sides[1].DX());
    s=(c*d-a*f)/(b*f-c*c);
    t=(a*c-b*d)/(b*f-c*c);
    pro=F.m_sides[0].I()+s*F.m_sides[0].DX()+t*F.m_sides[1].DX();
    for(i=0;i<F.m_Ns;i++) {
        if (dot(cross(F.m_sides[i].DX(),pro-F.m_sides[i].I()),nor)<0) inside=false;
    }
    if (inside) return CEdge3d(v,pro);
    
    CEdge3d temp,def;
    temp=distance_point_Edge3d(v,F.m_sides[0]);
    def=temp;
    lt=temp.L();
    ld=lt;
    for(i=1;i<F.m_Ns;i++) {
        temp=distance_point_Edge3d(v,F.m_sides[i]);
        lt=temp.L();
        if(lt<ld){
            def=temp;
            ld=lt;
        }
    }
    return def;
    
    
}

void CFace3d::translation(Vec3 & a) {
    int i;
    for (i=0;i<m_Ns;i++) {
        m_sides[i].translate(a);
    }
}

CFace3d rotation_face(const CFace3d &F,Quaternion & Q,Vec3 & v) {
    int i;
    //CFace3d temp(F.m_Ns);
    CFace3d temp;
    //for (i=0;i<F.m_Ns;i++) {
      //  temp.m_sides[i]=F.m_sides[i];
    //}
    temp=F;
    //cout << temp.m_sides << " " <<F.m_sides<< endl;
    //F.m_sides[0].I().show_on_screen(); cout << "1"<<endl;
    //cout << &temp << " " <<&F<<endl;
    temp.rotation(Q,v);
    //F.m_sides[0].I().show_on_screen(); cout << "2"<<endl;
    return temp;
}

void CFace3d::rotation(Quaternion & Q,Vec3 & v) {
    int i;
    for (i=0;i<m_Ns;i++) {
        m_sides[i].rotation(Q,v);
    }
}
