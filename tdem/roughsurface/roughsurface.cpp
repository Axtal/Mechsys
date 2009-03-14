#include <cstdlib>
#include <iostream>
#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<string>
#include "../../lib/dem/Graph.h"
extern Quaternion normalize_rotation(double,Vec3 &);


using namespace std;

int main ()
{
    //inicialmente vel=0.01
    srand(100000);
    double t,dt=1.e-3,k=100,vel=0.0001,y0,s,vy,v0=0.01,dz=0.5;
    int i,j,paso=100;
    Vec3 r,v,ome,I(1,0,0),J(0,1,0),K(0,0,1),Z(0,0,0),p;
    ofstream fric("fric.txt");
    Quaternion q(1,0,0,0);
    //Recuerde den 200 radio 0.1 largo 0.2 dio muy bien
    q=normalize_rotation(M_PI,I);
    r=1*K;
    CParticle rough1(3,3,200,0.05,0.2,r,Z,q);
    q=Quaternion(1,0,0,0);
    CParticle rough2(3,3,200,0.05,0.2,Z,Z,q);
    CInteracton In(rough1,rough2);
    In.setmu(0.);
    //In.setgn(80);
    //In.setgt(40);
    rough1.start(dt);
    rough2.start(dt);
    for(t=0,j=0;t<=20000;t+=dt,j++) {
        r=-400*K;
        vy=rough1.getv().Y();
        if(j==int(2/dt)) {
            y0=rough1.getr().Y();
            cout << y0 << endl;
        }
        if(t>2) {
            s=vel*(t-2)-(rough1.getr().Y()-y0);
            r=r+k*s*J;
            //p=dz*sin(2*M_PI*(t-2)/1)*K;
            //rough2=CParticle(3,3,200,0.05,0.2,p,Z,q);
        }
        rough1.startforce(r);
        In.calcForce(dt);
        rough1.setT();
        rough1.translation(dt);
        rough1.rotation(dt);
        if (j%paso==0) {
            CGraph gp("rough.pov");
            gp.Graph_Particle(rough1,"Gray");
            gp.Graph_Particle(rough2,"Gray");
            gp.close();
            gp.Animation(j/paso);
            if(t>2) fric << t << " " << k*s << " "<<rough1.getr().Y()<< " "<<vy<< " " << vel << " "<<dz*sin(2*M_PI*(t-2)/10)<<endl;
        }
    }
}
/*
int main () {
    srand(100000);
    double t,dt=1.e-3,k=100,vel=0.01,y0,s,vy,v0=0.01;
    int i,j,paso=100;
    Vec3 r,v,ome,I(1,0,0),J(0,1,0),K(0,0,1),Z(0,0,0);
    volume=0.3;
    ra=pow(volume/(2*M_PI*(-(1./3.)+(1./aspect_ratio))),1./3.);
    l=2*(-ra+ra/aspect_ratio);
    //ofstream fric("fric.txt");
    Quaternion q(1,0,0,0);
    //Recuerde den 200 radio 0.1 largo 0.2 dio muy bien
    //q=normalize_rotation(M_PI,I);
    r=Z;
    CParticle rough1(CParticle(l,1.,ra,r,v,v,q,s);
    CGraph gp("rough.pov");
    gp.Graph_Particle(rough1,"Gray");
    gp.close();
}*/