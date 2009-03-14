#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<stdlib.h>
#include<map>
#include<string>
#include "../../lib/dem/Graph.h"
#define N 288

extern Quaternion normalize_rotation(double,Vec3 &);



int main () {
    srand(100000);
    int i,j,paso=1000;
    char b[10];
    double t,dt=1.e-4,A,H,W,Ftop,Fside,Epot,Work,Krot,Ktra,mu,aspect_ratio,ra,l,volume,minv,maxv,pm,p0y,p0z,py,pz,kn,kt,tcon=5,tload=10;
    Vec3 r,v,ome,I(1,0,0),J(0,1,0),K(0,0,1),Z(0,0,0);
    string s("rice");  //!rice, tetra or yermis
    ofstream BI("biaxial.txt"),FN,EN("energy.txt"),GM("granulometry.txt");
    ifstream infile("parameters.txt");
    infile >> kn; 	infile.ignore(200,'\n');
    infile >> kt; 	infile.ignore(200,'\n');
    infile >> mu;     infile.ignore(200,'\n');
    infile >> aspect_ratio;     infile.ignore(200,'\n');
    infile >> pm;     infile.ignore(200,'\n');
    infile >> dt;  infile.ignore(200,'\n');
    infile >> minv;     infile.ignore(200,'\n'); 
    infile >> maxv;     infile.ignore(200,'\n');
    paso=int(0.1/dt);
    cout << kt<< " " <<mu<< " "<<endl;
    
    Quaternion q(1,0,0,0);
    CParticle aron[N];
    CParticle wall[6];
    r=-3*K;
    wall[0]=CParticle(10,10,0.1,r,q);
    r=(1.6*(N/9)-2)*K;
    wall[5]=CParticle(10,10,0.1,r,q);
    r=2.5*J;
    q=normalize_rotation(M_PI/2,I);
    wall[1]=CParticle(10,20,0.1,r,q);
    r=-2.5*J;
    wall[2]=CParticle(10,20,0.1,r,q);
    r=2.5*I;
    q=normalize_rotation(M_PI/2,J);
    wall[3]=CParticle(20,10,0.1,r,q);
    r=-2.5*I;
    wall[4]=CParticle(20,10,0.1,r,q);
    map<pair<int,int>,CInteracton> Inaron;
    map<pair<int,int>,CInteracton> Inawall;
    map<pair<int,int>,CInteracton>::iterator it;
    for (i=0;i<N-1;i++) {
          for(j=i+1;j<N;j++) {
            pair<int,int> indexes     = make_pair(i,j);
            CInteracton   interaction = CInteracton(aron[i],aron[j]);
            Inaron.insert(make_pair(indexes,interaction));
        }
    }
    for(i=0;i<N;i++) {
        for(j=0;j<6;j++) {
            pair<int,int> indexes     = make_pair(i,j);
            CInteracton   interaction = CInteracton(aron[i],wall[j]);
            Inawall.insert(make_pair(indexes,interaction));
        }
    }
    for(j=0;j<6;j++) {
        wall[j].start(dt);
    }
    for (i=0;i<N;i++) {
	    volume=(maxv-minv)*(1.*rand())/RAND_MAX+minv;
        ra=pow(volume/(2*M_PI*(-(1./3.)+(1./aspect_ratio))),1./3.);
        l=2*(-ra+ra/aspect_ratio); 
        ome=Vec3((1.*rand())/RAND_MAX,(1.*rand())/RAND_MAX,(1.*rand())/RAND_MAX);
        q=normalize_rotation(2*M_PI*(1.*rand())/RAND_MAX,ome);
        ome=0.5*Vec3((1.*rand())/RAND_MAX,(1.*rand())/RAND_MAX,(1.*rand())/RAND_MAX);
        j=int(sqrt(9));
        r=1.6*(i%j-j/2)*I+1.6*((i/j)%j-j/2)*J+(-2+1.6*(i/9))*K;
        v=Vec3(0,0,0);
        aron[i]=CParticle(l,1.,ra,r,v,v,q,s);//length_edge, mass, sphero-radius, position, velocity, ang_velocity, quaternion, particleType 
        aron[i].start(dt);
        GM <<i<<" "<<volume<<endl;
    }
    GM.close();
    for(it=Inaron.begin();it!=Inaron.end();it++) {
         (*it).second.setmu(0.01);
	     (*it).second.setkn(kn);
	     (*it).second.setkt(kt);
    }
    for(it=Inawall.begin();it!=Inawall.end();it++) {
         (*it).second.setmu(0.01);
	     (*it).second.setkn(kn);
	     (*it).second.setkt(kt);
    }
    
    //!md loop
    Work=0;
    for(t=0,j=0;t<=310;t+=dt,j++) {
        H=wall[5].getr().Z()+3;
        W=wall[1].getr().Y()-wall[2].getr().Y();
        
        if(fabs(t-tcon)<0.1*dt) {
            p0z=fabs(wall[5].getF().Z()+wall[5].getm()*9.8)/(5*W);
            p0y=(fabs(wall[1].getF().Y())+fabs(wall[2].getF().Y()))/(10*H);
        }
        if((fabs(t-0.5*(tload+tcon))<0.1*dt)||(t>=0.5*(tload+tcon))) {
            py=pm;
            pz=pm;
        }
        else {
            py=(pm-p0y)*(t-tcon)/(0.5*(tload-tcon))+p0y;
            pz=(pm-p0z)*(t-tcon)/(0.5*(tload-tcon))+p0z;
        }
        
        if(fabs(t-tload)<0.1*dt) {
            for(it=Inaron.begin();it!=Inaron.end();it++) {
                (*it).second.setmu(mu);
            }
            
            
            for(it=Inawall.begin();it!=Inawall.end();it++) {
                (*it).second.setmu(mu);
            }
            
            q=Quaternion(1,0,0,0);
            v=(-H/1000)*K;
            cout << v<< endl;
            wall[5].setv(v);
            wall[5].start(dt);
        }
        
        if ((fabs(t-tcon)<0.1*dt)||(t>=tcon)) {
            if((fabs(t-tload)<0.1*dt)||(t>=tload)) {
                wall[5].startforce(Z);
            }
            else {
                r=-pz*W*5*K;
                wall[5].startforce(r);
            }
            r=(-py*H*5)*J;
            wall[1].startforce(r);
            r=py*H*5*J;
            wall[2].startforce(r);
            r=Z;
        }
        else {
            wall[5].startforce((-wall[5].getm()*9.8)*K);
            wall[1].startforce(Z);
            wall[2].startforce(Z);
            r=-9.8*K;
        }
        
        //cout << t << " " << wall[5].getF() << " ";
        for(i=0;i<N;i++) {
            aron[i].startforce(aron[i].getm()*r);
        }
        
        for(it=Inawall.begin();it!=Inawall.end();it++) {
            (*it).second.calcForce(dt);
        }
        for(it=Inaron.begin();it!=Inaron.end();it++) {
            (*it).second.calcForce(dt);
        }
        for(i=0;i<N;i++) {
            aron[i].translation(dt);
            aron[i].rotation(dt);
        }
        //cout << wall[5].getF()<<" "<<W<<" "<<H<<endl;
        
        
        if ((fabs(t-tcon)<0.1*dt)||(t>=tcon)) {
            
            r=wall[1].getF().Y()*J;
            wall[1].startforce(r);
            wall[1].translation(dt);
            
            r=wall[2].getF().Y()*J;
            wall[2].startforce(r);
            wall[2].translation(dt);
            
        }
        
        if ((fabs(t-tload)<0.1*dt)||(t>=tload)) {
            Ftop=wall[5].getF().Z();
            wall[5].startforce(Z);
        }
        else {
            Ftop=wall[5].getF().Z()+pz*W*5;
        }
        
        
        r=wall[5].getF().Z()*K;
        wall[5].startforce(r);
        wall[5].translation(dt);
        
        if (j%paso==0) {
            CGraph gp("poli.pov");
            for(i=0;i<N;i++) {
                gp.Graph_Particle(aron[i],"Gray");
            }
            for(i=0;i<6;i++) {
                gp.Graph_Particle(wall[i],"Blue");
            }
            gp.close();
            gp.Animation(j/paso);
            
            
            
            
            
            if ((fabs(t-tload)<0.1*dt)||(t>=tload)) BI << t <<" "<<H<<" "<<W<<" "<<Ftop<< " "<< py*H*5 << endl;
            
            
        }
        
        
    
    }
    BI.close();
    return 0;
}