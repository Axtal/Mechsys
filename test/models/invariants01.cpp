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

// STL
#include <iostream>
#include <cmath>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;
using Util::SQ2;
using Util::SQ6;
using Util::PI;

int main(int argc, char **argv) try
{
    int tst = 1;
    if (argc>1) tst = atoi(argv[1]);

    Vec_t sig(4);
    switch (tst)
    {
        case 1: { sig =   0.375,  -0.25,   0.375,  0.0*SQ2;  break; }
        case 2: { sig = -90.0,   -120.0, -60.0,   30.0*SQ2;  break; }
        case 3: { sig = -90.0,    -90.0, -90.0,    0.0*SQ2;  break; }
        case 4: { sig = -90.0,    -20.0, -20.0,    0.0*SQ2;  break; }
        default: throw new Fatal("main: Test = %d is not available",tst);
    }

    Vec3_t L, dpdL, dqdL, dtdL, dev_L, Ii(1.0,1.0,1.0), Lb;
    Vec_t  dev_sig(4), dI1(4),dI2(4),dI3(4), P0,P1,P2;
    double s1,s2,s3, p,q,t, p1,q1,q2,q3,t1;
    double I1,I2,I3, I1b,I2b,I3b;
    double s1b,s2b,s3b;
    EigenProj (sig, L, P0,P1,P2);
    OctInvs   (sig, p, q, t);
    OctInvs   (L,   p1,q1,t1, dpdL,dqdL,dtdL);
    CharInvs  (sig, I1,I2,I3, dI1,dI2,dI3);
    pqth2L    (p,q, asin(t)/3.0, Lb);
    Dev       (sig, dev_sig);

    s1 =L (0);  s2 =L (1);  s3= L (2);
    s1b=Lb(0);  s2b=Lb(1);  s3b=Lb(2);

    dev_L = L - ((L(0)+L(1)+L(2))/3.0)*Ii;
    q2    = Norm(dev_L);
    q3    = Norm(dev_sig);
    I1b   = s1+s2+s3;
    I2b   = s1*s2+s2*s3+s3*s1;
    I3b   = s1*s2*s3;

    String buf;
    buf.Printf("sig                = [%g, %g, %g, %g]  \n",sig(0),sig(1),sig(2),sig(3));                 cout<<buf;
    buf.Printf("dev_sig            = [%g, %g, %g, %g]  \n",dev_sig(0),dev_sig(1),dev_sig(2),dev_sig(3)); cout<<buf;
    buf.Printf("theta              = %g                \n",60.0*asin(t)/PI);                             cout<<buf;
    buf.Printf("p,  q,  t          = %g, %g, %g        \n",p,  q,  t);                                   cout<<buf;
    buf.Printf("p1, q1, t1, q2, q3 = %g, %g, %g, %g, %g\n",p1, q1, t1,   q2, q3);                        cout<<buf;
    buf.Printf("s1,  s2,  s3       = %g, %g, %g        \n",s1,  s2,  s3);                                cout<<buf;
    buf.Printf("s1b, s2b, s3b      = %g, %g, %g        \n",s1b, s2b, s3b);                               cout<<buf;
    buf.Printf("I1 , I2,  I3       = %g, %g, %g        \n",I1 , I2 , I3);                                cout<<buf;
    buf.Printf("I1b, I2b, I3b      = %g, %g, %g        \n",I1b, I2b, I3b);                               cout<<buf;
    buf.Printf("dI1dsig            = [%g, %g, %g, %g]  \n",dI1(0),dI1(1),dI1(2),dI1(3));                 cout<<buf;
    buf.Printf("dI2dsig            = [%g, %g, %g, %g]  \n",dI2(0),dI2(1),dI2(2),dI2(3));                 cout<<buf;
    buf.Printf("dI3dsig            = [%g, %g, %g, %g]  \n",dI3(0),dI3(1),dI3(2),dI3(3));                 cout<<buf;

    Util::Sort(s1,s2,s3);
    Util::Sort(s1b,s2b,s3b);
    Array<double> err(11);
    err = fabs(p-p1  ),  fabs(q-q1  ),  fabs(t-t1  ),  fabs(q-q2),  fabs(q-q3),
          fabs(s1-s1b),  fabs(s2-s2b),  fabs(s3-s3b),
          fabs(I1-I1b),  fabs(I2-I2b),  fabs(I3-I3b);

    buf.Printf("ERROR(p,q,t)       = %g, %g, %g, %g, %g\n",err[0], err[1], err[2], err[3], err[4]);      cout<<buf;
    buf.Printf("ERROR(sk-skb)      = %g, %g, %g        \n",err[5], err[6], err[7]);                      cout<<buf;
    buf.Printf("ERROR(Ik)          = %g, %g, %g        \n",err[8], err[9], err[10]);                     cout<<buf;

    if (err.TheMax()>1.0e-8) return 1;
    else                     return 0;
}
MECHSYS_CATCH
