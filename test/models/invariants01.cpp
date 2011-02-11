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
#ifdef HAS_TENSORS
  using namespace TensorsLib;
#endif

int main(int argc, char **argv) try
{
    int tst = 5;
    if (argc>1) tst = atoi(argv[1]);

    Vec_t sig;
    switch (tst)
    {
        case 1: { sig.change_dim(4); sig =   0.375,  -0.25,   0.375,  0.0*SQ2;  break; }
        case 2: { sig.change_dim(4); sig = -90.0,   -120.0, -60.0,   30.0*SQ2;  break; }
        case 3: { sig.change_dim(4); sig = -90.0,    -90.0, -90.0,    0.0*SQ2;  break; }
        case 4: { sig.change_dim(4); sig = -90.0,    -20.0, -20.0,    0.0*SQ2;  break; }
        case 5: { sig.change_dim(6); sig = -90.0,    -20.0, -20.0,    1.0*SQ2, 2.*SQ2, 1.5*SQ2;  break; }
        default: throw new Fatal("main: Test = %d is not available",tst);
    }
    size_t ncp = size(sig);

    Vec3_t L, v0,v1,v2, dpdL, dqdL;
    Vec3_t dtdL;
    Vec3_t dev_L, Ii(1.0,1.0,1.0), Lb;
    Vec_t  dev_sig(ncp), dI1(ncp),dI2(ncp),dI3(ncp), P0,P1,P2;
    double s1,s2,s3, p,q,t, p1,q1,q2,q3,t1;
    double I1,I2,I3, I1b,I2b,I3b;
    double s1b,s2b,s3b;
    EigenProj (sig, L, v0,v1,v2, P0,P1,P2);
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

    // derivative of invariants
    Vec3_t S((2.*L(0)-L(1)-L(2))/3., 
             (2.*L(1)-L(2)-L(0))/3.,
             (2.*L(2)-L(0)-L(1))/3.);
    Vec3_t devSS((2.*S(1)*S(2) - S(2)*S(0) - S(0)*S(1))/3.,
                 (2.*S(2)*S(0) - S(1)*S(2) - S(0)*S(1))/3.,
                 (2.*S(0)*S(1) - S(1)*S(2) - S(2)*S(0))/3.);
    Vec3_t dtdL2;
    dtdL2 = (-3./pow(q,3.))*(q*t*S + SQ6*devSS);

    double p4,q4,t4, p5,q5,t5;
    Mat3_t dpqthdL, dLdpqth;
    Vec3_t Lsorted(L);
    Util::Sort   (Lsorted(0), Lsorted(1), Lsorted(2));
    OctDerivs    (Lsorted, p4,q4,t4, dpqthdL);
    InvOctDerivs (Lsorted, p5,q5,t5, dLdpqth);
    cout << "dpqthdL = \n" << PrintMatrix(dpqthdL);
    cout << "dLdpqth = \n" << PrintMatrix(dLdpqth);
    Mat3_t res;
    res = product(dpqthdL,dLdpqth);
    cout << "dpqthdL*dLdpqth = \n" << PrintMatrix(res);
    bool err_deriv = false;
    for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
    {
        if (i!=j && fabs(res(i,j)>1.0e-10)) err_deriv = true;
    }
    if (err_deriv) cout << TERM_RED << "Error: dpqthdL*dLdpqth is not identity" << TERM_RST << endl;

    // inverse
    Vec_t invSig;
    Mat3_t m_sig,m_invSig;
    Inv     (sig,    invSig);
    Ten2Mat (sig,    m_sig);
    Ten2Mat (invSig, m_invSig);
    res = product(m_sig,m_invSig);
    cout << "\nsig*invSig = \n" << PrintMatrix(res);
    bool err_inv = false;
    for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
    {
        if (i!=j && fabs(res(i,j)>1.0e-15)) err_inv = true;
    }
    if (err_inv) cout << TERM_RED << "Error: inverse failed" << TERM_RST << endl;

    // tensors
    bool err_zer   = false;
    bool err_t2ten = false;
    bool err_ten2t = false;
    bool err_invt  = false;
    bool err_mat   = false;
#ifdef HAS_TENSORS
    // derivative of inverse
    if (ncp==6)
    {
        Vec_t             Atmp;
        Tensor2<double,3> A, iA;
        Ten2Tensor (sig, A);
        Tensor2Ten (A, Atmp, ncp);
        Inv        (A, iA);
        cout << "\nsig ="  << PrintVector(sig);
        cout << "sig ="    << PrintVector(Atmp);
        cout << "sig =\n"  << A << endl;
        cout << "iSig ="   << PrintVector(invSig);
        cout << "iSig =\n" << iA << endl;
        for (size_t i=0; i<3; ++i)
        for (size_t j=0; j<3; ++j)
        {
            int    k = Tensor2ToMandel_i[i][j];
            double c = (i==j ? 1.0 : SQ2);
            if (k>static_cast<int>(ncp-1) && i!=j) // 2D
            {
                // check if off diagonal components are zero
                if (fabs(A[i][j])>1.0e-15) err_zer = true;
            }
            else
            {
                if (fabs(A[i][j]*c - sig(k))>1.0e-15) err_t2ten = true;
            }
        }
        for (size_t i=0; i<ncp; ++i)
        {
            int    k = MandelToTensor2_i[i];
            int    l = MandelToTensor2_j[i];
            double c = (k==l ? 1.0 : SQ2);
            if (fabs(sig(i)    - Atmp(i)   )>1.0e-15) err_ten2t = true;
            if (fabs(invSig(i) - iA[k][l]*c)>1.0e-15) err_invt  = true;
        }
        if (err_zer)   cout << TERM_RED << "Error: Off diagonal components of A in 2D must be zero"   << TERM_RST << endl;
        if (err_t2ten) cout << TERM_RED << "Error: Ten2Tensor failed"                                 << TERM_RST << endl;
        if (err_ten2t) cout << TERM_RED << "Error: Tensor2Ten failed"                                 << TERM_RST << endl;
        if (err_invt)  cout << TERM_RED << "Error: Inverse of Tensor2 is not equal to inverse of Ten" << TERM_RST << endl;
        Mat_t             diSigdSig, diAdAmat;
        Tensor4<double,3> diAdA;
        DerivInv (sig, invSig, diSigdSig);
        diAdA = ((iA^iA) + (iA|iA)) * (-0.5);
        Tensor2Ten (diAdA, diAdAmat, ncp);
        cout << endl;
        cout << "diSigdSig =\n" << PrintMatrix(diSigdSig) << endl;
        cout << "diSigdSig =\n" << PrintMatrix(diAdAmat)  << endl;
        if (CompareMatrices(diSigdSig, diAdAmat)>1.0e-15) err_mat = true;
        if (err_mat) cout << TERM_RED << "Error: diSigdSig is not equal to diAdAmat" << TERM_RST << endl;
        /*
        Tensor4<double,3> One;
        for (size_t i=0; i<3; ++i)
        for (size_t j=0; j<3; ++j)
        for (size_t k=0; k<3; ++k)
        for (size_t l=0; l<3; ++l)
            One[i][j][k][l] = 1.0;+i+j+k+l;
        Mat_t one;
        Tensor2Ten (One, one, 6);
        cout << PrintMatrix(one) << endl;
        */
    }
#endif

    printf("\n");
    printf("sig                = [%g, %g, %g, %g]  \n",sig(0),sig(1),sig(2),sig(3));
    printf("dev_sig            = [%g, %g, %g, %g]  \n",dev_sig(0),dev_sig(1),dev_sig(2),dev_sig(3));
    printf("theta              = %g                \n",60.0*asin(t)/PI);
    printf("p,  q,  t          = %g, %g, %g        \n",p,  q,  t);
    printf("p1, q1, t1, q2, q3 = %g, %g, %g, %g, %g\n",p1, q1, t1,   q2, q3);
    printf("s1,  s2,  s3       = %g, %g, %g        \n",s1,  s2,  s3);
    printf("s1b, s2b, s3b      = %g, %g, %g        \n",s1b, s2b, s3b);
    printf("I1 , I2,  I3       = %g, %g, %g        \n",I1 , I2 , I3);
    printf("I1b, I2b, I3b      = %g, %g, %g        \n",I1b, I2b, I3b);
    printf("dI1dsig            = [%g, %g, %g, %g]  \n",dI1(0),dI1(1),dI1(2),dI1(3));
    printf("dI2dsig            = [%g, %g, %g, %g]  \n",dI2(0),dI2(1),dI2(2),dI2(3));
    printf("dI3dsig            = [%g, %g, %g, %g]  \n",dI3(0),dI3(1),dI3(2),dI3(3));
    printf("dpdL               = [%g, %g, %g]      \n",dpdL (0),dpdL (1),dpdL (2));
    printf("dqdL               = [%g, %g, %g]      \n",dqdL (0),dqdL (1),dqdL (2));
    printf("dtdL               = [%g, %g, %g]      \n",dtdL (0),dtdL (1),dtdL (2));
    printf("dtdL2              = [%g, %g, %g]      \n",dtdL2(0),dtdL2(1),dtdL2(2));
    printf("dev_L              = [%g, %g, %g]      \n",dev_L(0),dev_L(1),dev_L(2));
    printf("S                  = [%g, %g, %g]      \n",S    (0),S    (1),S    (2));

    Util::Sort(s1,s2,s3);
    Util::Sort(s1b,s2b,s3b);
    Array<double> err(11);
    err = fabs(p-p1  ),  fabs(q-q1  ),  fabs(t-t1  ),  fabs(q-q2),  fabs(q-q3),
          fabs(s1-s1b),  fabs(s2-s2b),  fabs(s3-s3b),
          fabs(I1-I1b),  fabs(I2-I2b),  fabs(I3-I3b);

    printf("ERROR(p,q,t)       = %g, %g, %g, %g, %g\n",err[0], err[1], err[2], err[3], err[4]);
    printf("ERROR(sk-skb)      = %g, %g, %g        \n",err[5], err[6], err[7]);
    printf("ERROR(Ik)          = %g, %g, %g        \n",err[8], err[9], err[10]);

    printf("\n ----------------------------- Model functions --------------------------------\n");
    double E   = 3000.0;
    double nu  = 0.25;
    double K   = Calc_K (E, nu);
    double G   = Calc_G (E, nu);
    double G1  = Calc_G_(K, nu);
    double E1  = Calc_E (K, G);
    double E2  = Calc_E_(K, nu);
    double nu1 = Calc_nu (K, G);
    printf("E  = %g,   E1  = %g,   E2 = %g\n", E,  E1, E2);
    printf("nu = %g,   nu1 = %g\n", nu, nu1);
    printf("K  = %g\n", K);
    printf("G  = %g,   G1  = %g\n", G, G1);
    bool err_E = (fabs(E-E1)>0.0) || (fabs(E-E2)>0.0);
    bool err_nu = (fabs(nu-nu1)>0.0);
    bool err_G  = (fabs(G-G1)>0.0);

    if (err_G)               return 1;
    if (err_E)               return 1;
    if (err_nu)              return 1;
    if (err_deriv)           return 1;
    if (err_inv)             return 1;
    if (err_zer)             return 1;
    if (err_t2ten)           return 1;
    if (err_ten2t)           return 1;
    if (err_invt)            return 1;
    if (err_mat)             return 1;
    if (err.TheMax()>1.0e-8) return 1;
    else                     return 0;
}
MECHSYS_CATCH
