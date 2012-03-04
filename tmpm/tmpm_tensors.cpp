/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso                                *
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

/* ttensors - Copyright (C) 2007 Dorival de Moraes Pedroso */

// STL
#include <iostream>
#include <cmath> // for fabs

// Local
#include <mechsys/mpm/defs.h>
#include <mechsys/mpm/fmtnum.h>
#include <mechsys/mpm/tensors.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using MPM::Vector3D;
using MPM::STensor2;
using MPM::STensor4;
using MPM::ATensor2;
using MPM::ATensor4;
using MPM::SQ2;

double Error (Vector3D const & v, Vector3D const & c)
{
	double err = 0.0;
	for (int i=0; i<3; ++i)
		err += fabs(v(i)-c(i));
	return err;
}

double Error (STensor2 const & T, STensor2 const & C)
{
	double err = 0.0;
	for (int i=0; i<6; ++i)
		err += fabs(T(i)-C(i));
	return err;
}

double Error (STensor4 const & T, STensor4 const & C)
{
	double err = 0.0;
	for (int i=0; i<6; ++i)
	for (int j=0; j<6; ++j)
		err += fabs(T(i,j)-C(i,j));
	return err;
}

double Error (ATensor2 const & T, ATensor2 const & C)
{
	double err = 0.0;
	for (int i=0; i<9; ++i)
		err += fabs(T(i)-C(i));
	return err;
}

int main(int argc, char **argv) try
{
	cout << "\n ~~~~~~~~~~~ Operations ~~~~~~~~~~~ " << endl;

/** Add scaled tensors Z = a*X + b*Y. */
	{
		STensor4 X,Y,Z,cZ;
		X = 0.,1.,2.,3.,4.,5.,
		    5.,4.,3.,2.,1.,0.,
		    0.,1.,2.,3.,4.,5.,
		    5.,4.,3.,2.,1.,0.,
		    0.,1.,2.,3.,4.,5.,
		    5.,4.,3.,2.,1.,0.;
		Y = X;
		cZ =  0.0 ,  3.0 , 6.0 , 9.0 , 12.0 , 15.0, 
		     15.0 , 12.0 , 9.0 , 6.0 ,  3.0 ,  0.0, 
		      0.0 ,  3.0 , 6.0 , 9.0 , 12.0 , 15.0, 
		     15.0 , 12.0 , 9.0 , 6.0 ,  3.0 ,  0.0, 
		      0.0 ,  3.0 , 6.0 , 9.0 , 12.0 , 15.0, 
		     15.0 , 12.0 , 9.0 , 6.0 ,  3.0 ,  0.0;
        MPM::AddScaled (1.0, X, 2.0, Y,  Z);
		cout << "AddScaled    error = " << Error(Z,cZ) << endl;
	}

/** Symmetric of a tensor S = a*sym(T) = a * 0.5*(T + trn(T)). */
	{
		ATensor2 T;
		STensor2 S,cS;
		T  = 1., 2., 3., 4. , 5. , 6. , 7. , 8. , 9.;
		cS = 1., 2., 3., 5.5*SQ2, 6.5*SQ2, 7.5*SQ2;
        MPM::Sym (1.0, T,  S);
		cout << "Sym          error = " << Error(S,cS) << endl;
	}

/** Determinant. */
	{
		ATensor2 T;
		T  = 1., 2., 3., 4. , 5. , 6. , 7. , 8. , 9.;
		cout << "Det          error = " << fabs(MPM::Det(T)-290.0) << endl;
	}

	cout << "\n ~~~~~~~~~~~ Asymmetric ~~~~~~~~~~~ " << endl;

//////////// Only assignment

/** 1) Dot product: v=A*u  =>  v(i)=A(i,k)*u(k) */
	{
		ATensor2 A;
		Vector3D u;
		Vector3D v,cv;
		A =  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.;
		u =  1. , 2. , 3.;
		cv = 27., 26., 34.;
        MPM::Dot (A, u,  v);
		cout << "1) Dot       error = " << Error(v,cv) << endl;
	}

/** 2) Dot product: v=u*A  =>  v(i)=u(k)*A(k,i) */
	{
		ATensor2 A;
		Vector3D u;
		Vector3D v,cv;
		A =  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.;
		u =  1. , 2. , 3.;
		cv = 42., 32., 25.;
        MPM::Dot (u, A,  v);
		cout << "2) Dot       error = " << Error(v,cv) << endl;
	}

/** 3) Dot product: y = T*x  =>  y = T:x  =>  y(i,j)=T(i,j,k,l)*x(k,l) */
	{
		ATensor4 T;
		ATensor2 x,y,cy;
		T =  1. ,  2. ,  3. ,  4. ,  5. ,  6. ,  1. ,  2. ,  3. ,
		     6. ,  5. ,  4. ,  3. ,  2. ,  1. ,  6. ,  5. ,  4. ,
		     7. ,  8. ,  9. , 10. , 11. , 12. ,  7. ,  8. ,  9. ,
		    12. , 11. , 10. ,  9. ,  8. ,  7. , 12. , 11. , 10. ,
		     1. ,  2. ,  3. ,  4. ,  5. ,  6. ,  1. ,  2. ,  3. ,
		     6. ,  5. ,  4. ,  3. ,  2. ,  1. ,  6. ,  5. ,  4. ,
		     7. ,  8. ,  9. , 10. , 11. , 12. ,  7. ,  8. ,  9. ,
		    12. , 11. , 10. ,  9. ,  8. ,  7. , 12. , 11. , 10. ,
		     1. ,  2. ,  3. ,  4. ,  5. ,  6. ,  1. ,  2. ,  3. ;
		x =  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.;
		cy = 141.,174.,411.,444.,141.,174.,411.,444.,141.;
        MPM::Dot (T, x,  y);
		cout << "3) Dot       error = " << Error(y,cy) << endl;
	}

/** 4) Dot product: C=A*B  =>  C(i,j)=A(i,k)*B(k,j) */
	{
		ATensor2 A,B,C,cC;
		A =  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.;
		B =  9. , 8. , 7. , 6. , 5. , 4. , 3. , 2. , 1.;
		cC = 27., 68., 97., 50., 73., 66., 74., 124., 108.;
        MPM::Dot (A, B,  C);
		cout << "4) Dot       error = " << Error(C,cC) << endl;
	}

/** a) Dyadic product: A = u dyad v  =>  A(i,j)=u(i)*v(j) */
	{
		Vector3D u,v;
		ATensor2 A,cA;
		u = 1., 2., 3.;
		v = 4., 5., 6.;
		cA = 4., 10., 18., 5., 12., 6., 8., 15., 12.;
        MPM::Dyad (u, v,  A);
		cout << "a) Dyad      error = " << Error(A,cA) << endl;
	}

/** b) Dyadic with Dot: R=v dyad (F*u)  =>  R(i,j)=v(i)*F(j,k)*u(k) */
	{
		Vector3D v,u;
		ATensor2 F,R,cR;
		u = 1., 2., 3.;
		v = 4., 5., 6.;
		F = 4., 10., 18., 5., 12., 6., 8., 15., 12.;
		cR = 128., 320., 576., 256., 480., 384., 160., 384., 192.;
        MPM::DyadDot (v, F, u,  R);
		cout << "b) DyadDot   error = " << Error(R,cR) << endl;
	}

/** c) Dyadic with Dot: R=(v dyad u)*F  =>  R(i,j)=v(i)*u(k)*F(k,j) */
	{
		Vector3D v,u;
		ATensor2 F,R,cR;
		u = 1., 2., 3.;
		v = 4., 5., 6.;
		F = 4., 10., 18., 5., 12., 6., 8., 15., 12.;
		cR = 224., 350., 504., 280., 420., 336., 280., 420., 336.;
        MPM::DyadDot (v, u, F,  R);
		cout << "c) DyadDot   error = " << Error(R,cR) << endl;
	}

//////////////// With update

/** 4) Dot product (with update): v += A*u  =>  v(i)+=A(i,k)*u(k) */
	{
		ATensor2 A;
		Vector3D u;
		Vector3D v,cv;
		A =  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.;
		u =  1. , 2. , 3.;
		v = 0.1;
		cv = 27.1, 26.1, 34.1;
        MPM::DotUp (A, u,  v);
		cout << "4) DotUp     error = " << Error(v,cv) << endl;
	}
 
/** 5) Dot product (with update): v += u*A  =>  v(i)+=u(k)*A(k,i) */
	{
		ATensor2 A;
		Vector3D u;
		Vector3D v,cv;
		A =  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.;
		u =  1. , 2. , 3.;
		cv = 42.1, 32.1, 25.1;
		v = 0.1;
        MPM::DotUp (u, A,  v);
		cout << "5) DotUp     error = " << Error(v,cv) << endl;
	}

/** 6) Dot product (with update): y += T*x  =>  y += T:x  =>  y(i,j)+=T(i,j,k,l)*x(k,l) */
	{
		ATensor4 T;
		ATensor2 x,y,cy;
		T =  1. ,  2. ,  3. ,  4. ,  5. ,  6. ,  1. ,  2. ,  3. ,
		     6. ,  5. ,  4. ,  3. ,  2. ,  1. ,  6. ,  5. ,  4. ,
		     7. ,  8. ,  9. , 10. , 11. , 12. ,  7. ,  8. ,  9. ,
		    12. , 11. , 10. ,  9. ,  8. ,  7. , 12. , 11. , 10. ,
		     1. ,  2. ,  3. ,  4. ,  5. ,  6. ,  1. ,  2. ,  3. ,
		     6. ,  5. ,  4. ,  3. ,  2. ,  1. ,  6. ,  5. ,  4. ,
		     7. ,  8. ,  9. , 10. , 11. , 12. ,  7. ,  8. ,  9. ,
		    12. , 11. , 10. ,  9. ,  8. ,  7. , 12. , 11. , 10. ,
		     1. ,  2. ,  3. ,  4. ,  5. ,  6. ,  1. ,  2. ,  3. ;
		x =  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.;
		cy = 141.1,174.1,411.1,444.1,141.1,174.1,411.1,444.1,141.1;
		y = 0.1;
        MPM::DotUp (T, x,  y);
		cout << "6) DotUp     error = " << Error(y,cy) << endl;
	}

/** d) Dyadic product (with update): A += u dyad v  =>  A(i,j)+=u(i)*v(j) */
	{
		Vector3D u,v;
		ATensor2 A,cA;
		u = 1., 2., 3.;
		v = 4., 5., 6.;
		cA = 4.1, 10.1, 18.1, 5.1, 12.1, 6.1, 8.1, 15.1, 12.1;
		A = 0.1;
        MPM::DyadUp (u, v,  A);
		cout << "d) DyadUp    error = " << Error(A,cA) << endl;
	}

/** e) Dyadic with Dot (with update): R += v dyad (F*u)  =>  R(i,j)+=v(i)*F(j,k)*u(k) */
	{
		Vector3D v,u;
		ATensor2 F,R,cR;
		u = 1., 2., 3.;
		v = 4., 5., 6.;
		F = 4., 10., 18., 5., 12., 6., 8., 15., 12.;
		R = 0.1;
		cR = 128.1, 320.1, 576.1, 256.1, 480.1, 384.1, 160.1, 384.1, 192.1;
        MPM::DyadDotUp (v, F, u,  R);
		cout << "e) DyadDotUp error = " << Error(R,cR) << endl;
	}

/** f) Dyadic with Dot (with update): R += (v dyad u)*F  =>  R(i,j)+=v(i)*u(k)*F(k,j) */
	{
		Vector3D v,u;
		ATensor2 F,R,cR;
		u = 1., 2., 3.;
		v = 4., 5., 6.;
		F = 4., 10., 18., 5., 12., 6., 8., 15., 12.;
		R = 0.1;
		cR = 224.1, 350.1, 504.1, 280.1, 420.1, 336.1, 280.1, 420.1, 336.1;
        MPM::DyadDotUp (v, u, F,  R);
		cout << "f) DyadDotUp error = " << Error(R,cR) << endl;
	}

	cout << "\n ~~~~~~~~~~~ Symmetric ~~~~~~~~~~~ " << endl;

//////////// Only assignment

/** 7) Dot product: v = A*u  =>  v(i)=A(i,k)*u(k) */
	{
		STensor2 A;
		Vector3D u;
		Vector3D v,cv;
		A =  1. , 2. , 3. , 4.*SQ2 , 5.*SQ2 , 6.*SQ2;
		u =  1. , 2. , 3.;
		cv = 27., 23., 25.;
        MPM::Dot (A, u,  v);
		cout << "7) Dot       error = " << Error(v,cv) << endl;
	}

/** 8) Dot product: v = u*A  =>  v(i)=u(k)*A(k,i) */
	{
		STensor2 A;
		Vector3D u;
		Vector3D v,cv;
		A =  1. , 2. , 3. , 4.*SQ2 , 5.*SQ2 , 6.*SQ2;
		u =  1. , 2. , 3.;
		cv = 27., 23., 25.;
        MPM::Dot (u, A,  v);
		cout << "8) Dot       error = " << Error(v,cv) << endl;
	}

/** 9) Dot product: y = T*x  =>  y = T:x  =>  y(i,j)=T(i,j,k,l)*x(k,l) */
	{
		STensor4 T;
		STensor2 x,y,cy;
		T =  1.     ,  2.     ,  3.     ,  4.*SQ2 ,  5.*SQ2 ,  6.*SQ2,
		     6.     ,  5.     ,  4.     ,  3.*SQ2 ,  2.*SQ2 ,  1.*SQ2,
		     7.     ,  8.     ,  9.     , 10.*SQ2 , 11.*SQ2 , 12.*SQ2,
		    12.*SQ2 , 11.*SQ2 , 10.*SQ2 ,  9.*2   ,  8.*2   ,  7.*2,
		     1.*SQ2 ,  2.*SQ2 ,  3.*SQ2 ,  4.*2   ,  5.*2   ,  6.*2,
		     6.*SQ2 ,  5.*SQ2 ,  4.*SQ2 ,  3.*2   ,  2.*2   ,  1.*2;
		x =  1. , 2. , 3. , 4.*SQ2 , 5.*SQ2 , 6.*SQ2;
		cy = 168., 84., 384., 300.*SQ2, 168.*SQ2, 84.*SQ2;
        MPM::Dot (T, x,  y);
		cout << "9) Dot       error = " << Error(y,cy) << endl;
	}

//////////////// With update

/** 10) Dot product (with update): v += A*u  =>  v(i)+=A(i,k)*u(k) */
	{
		STensor2 A;
		Vector3D u;
		Vector3D v,cv;
		A =  1. , 2. , 3. , 4.*SQ2 , 5.*SQ2 , 6.*SQ2;
		u =  1. , 2. , 3.;
		v = 0.1;
		cv = 27.1, 23.1, 25.1;
        MPM::DotUp (A, u,  v);
		cout << "10) DotUp    error = " << Error(v,cv) << endl;
	}

/** 11) Dot product (with update): v += u*A  =>  v(i)+=u(k)*A(k,i) */
	{
		STensor2 A;
		Vector3D u;
		Vector3D v,cv;
		A =  1. , 2. , 3. , 4.*SQ2 , 5.*SQ2 , 6.*SQ2;
		u =  1. , 2. , 3.;
		v = 0.1;
		cv = 27.1, 23.1, 25.1;
        MPM::DotUp (u, A,  v);
		cout << "11) DotUp    error = " << Error(v,cv) << endl;
	}

/** 12) Dot product: y += T*x  =>  y += T:x  =>  y(i,j)+=T(i,j,k,l)*x(k,l) */
	{
		STensor4 T;
		STensor2 x,y,cy;
		T =  1.     ,  2.     ,  3.     ,  4.*SQ2 ,  5.*SQ2 ,  6.*SQ2,
		     6.     ,  5.     ,  4.     ,  3.*SQ2 ,  2.*SQ2 ,  1.*SQ2,
		     7.     ,  8.     ,  9.     , 10.*SQ2 , 11.*SQ2 , 12.*SQ2,
		    12.*SQ2 , 11.*SQ2 , 10.*SQ2 ,  9.*2   ,  8.*2   ,  7.*2,
		     1.*SQ2 ,  2.*SQ2 ,  3.*SQ2 ,  4.*2   ,  5.*2   ,  6.*2,
		     6.*SQ2 ,  5.*SQ2 ,  4.*SQ2 ,  3.*2   ,  2.*2   ,  1.*2;
		x =  1. , 2. , 3. , 4.*SQ2 , 5.*SQ2 , 6.*SQ2;
		y = 0.1, 0.1, 0.1, 0.1*SQ2, 0.1*SQ2, 0.1*SQ2;
		cy = 168.1, 84.1, 384.1, 300.1*SQ2, 168.1*SQ2, 84.1*SQ2;
        MPM::DotUp (T, x,  y);
		cout << "12) DotUp    error = " << Error(y,cy) << endl;
	}

	return 0;
}
MECHSYS_CATCH
