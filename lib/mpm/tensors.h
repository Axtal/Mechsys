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

/* Tensors - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_TENSORS_H
#define MPM_TENSORS_H

// STL
#include <cmath>  // for sqrt ...
#include <cfloat> // for DBL_EPSILON
#include <sstream>

// MechSys
#include <mechsys/linalg/matvec.h>

// Local
#include <mechsys/mpm/defs.h>
#include <mechsys/mpm/jacobirot.h>
#include <mechsys/mpm/fmtnum.h>

namespace MPM {

/* Second order tensors:
 *      y             Sy
 *      |             |
 *      |__x      ____|____________
 *    ,'        ,'    |          ,'|
 *   z        ,'     \|/       ,'  |
 *          ,'  <-----'---   ,'    | 
 *        ,'   Sxy         ,'      |
 *      ,'_______________,'   |    |
 *      |                |    |    |
 *      |                |    | <------- Sx
 *      |        __,     |   \|    |
 *      |        ,'|     |    '   ,'
 *      |      ,'        |  Sxy ,'
 *      |    ,'          |   ,'
 *      |  Sz            | ,'
 *      |________________|'
 *
 * (Symmetric) Tensor components (Mandel's basis)
 *   / Sx   Sxy  Sxz \
 *   | Sxy  Sy   Syz | => { Sx, Sy, Sz, SQ2*Sxy, SQ2*Syz, SQ2*Sxz }
 *   \ Sxz  Syz  Sz  /
 *
 * (Assymmetric) Tensor components
 *   / Sx   Sxy  Sxz \
 *   | Syx  Sy   Syz | => { Sx, Sy, Sz, Sxy, Syz, Sxz, Syx, Szy, Szx }
 *   \ Szx  Szy  Sz  /
 */


/////////////////////////////////////////////////////////////////////////////////////////// Entities /////

/** Second order identity tensor (symmetric/Mandel's basis). */
STensor2 SymI;

/** Forth order identity tensor (symmetric/Mandel's basis). */
STensor4 SymII;

/** Forth order tensor (symmetric/Mandel's basis). */
STensor4 SymIdyI;

/** Forth order symmetric-deviatoric tensor (symmetric/Mandel's basis). */
STensor4 SymPsd;

/** Forth order isotropic tensor (symmetric/Mandel's basis). */
STensor4 SymPiso;


///////////////////////////////////////////////////////////////////////////////////////// Operations /////

/** Euclidian norm */

inline double NormV (Vector3D const & v)
{
	return sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2));
}

inline double NormT (STensor2 const & T)
{
	return sqrt(T(0)*T(0) + T(1)*T(1) + T(2)*T(2) + T(3)*T(3) + T(4)*T(4) + T(5)*T(5));
}

/** Add scaled tensors Z = a*X + b*Y. */
inline void AddScaled (double const & a, STensor4 const & X,
                       double const & b, STensor4 const & Y,  STensor4 & Z)
{
	Z(0,0)=a*X(0,0)+b*Y(0,0); Z(0,1)=a*X(0,1)+b*Y(0,1); Z(0,2)=a*X(0,2)+b*Y(0,2); Z(0,3)=a*X(0,3)+b*Y(0,3); Z(0,4)=a*X(0,4)+b*Y(0,4); Z(0,5)=a*X(0,5)+b*Y(0,5);
	Z(1,0)=a*X(1,0)+b*Y(1,0); Z(1,1)=a*X(1,1)+b*Y(1,1); Z(1,2)=a*X(1,2)+b*Y(1,2); Z(1,3)=a*X(1,3)+b*Y(1,3); Z(1,4)=a*X(1,4)+b*Y(1,4); Z(1,5)=a*X(1,5)+b*Y(1,5);
	Z(2,0)=a*X(2,0)+b*Y(2,0); Z(2,1)=a*X(2,1)+b*Y(2,1); Z(2,2)=a*X(2,2)+b*Y(2,2); Z(2,3)=a*X(2,3)+b*Y(2,3); Z(2,4)=a*X(2,4)+b*Y(2,4); Z(2,5)=a*X(2,5)+b*Y(2,5);
	Z(3,0)=a*X(3,0)+b*Y(3,0); Z(3,1)=a*X(3,1)+b*Y(3,1); Z(3,2)=a*X(3,2)+b*Y(3,2); Z(3,3)=a*X(3,3)+b*Y(3,3); Z(3,4)=a*X(3,4)+b*Y(3,4); Z(3,5)=a*X(3,5)+b*Y(3,5);
	Z(4,0)=a*X(4,0)+b*Y(4,0); Z(4,1)=a*X(4,1)+b*Y(4,1); Z(4,2)=a*X(4,2)+b*Y(4,2); Z(4,3)=a*X(4,3)+b*Y(4,3); Z(4,4)=a*X(4,4)+b*Y(4,4); Z(4,5)=a*X(4,5)+b*Y(4,5);
	Z(5,0)=a*X(5,0)+b*Y(5,0); Z(5,1)=a*X(5,1)+b*Y(5,1); Z(5,2)=a*X(5,2)+b*Y(5,2); Z(5,3)=a*X(5,3)+b*Y(5,3); Z(5,4)=a*X(5,4)+b*Y(5,4); Z(5,5)=a*X(5,5)+b*Y(5,5);
}


////////////////////////////////////////////////////////////////////////////////////////// Functions /////

/** Symmetric of a tensor S = a*sym(T) = a * 0.5*(T + trn(T)). */
inline void Sym (double a, ATensor2 const & T,  STensor2 & S)
{
	S(0) = a*T(0);
	S(1) = a*T(1);
	S(2) = a*T(2);
	S(3) = a*0.5*(T(3)+T(6))*SQ2;
	S(4) = a*0.5*(T(4)+T(7))*SQ2;
	S(5) = a*0.5*(T(5)+T(8))*SQ2;
}

/** Determinant. */
inline double Det (ATensor2 const & B)
{
	return B(0)*(B(1)*B(2)-B(4)*B(7))-B(3)*(B(6)*B(2)-B(4)*B(8))+B(5)*(B(6)*B(7)-B(1)*B(8));
}

/** Determinant. */
inline double Det (STensor2 const & T)
{
	return T(0)*T(1)*T(2) + T(3)*T(4)*T(5)/SQ2 - T(0)*T(4)*T(4)/2.0 - T(1)*T(5)*T(5)/2.0 - T(2)*T(3)*T(3)/2.0;
}

/** Cambridge's p mean stress invariant. */
inline double Cam_p (STensor2 const & S)
{
	return (S(0)+S(1)+S(2))/3.0;
}

/** Cambridge's q deviatoric stress invariant. */
inline double Cam_q (STensor2 const & S)
{
	return sqrt((   (S(0)-S(1))*(S(0)-S(1))
	              + (S(1)-S(2))*(S(1)-S(2))
	              + (S(2)-S(0))*(S(2)-S(0)) + 3.0*(S(3)*S(3)+S(4)*S(4)+S(5)*S(5)) )/2.0);
}

inline void Stress_p_q (STensor2 const & Sig, double & p, double & q)
{
	p = (Sig(0)+Sig(1)+Sig(2))/3.0;
	q = sqrt((   (Sig(0)-Sig(1))*(Sig(0)-Sig(1))
			   + (Sig(1)-Sig(2))*(Sig(1)-Sig(2))
			   + (Sig(2)-Sig(0))*(Sig(2)-Sig(0)) + 3.0*(Sig(3)*Sig(3) + Sig(4)*Sig(4) + Sig(5)*Sig(5)) )/2.0);
}

inline void Eigenvals (STensor2 const & T, double L[3])
{
	// Calculate eigenvalues and eigenvectors
	if (JacobiRot(T,L)==-1)
	   throw new Fatal("Eigenvals: Jacobi rotation did not converge. T=<%.6f %.6f %.6f %.6f %.6f %.6f>",T(0),T(1),T(2),T(3)/SQ2,T(4)/SQ2,T(5)/SQ2);
}

inline void Eigenvp (STensor2 const & T, double L[3], STensor2 P[3])
{
	double V0[3]; // eigenvector
	double V1[3]; // eigenvector
	double V2[3]; // eigenvector

	// Calculate eigenvalues and eigenvectors
	if (JacobiRot(T,V0,V1,V2,L)==-1)
	   throw new Fatal("Eigenvp: Jacobi rotation did not converge. T=<%.6f %.6f %.6f %.6f %.6f %.6f>",T(0),T(1),T(2),T(3)/SQ2,T(4)/SQ2,T(5)/SQ2);

	// Calculate eigenprojectors (Pi = Vi dyad Vi)
	// Eigenprojector 0
	P[0](0) = V0[0] * V0[0];
	P[0](1) = V0[1] * V0[1];
	P[0](2) = V0[2] * V0[2];
	P[0](3) = V0[0] * V0[1] * SQ2;
	P[0](4) = V0[1] * V0[2] * SQ2;
	P[0](5) = V0[0] * V0[2] * SQ2;
	// Eigenprojector 1
	P[1](0) = V1[0] * V1[0];
	P[1](1) = V1[1] * V1[1];
	P[1](2) = V1[2] * V1[2];
	P[1](3) = V1[0] * V1[1] * SQ2;
	P[1](4) = V1[1] * V1[2] * SQ2;
	P[1](5) = V1[0] * V1[2] * SQ2;
	// Eigenprojector 2
	P[2](0) = V2[0] * V2[0];
	P[2](1) = V2[1] * V2[1];
	P[2](2) = V2[2] * V2[2];
	P[2](3) = V2[0] * V2[1] * SQ2;
	P[2](4) = V2[1] * V2[2] * SQ2;
	P[2](5) = V2[0] * V2[2] * SQ2;
}

inline void Stress_psmp_qsmp (STensor2 const & Sig, double & psmp, double & qsmp)
{
	// Principal values
	double L[3];
	Eigenvals (Sig, L);

	// Invariants
	if (L[0]>0.0 && L[1]>0.0 && L[2]>0.0)
	{
		// Characteristic invariants
		double Is[4];
		Is[0] = L[0]+L[1]+L[2];
		Is[1] = L[0]*L[1]+L[1]*L[2]+L[2]*L[0];
		Is[2] = L[0]*L[1]*L[2];
		Is[3] = sqrt(L[0])+sqrt(L[1])+sqrt(L[2]);

		// psmp and qsmp invariants
		if (Is[1]>0.0)
		{
			psmp = Is[3]*sqrt(Is[2]/Is[1]);
			double del = L[0]*L[0]+L[1]*L[1]+L[2]*L[2]-psmp*psmp;
			if (fabs(del)<DBL_EPSILON*100.0) del = 0.0;
			if (del>=0.0) qsmp = sqrt(del);
			else qsmp=0.0; //throw new Fatal(_("Stress_psmp_qsmp: del=s1^2+s2^2+s3^2-psmp^2=%f must be positive. s1=%f, s2=%f, s3=%f, psmp=%f"),del,L[0],L[1],L[2],psmp);
		}
		else throw new Fatal("Stress_psmp_qsmp: I2=%f must be positive",Is[1]);
	}
	else throw new Fatal("Stress_psmp_qsmp: The principal stresses must be positive. s1=%f, s2=%f, s3=%f",L[0],L[1],L[2]);
}

inline void Derivs_psmp_qsmp (STensor2 const & Sig, double & psmp, double & qsmp, double DpsmpDs[3], double DqsmpDs[3], STensor2 Proj[3])
{
	// Eigenvalues and Eigenprojectors
	double L[3];
	Eigenvp (Sig, L, Proj);

	// Invariants and derivatives
	if (L[0]>0.0 && L[1]>0.0 && L[2]>0.0)
	{
		// Characteristic invariants
		double Is[4];
		Is[0] = L[0]+L[1]+L[2];
		Is[1] = L[0]*L[1]+L[1]*L[2]+L[2]*L[0];
		Is[2] = L[0]*L[1]*L[2];
		Is[3] = sqrt(L[0])+sqrt(L[1])+sqrt(L[2]);

		// Derivatives and invariants
		if (Is[1]>0.0)
		{
			// Invariants
			psmp = Is[3]*sqrt(Is[2]/Is[1]);
			double del = L[0]*L[0]+L[1]*L[1]+L[2]*L[2]-psmp*psmp;
			if (fabs(del)<DBL_EPSILON*100.0) del = 0.0;
			if (del>=0.0) qsmp = sqrt(del);
			else throw new Fatal("Derivs_psmp_qsmp: del=s1^2+s2^2+s3^2-psmp^2=%f must be positive. s1=%f, s2=%f, s3=%f, psmp=%f",del,L[0],L[1],L[2],psmp);

			// Derivatives
			double dI1ds[3], dI2ds[3], dI3ds[3];
			for (int k=0; k<3; ++k)
			{
				// Derivatives of the characteristic invariants
				dI1ds[k] = Is[0]-L[k];
				dI2ds[k] = Is[2]/L[k];
				dI3ds[k] = 0.5/sqrt(L[k]);

				// Derivatives of the psmp and qsmp invariants
				DpsmpDs[k] = (-0.5*Is[3]*sqrt(Is[2])/pow(Is[1],1.5))*dI1ds[k] + (0.5*Is[3]/sqrt(Is[1]*Is[2]))*dI2ds[k] + sqrt(Is[2]/Is[1])*dI3ds[k];
				if (qsmp>DBL_EPSILON*100.0) DqsmpDs[k] = (1.0/qsmp)*(L[k]-psmp*DpsmpDs[k]);
				else                        DqsmpDs[k] = 0.0;
			}
		}
		else throw new Fatal("Derivs_psmp_qsmp: I2=%f must be positive",Is[1]);
	}
	else throw new Fatal("Derivs_psmp_qsmp: The principal stresses must be positive. s1=%f, s2=%f, s3=%f",L[0],L[1],L[2]);
}

inline void Stress_ssmp_tsmp (STensor2 const & Sig, double & ssmp, double & tsmp)
{
	// Principal values
	double L[3];
	Eigenvals (Sig, L);

	// Invariants
	if (L[0]>0.0 && L[1]>0.0 && L[2]>0.0)
	{
		// Characteristic invariants
		double Is[4];
		Is[0] = L[0]+L[1]+L[2];
		Is[1] = L[0]*L[1]+L[1]*L[2]+L[2]*L[0];
		Is[2] = L[0]*L[1]*L[2];
		Is[3] = sqrt(L[0])+sqrt(L[1])+sqrt(L[2]);

		// ssmp and tsmp invariants
		if (Is[1]>0.0)
		{
			ssmp = 3.0*Is[2]/Is[1];
			double del = Is[0]*Is[2]/Is[1]-ssmp*ssmp;
			if (fabs(del)<DBL_EPSILON*100.0) del = 0.0;
			if (del>=0.0) tsmp = sqrt(del);
			else throw new Fatal("Stress_ssmp_tsmp: del=I1*I3/I2-ssmp*ssmp=%f must be positive. I1=%f, I2=%f, I3=%f, ssmp=%f",del,Is[0],Is[1],Is[2],ssmp);
		}
		else throw new Fatal("Stress_ssmp_tsmp: I2=%f must be positive",Is[1]);
	}
	else throw new Fatal("Stress_ssmp_tsmp: The principal stresses must be positive. s1=%f, s2=%f, s3=%f",L[0],L[1],L[2]);
}

inline void Derivs_ssmp_tsmp (STensor2 const & Sig, double & ssmp, double & tsmp, double DssmpDs[3], double DtsmpDs[3], STensor2 Proj[3])
{
	// Eigenvalues and Eigenprojectors
	double L[3];
	Eigenvp (Sig, L, Proj);

	// Invariants and derivatives
	if (L[0]>0.0 && L[1]>0.0 && L[2]>0.0)
	{
		// Characteristic invariants
		double Is[4];
		Is[0] = L[0]+L[1]+L[2];
		Is[1] = L[0]*L[1]+L[1]*L[2]+L[2]*L[0];
		Is[2] = L[0]*L[1]*L[2];
		Is[3] = sqrt(L[0])+sqrt(L[1])+sqrt(L[2]);

		// Derivatives and invariants
		if (Is[1]>0.0)
		{
			// Invariants
			ssmp = 3.0*Is[2]/Is[1];
			double del = Is[0]*Is[2]/Is[1]-ssmp*ssmp;
			if (fabs(del)<DBL_EPSILON*100.0) del = 0.0;
			if (del>=0.0) tsmp = sqrt(del);
			else throw new Fatal("Derivs_ssmp_tsmp: del=I1*I3/I2-ssmp*ssmp=%f must be positive. I1=%f, I2=%f, I3=%f, ssmp=%f",del,Is[0],Is[1],Is[2],ssmp);

			// Derivatives
			double dI0ds[3], dI1ds[3], dI2ds[3];
			for (int k=0; k<3; ++k)
			{
				// Derivatives of the characteristic invariants
				dI0ds[k] = 1.0;
				dI1ds[k] = Is[0]-L[k];
				dI2ds[k] = Is[2]/L[k];

				// Derivatives of the ssmp and tsmp invariants
				DssmpDs[k] = (-3.0*Is[2]/(Is[1]*Is[1]))*dI1ds[k] + (3.0/Is[1])*dI2ds[k];
				if (tsmp>DBL_EPSILON*100.0) DtsmpDs[k] = (0.5/tsmp)*((Is[2]/Is[1])*dI0ds[k]-(Is[0]*Is[2]/(Is[1]*Is[1]))*dI1ds[k]+(Is[0]/Is[1])*dI2ds[k]-2.0*ssmp*DssmpDs[k]);
				else                        DtsmpDs[k] = 0.0;
			}
		}
		else throw new Fatal("Derivs_ssmp_tsmp: I2=%f must be positive",Is[1]);
	}
	else throw new Fatal("Derivs_ssmp_tsmp: The principal stresses must be positive. s1=%f, s2=%f, s3=%f",L[0],L[1],L[2]);
}

inline void Derivs_p_q (STensor2 const & Sig, double & p, double & q, double DpDs[3], double DqDs[3], STensor2 Proj[3])
{
	// Eigenvalues and Eigenprojectors
	double L[3];
	Eigenvp (Sig, L, Proj);

	// Invariants
	p = (L[0]+L[1]+L[2])/3.0;
	q = sqrt((L[0]-L[1])*(L[0]-L[1])+(L[1]-L[2])*(L[1]-L[2])+(L[2]-L[0])*(L[2]-L[0]))/SQ2;

	// Derivatives
	for (int k=0; k<3; ++k)
	{
		DpDs[k] = 1.0/3.0;
		if (q>DBL_EPSILON*100) DqDs[k] = 3.0*(L[k]-p)/(2.0*q);
		else                   DqDs[k] = 0.0;
	}
}

inline void Stress_p_q_t (STensor2 const & Sig, double & p, double & q, double & t)
{
	Stress_p_q (Sig,p,q);
	STensor2 S;  S = Sig - SymI * p;
	if (q<=DBL_EPSILON*100.0) t = -1.0; // 3th==-90 => th=30
	else                      t = 27.0*Det(S)/(2.0*q*q*q);
	if (t>= 1.0)              t =  1.0;
	if (t<=-1.0)              t = -1.0;
}

inline void Derivs_p_q_t (STensor2 const & Sig, double & p, double & q, double & t, double DpDs[3], double DqDs[3], double DtDs[3], STensor2 Proj[3])
{
	// Eigenvalues and Eigenprojectors
	double L[3];
	Eigenvp (Sig, L, Proj);

	// Invariants
	p = (L[0]+L[1]+L[2])/3.0;
	q = sqrt((L[0]-L[1])*(L[0]-L[1])+(L[1]-L[2])*(L[1]-L[2])+(L[2]-L[0])*(L[2]-L[0]))/SQ2;
	double S[3]; for (int k=0; k<3; ++k) S[k] = L[k]-p;
	if (q>DBL_EPSILON*100.0) t = 27.0*S[0]*S[1]*S[2]/(2.0*q*q*q);
	else                     t = 0.0;

	// Derivatives
	double B[3];
	B[0] = L[2]-L[1];
	B[1] = L[0]-L[2];
	B[2] = L[1]-L[0];
	double l = (L[0]-L[1])*(L[1]-L[2])*(L[2]-L[0]);
	for (int k=0; k<3; ++k)
	{
		DpDs[k] = 1.0/3.0;
		if (q>DBL_EPSILON*100)
		{
			DqDs[k] = 3.0*S[k]/(2.0*q);
			DtDs[k] = 27.0*l*B[k]/(4.0*pow(q,5.0));
		}
		else
		{
			DqDs[k] = 0.0;
			DtDs[k] = 0.0;
		}
	}
}

inline void Strain_Ev_Ed (STensor2 const & Eps, double & Ev, double & Ed)
{
	Ev = Eps(0)+Eps(1)+Eps(2);
	Ed = sqrt(2.0*(   (Eps(0)-Eps(1))*(Eps(0)-Eps(1))
				    + (Eps(1)-Eps(2))*(Eps(1)-Eps(2))
					+ (Eps(2)-Eps(0))*(Eps(2)-Eps(0)) + 3.0*(Eps(3)*Eps(3) + Eps(4)*Eps(4) + Eps(5)*Eps(5)) ))/3.0;
}

/** Output stress and strain values and invariants for use with plotn in plot.R - HEADER. */
inline void OutputHeader (std::ostream & Os, bool NewLine=true)
{
	Os << _8s<< "Sx" << _8s<< "Sy" << _8s<< "Sz" << _8s<< "Sxy" << _8s<< "Syz" << _8s<< "Sxz";
	Os << _8s<< "Ex" << _8s<< "Ey" << _8s<< "Ez" << _8s<< "Exy" << _8s<< "Eyz" << _8s<< "Exz";
	Os << _8s<< "p"  << _8s<< "q"  << _8s<< "Ev" << _8s<< "Ed";
	if (NewLine) Os << std::endl;
}

/** Output stress and strain values and invariants for use with plotn in plot.R. */
inline void Output (STensor2 const & Sig, STensor2 const & Eps,  std::ostream & Os, bool NewLine=true)
{
	double p,q,ev,ed;
	Stress_p_q   (Sig,p,q);
	Strain_Ev_Ed (Eps,ev,ed);
	Os << _8s<< Sig(0)       << _8s<< Sig(1)       << _8s<< Sig(2)        << _8s<< Sig(3)/SQ2       << _8s<< Sig(4)/SQ2       << _8s<< Sig(5)/SQ2;
	Os << _8s<< Eps(0)*100.0 << _8s<< Eps(1)*100.0 << _8s<< Eps(2)*100.0  << _8s<< Eps(3)*100.0/SQ2 << _8s<< Eps(4)*100.0/SQ2 << _8s<< Eps(5)*100.0/SQ2;
	Os << _8s<< p            << _8s<< q            << _8s<< ev*100.0      << _8s<< ed*100.0;
	if (NewLine) Os << std::endl;
}

/** Converts hydrostatic coordinates to sigma (principal) coord (A in radians). */
inline void Hid2Sig (Vector3D const & RAZ, Vector3D & XYZ)
{
	double SA = -RAZ(0)*sin(RAZ(1));
	double SB =  RAZ(0)*cos(RAZ(1));
	XYZ = RAZ(2) + 2.0*SB/SQ6         ,
	      RAZ(2) -     SB/SQ6 - SA/SQ2,
	      RAZ(2) -     SB/SQ6 + SA/SQ2;
}


/////////////////////////////////////////////////////////////////////// 3D

///////////////////////////////////////// Asymmetric

//////////// Only assignment

/** 1) Dot product: v=A*u  =>  v(i)=A(i,k)*u(k) */
inline void Dot (ATensor2 const & A, Vector3D const & u,  Vector3D & v)
{
	v(0) = A(0)*u(0)+A(3)*u(1)+A(5)*u(2);
	v(1) = A(6)*u(0)+A(1)*u(1)+A(4)*u(2);
	v(2) = A(8)*u(0)+A(7)*u(1)+A(2)*u(2);
}
 
/** 2) Dot product: v=u*A  =>  v(i)=u(k)*A(k,i) */
inline void Dot (Vector3D const & u, ATensor2 const & A,  Vector3D & v)
{
	v(0) = u(0)*A(0)+u(1)*A(6)+u(2)*A(8);
	v(1) = u(0)*A(3)+u(1)*A(1)+u(2)*A(7);
	v(2) = u(0)*A(5)+u(1)*A(4)+u(2)*A(2);
}

/** 3) Dot product: y = T*x  =>  y = T:x  =>  y(i,j)=T(i,j,k,l)*x(k,l) */
inline void Dot (ATensor4 const & T, ATensor2 const & x,  ATensor2 & y)
{
	y(0) = T(0,0)*x(0)+T(0,3)*x(3)+T(0,5)*x(5)+T(0,6)*x(6)+T(0,1)*x(1)+T(0,4)*x(4)+T(0,8)*x(8)+T(0,7)*x(7)+T(0,2)*x(2);
	y(3) = T(3,0)*x(0)+T(3,3)*x(3)+T(3,5)*x(5)+T(3,6)*x(6)+T(3,1)*x(1)+T(3,4)*x(4)+T(3,8)*x(8)+T(3,7)*x(7)+T(3,2)*x(2);
	y(5) = T(5,0)*x(0)+T(5,3)*x(3)+T(5,5)*x(5)+T(5,6)*x(6)+T(5,1)*x(1)+T(5,4)*x(4)+T(5,8)*x(8)+T(5,7)*x(7)+T(5,2)*x(2);
	y(6) = T(6,0)*x(0)+T(6,3)*x(3)+T(6,5)*x(5)+T(6,6)*x(6)+T(6,1)*x(1)+T(6,4)*x(4)+T(6,8)*x(8)+T(6,7)*x(7)+T(6,2)*x(2);
	y(1) = T(1,0)*x(0)+T(1,3)*x(3)+T(1,5)*x(5)+T(1,6)*x(6)+T(1,1)*x(1)+T(1,4)*x(4)+T(1,8)*x(8)+T(1,7)*x(7)+T(1,2)*x(2);
	y(4) = T(4,0)*x(0)+T(4,3)*x(3)+T(4,5)*x(5)+T(4,6)*x(6)+T(4,1)*x(1)+T(4,4)*x(4)+T(4,8)*x(8)+T(4,7)*x(7)+T(4,2)*x(2);
	y(8) = T(8,0)*x(0)+T(8,3)*x(3)+T(8,5)*x(5)+T(8,6)*x(6)+T(8,1)*x(1)+T(8,4)*x(4)+T(8,8)*x(8)+T(8,7)*x(7)+T(8,2)*x(2);
	y(7) = T(7,0)*x(0)+T(7,3)*x(3)+T(7,5)*x(5)+T(7,6)*x(6)+T(7,1)*x(1)+T(7,4)*x(4)+T(7,8)*x(8)+T(7,7)*x(7)+T(7,2)*x(2);
	y(2) = T(2,0)*x(0)+T(2,3)*x(3)+T(2,5)*x(5)+T(2,6)*x(6)+T(2,1)*x(1)+T(2,4)*x(4)+T(2,8)*x(8)+T(2,7)*x(7)+T(2,2)*x(2);
}

/** 4) Dot product: C=A*B  =>  C(i,j)=A(i,k)*B(k,j) */
inline void Dot (ATensor2 const & A, ATensor2 const & B,  ATensor2 & C)
{
	C(0) = A(0)*B(0)+A(3)*B(6)+A(5)*B(8);
	C(3) = A(0)*B(3)+A(3)*B(1)+A(5)*B(7);
	C(5) = A(0)*B(5)+A(3)*B(4)+A(5)*B(2);
	C(6) = A(6)*B(0)+A(1)*B(6)+A(4)*B(8);
	C(1) = A(6)*B(3)+A(1)*B(1)+A(4)*B(7);
	C(4) = A(6)*B(5)+A(1)*B(4)+A(4)*B(2);
	C(8) = A(8)*B(0)+A(7)*B(6)+A(2)*B(8);
	C(7) = A(8)*B(3)+A(7)*B(1)+A(2)*B(7);
	C(2) = A(8)*B(5)+A(7)*B(4)+A(2)*B(2);
}

/** a) Dyadic product: A = u dyad v  =>  A(i,j)=u(i)*v(j) */
inline void Dyad (Vector3D const & u, Vector3D const & v,  ATensor2 & A)
{
	A(0) = u(0)*v(0);
	A(3) = u(0)*v(1);
	A(5) = u(0)*v(2);
	A(6) = u(1)*v(0);
	A(1) = u(1)*v(1);
	A(4) = u(1)*v(2);
	A(8) = u(2)*v(0);
	A(7) = u(2)*v(1);
	A(2) = u(2)*v(2);
}

/** b) Dyadic with Dot: R=v dyad (F*u)  =>  R(i,j)=v(i)*F(j,k)*u(k) */
inline void DyadDot (Vector3D const & v, ATensor2 const & F, Vector3D const & u,  ATensor2 & R)
{
	R(0) = v(0)*F(0)*u(0)+v(0)*F(3)*u(1)+v(0)*F(5)*u(2);
	R(3) = v(0)*F(6)*u(0)+v(0)*F(1)*u(1)+v(0)*F(4)*u(2);
	R(5) = v(0)*F(8)*u(0)+v(0)*F(7)*u(1)+v(0)*F(2)*u(2);
	R(6) = v(1)*F(0)*u(0)+v(1)*F(3)*u(1)+v(1)*F(5)*u(2);
	R(1) = v(1)*F(6)*u(0)+v(1)*F(1)*u(1)+v(1)*F(4)*u(2);
	R(4) = v(1)*F(8)*u(0)+v(1)*F(7)*u(1)+v(1)*F(2)*u(2);
	R(8) = v(2)*F(0)*u(0)+v(2)*F(3)*u(1)+v(2)*F(5)*u(2);
	R(7) = v(2)*F(6)*u(0)+v(2)*F(1)*u(1)+v(2)*F(4)*u(2);
	R(2) = v(2)*F(8)*u(0)+v(2)*F(7)*u(1)+v(2)*F(2)*u(2);
}

/** c) Dyadic with Dot: R=(v dyad u)*F  =>  R(i,j)=v(i)*u(k)*F(k,j) */
inline void DyadDot (Vector3D const & v, Vector3D const & u, ATensor2 const & F,  ATensor2 & R)
{
	R(0) = v(0)*u(0)*F(0)+v(0)*u(1)*F(6)+v(0)*u(2)*F(8);
	R(3) = v(0)*u(0)*F(3)+v(0)*u(1)*F(1)+v(0)*u(2)*F(7);
	R(5) = v(0)*u(0)*F(5)+v(0)*u(1)*F(4)+v(0)*u(2)*F(2);
	R(6) = v(1)*u(0)*F(0)+v(1)*u(1)*F(6)+v(1)*u(2)*F(8);
	R(1) = v(1)*u(0)*F(3)+v(1)*u(1)*F(1)+v(1)*u(2)*F(7);
	R(4) = v(1)*u(0)*F(5)+v(1)*u(1)*F(4)+v(1)*u(2)*F(2);
	R(8) = v(2)*u(0)*F(0)+v(2)*u(1)*F(6)+v(2)*u(2)*F(8);
	R(7) = v(2)*u(0)*F(3)+v(2)*u(1)*F(1)+v(2)*u(2)*F(7);
	R(2) = v(2)*u(0)*F(5)+v(2)*u(1)*F(4)+v(2)*u(2)*F(2);
}

//////////////// With update

/** 4) Dot product (with update): v += A*u  =>  v(i)+=A(i,k)*u(k) */
inline void DotUp (ATensor2 const & A, Vector3D const & u,  Vector3D & v)
{
	v(0) += A(0)*u(0)+A(3)*u(1)+A(5)*u(2);
	v(1) += A(6)*u(0)+A(1)*u(1)+A(4)*u(2);
	v(2) += A(8)*u(0)+A(7)*u(1)+A(2)*u(2);
}
 
/** 5) Dot product (with update): v += u*A  =>  v(i)+=u(k)*A(k,i) */
inline void DotUp (Vector3D const & u, ATensor2 const & A,  Vector3D & v)
{
	v(0) += u(0)*A(0)+u(1)*A(6)+u(2)*A(8);
	v(1) += u(0)*A(3)+u(1)*A(1)+u(2)*A(7);
	v(2) += u(0)*A(5)+u(1)*A(4)+u(2)*A(2);
}

/** 6) Dot product (with update): y += T*x  =>  y = T:x  =>  y(i,j)+=T(i,j,k,l)*x(k,l) */
inline void DotUp (ATensor4 const & T, ATensor2 const & x,  ATensor2 & y)
{
	y(0) += T(0,0)*x(0)+T(0,3)*x(3)+T(0,5)*x(5)+T(0,6)*x(6)+T(0,1)*x(1)+T(0,4)*x(4)+T(0,8)*x(8)+T(0,7)*x(7)+T(0,2)*x(2);
	y(3) += T(3,0)*x(0)+T(3,3)*x(3)+T(3,5)*x(5)+T(3,6)*x(6)+T(3,1)*x(1)+T(3,4)*x(4)+T(3,8)*x(8)+T(3,7)*x(7)+T(3,2)*x(2);
	y(5) += T(5,0)*x(0)+T(5,3)*x(3)+T(5,5)*x(5)+T(5,6)*x(6)+T(5,1)*x(1)+T(5,4)*x(4)+T(5,8)*x(8)+T(5,7)*x(7)+T(5,2)*x(2);
	y(6) += T(6,0)*x(0)+T(6,3)*x(3)+T(6,5)*x(5)+T(6,6)*x(6)+T(6,1)*x(1)+T(6,4)*x(4)+T(6,8)*x(8)+T(6,7)*x(7)+T(6,2)*x(2);
	y(1) += T(1,0)*x(0)+T(1,3)*x(3)+T(1,5)*x(5)+T(1,6)*x(6)+T(1,1)*x(1)+T(1,4)*x(4)+T(1,8)*x(8)+T(1,7)*x(7)+T(1,2)*x(2);
	y(4) += T(4,0)*x(0)+T(4,3)*x(3)+T(4,5)*x(5)+T(4,6)*x(6)+T(4,1)*x(1)+T(4,4)*x(4)+T(4,8)*x(8)+T(4,7)*x(7)+T(4,2)*x(2);
	y(8) += T(8,0)*x(0)+T(8,3)*x(3)+T(8,5)*x(5)+T(8,6)*x(6)+T(8,1)*x(1)+T(8,4)*x(4)+T(8,8)*x(8)+T(8,7)*x(7)+T(8,2)*x(2);
	y(7) += T(7,0)*x(0)+T(7,3)*x(3)+T(7,5)*x(5)+T(7,6)*x(6)+T(7,1)*x(1)+T(7,4)*x(4)+T(7,8)*x(8)+T(7,7)*x(7)+T(7,2)*x(2);
	y(2) += T(2,0)*x(0)+T(2,3)*x(3)+T(2,5)*x(5)+T(2,6)*x(6)+T(2,1)*x(1)+T(2,4)*x(4)+T(2,8)*x(8)+T(2,7)*x(7)+T(2,2)*x(2);
}

/** d) Dyadic product (with update): A += u dyad v  =>  A(i,j)+=u(i)*v(j) */
inline void DyadUp (Vector3D const & u, Vector3D const & v,  ATensor2 & A)
{
	A(0) += u(0)*v(0);
	A(3) += u(0)*v(1);
	A(5) += u(0)*v(2);
	A(6) += u(1)*v(0);
	A(1) += u(1)*v(1);
	A(4) += u(1)*v(2);
	A(8) += u(2)*v(0);
	A(7) += u(2)*v(1);
	A(2) += u(2)*v(2);
}

/** da) Scaled Dyadic product (with update): A += s *(u dyad v)  =>  A(i,j)+=s*u(i)*v(j) */
inline void ScDyadUp (double s, Vector3D const & u, Vector3D const & v,  ATensor2 & A)
{
	A(0) += s*u(0)*v(0);
	A(3) += s*u(0)*v(1);
	A(5) += s*u(0)*v(2);
	A(6) += s*u(1)*v(0);
	A(1) += s*u(1)*v(1);
	A(4) += s*u(1)*v(2);
	A(8) += s*u(2)*v(0);
	A(7) += s*u(2)*v(1);
	A(2) += s*u(2)*v(2);
}

/** e) Dyadic with Dot (with update): R += v dyad (F*u)  =>  R(i,j)+=v(i)*F(j,k)*u(k) */
inline void DyadDotUp (Vector3D const & v, ATensor2 const & F, Vector3D const & u,  ATensor2 & R)
{
	R(0) += v(0)*F(0)*u(0)+v(0)*F(3)*u(1)+v(0)*F(5)*u(2);
	R(3) += v(0)*F(6)*u(0)+v(0)*F(1)*u(1)+v(0)*F(4)*u(2);
	R(5) += v(0)*F(8)*u(0)+v(0)*F(7)*u(1)+v(0)*F(2)*u(2);
	R(6) += v(1)*F(0)*u(0)+v(1)*F(3)*u(1)+v(1)*F(5)*u(2);
	R(1) += v(1)*F(6)*u(0)+v(1)*F(1)*u(1)+v(1)*F(4)*u(2);
	R(4) += v(1)*F(8)*u(0)+v(1)*F(7)*u(1)+v(1)*F(2)*u(2);
	R(8) += v(2)*F(0)*u(0)+v(2)*F(3)*u(1)+v(2)*F(5)*u(2);
	R(7) += v(2)*F(6)*u(0)+v(2)*F(1)*u(1)+v(2)*F(4)*u(2);
	R(2) += v(2)*F(8)*u(0)+v(2)*F(7)*u(1)+v(2)*F(2)*u(2);
}

/** f) Dyadic with Dot (with update): R += (v dyad u)*F  =>  R(i,j)+=v(i)*u(k)*F(k,j) */
inline void DyadDotUp (Vector3D const & v, Vector3D const & u, ATensor2 const & F,  ATensor2 & R)
{
	R(0) += v(0)*u(0)*F(0)+v(0)*u(1)*F(6)+v(0)*u(2)*F(8);
	R(3) += v(0)*u(0)*F(3)+v(0)*u(1)*F(1)+v(0)*u(2)*F(7);
	R(5) += v(0)*u(0)*F(5)+v(0)*u(1)*F(4)+v(0)*u(2)*F(2);
	R(6) += v(1)*u(0)*F(0)+v(1)*u(1)*F(6)+v(1)*u(2)*F(8);
	R(1) += v(1)*u(0)*F(3)+v(1)*u(1)*F(1)+v(1)*u(2)*F(7);
	R(4) += v(1)*u(0)*F(5)+v(1)*u(1)*F(4)+v(1)*u(2)*F(2);
	R(8) += v(2)*u(0)*F(0)+v(2)*u(1)*F(6)+v(2)*u(2)*F(8);
	R(7) += v(2)*u(0)*F(3)+v(2)*u(1)*F(1)+v(2)*u(2)*F(7);
	R(2) += v(2)*u(0)*F(5)+v(2)*u(1)*F(4)+v(2)*u(2)*F(2);
}

////////////////////////////////////////// Symmetric

//////////// Only assignment

/** 7) Dot product: v = A*u  =>  v(i)=A(i,k)*u(k) */
inline void Dot (STensor2 const & A, Vector3D const & u,  Vector3D & v)
{
	v(0) = A(0)*u(0)     + A(3)*u(1)/SQ2 + A(5)*u(2)/SQ2;
	v(1) = A(3)*u(0)/SQ2 + A(1)*u(1)     + A(4)*u(2)/SQ2;
	v(2) = A(5)*u(0)/SQ2 + A(4)*u(1)/SQ2 + A(2)*u(2);
}

/** 8) Dot product: v = u*A  =>  v(i)=u(k)*A(k,i) */
inline void Dot (Vector3D const & u, STensor2 const & A,  Vector3D & v)
{
	v(0) = u(0)*A(0)     + u(1)*A(3)/SQ2 + u(2)*A(5)/SQ2;
	v(1) = u(0)*A(3)/SQ2 + u(1)*A(1)     + u(2)*A(4)/SQ2;
	v(2) = u(0)*A(5)/SQ2 + u(1)*A(4)/SQ2 + u(2)*A(2);
}

/** 9) Dot product: y = T*x  =>  y = T:x  =>  y(i,j)=T(i,j,k,l)*x(k,l) */
inline void Dot (STensor4 const & T, STensor2 const & x,  STensor2 & y)
{
	y(0) = T(0,0)*x(0)+T(0,1)*x(1)+T(0,2)*x(2)+T(0,3)*x(3)+T(0,4)*x(4)+T(0,5)*x(5);
	y(1) = T(1,0)*x(0)+T(1,1)*x(1)+T(1,2)*x(2)+T(1,3)*x(3)+T(1,4)*x(4)+T(1,5)*x(5);
	y(2) = T(2,0)*x(0)+T(2,1)*x(1)+T(2,2)*x(2)+T(2,3)*x(3)+T(2,4)*x(4)+T(2,5)*x(5);
	y(3) = T(3,0)*x(0)+T(3,1)*x(1)+T(3,2)*x(2)+T(3,3)*x(3)+T(3,4)*x(4)+T(3,5)*x(5);
	y(4) = T(4,0)*x(0)+T(4,1)*x(1)+T(4,2)*x(2)+T(4,3)*x(3)+T(4,4)*x(4)+T(4,5)*x(5);
	y(5) = T(5,0)*x(0)+T(5,1)*x(1)+T(5,2)*x(2)+T(5,3)*x(3)+T(5,4)*x(4)+T(5,5)*x(5);
}

/** 9x) Scaled Dot product: y = s*x*T  =>  y = s*x:T  =>  y(i,j)=s*x(k,l)*T(k,l,i,j) */
inline void ScDot (double s, STensor2 const & x, STensor4 const & T,  STensor2 & y)
{
	y(0) = s*(x(0)*T(0,0)+x(1)*T(1,0)+x(2)*T(2,0)+x(3)*T(3,0)+x(4)*T(4,0)+x(5)*T(5,0));
	y(1) = s*(x(0)*T(0,1)+x(1)*T(1,1)+x(2)*T(2,1)+x(3)*T(3,1)+x(4)*T(4,1)+x(5)*T(5,1));
	y(2) = s*(x(0)*T(0,2)+x(1)*T(1,2)+x(2)*T(2,2)+x(3)*T(3,2)+x(4)*T(4,2)+x(5)*T(5,2));
	y(3) = s*(x(0)*T(0,3)+x(1)*T(1,3)+x(2)*T(2,3)+x(3)*T(3,3)+x(4)*T(4,3)+x(5)*T(5,3));
	y(4) = s*(x(0)*T(0,4)+x(1)*T(1,4)+x(2)*T(2,4)+x(3)*T(3,4)+x(4)*T(4,4)+x(5)*T(5,4));
	y(5) = s*(x(0)*T(0,5)+x(1)*T(1,5)+x(2)*T(2,5)+x(3)*T(3,5)+x(4)*T(4,5)+x(5)*T(5,5));
}

//////////////// With update

/** 10) Dot product (with update): v += A*u  =>  v(i)+=A(i,k)*u(k) */
inline void DotUp (STensor2 const & A, Vector3D const & u,  Vector3D & v)
{
	v(0) += A(0)*u(0)     + A(3)*u(1)/SQ2 + A(5)*u(2)/SQ2;
	v(1) += A(3)*u(0)/SQ2 + A(1)*u(1)     + A(4)*u(2)/SQ2;
	v(2) += A(5)*u(0)/SQ2 + A(4)*u(1)/SQ2 + A(2)*u(2);
}

/** 10x) Scaled Dot product (with update): v += s*A*u  =>  v(i)+=s*A(i,k)*u(k) */
inline void ScDotUp (double s, STensor2 const & A, Vector3D const & u,  Vector3D & v)
{
	v(0) += s*( A(0)*u(0)     + A(3)*u(1)/SQ2 + A(5)*u(2)/SQ2 );
	v(1) += s*( A(3)*u(0)/SQ2 + A(1)*u(1)     + A(4)*u(2)/SQ2 );
	v(2) += s*( A(5)*u(0)/SQ2 + A(4)*u(1)/SQ2 + A(2)*u(2)     );
}

/** 10xx) Scaled Dot product (with update): v += s*A*u  =>  v(i)+=s*A(i,k)*u(k) */
inline void ScDotUp (double s, Vec_t const & A, Vector3D const & u,  Vector3D & v)
{
    if (size(A)==4)
    {
        v(0) += s*( A(0)*u(0)     + A(3)*u(1)/SQ2                 );
        v(1) += s*( A(3)*u(0)/SQ2 + A(1)*u(1)                     );
        v(2) += s*(                               + A(2)*u(2)     );
    }
    else
    {
        v(0) += s*( A(0)*u(0)     + A(3)*u(1)/SQ2 + A(5)*u(2)/SQ2 );
        v(1) += s*( A(3)*u(0)/SQ2 + A(1)*u(1)     + A(4)*u(2)/SQ2 );
        v(2) += s*( A(5)*u(0)/SQ2 + A(4)*u(1)/SQ2 + A(2)*u(2)     );
    }
}

/** 11) Dot product (with update): v += u*A  =>  v(i)+=u(k)*A(k,i) */
inline void DotUp (Vector3D const & u, STensor2 const & A,  Vector3D & v)
{
	v(0) += u(0)*A(0)     + u(1)*A(3)/SQ2 + u(2)*A(5)/SQ2;
	v(1) += u(0)*A(3)/SQ2 + u(1)*A(1)     + u(2)*A(4)/SQ2;
	v(2) += u(0)*A(5)/SQ2 + u(1)*A(4)/SQ2 + u(2)*A(2);
}

/** 11x) Scaled Dot product (with update): v += s*u*A  =>  v(i)+=u(k)*A(k,i) */
inline void ScDotUp (double s, Vector3D const & u, STensor2 const & A,  Vector3D & v)
{
	v(0) += s*( u(0)*A(0)     + u(1)*A(3)/SQ2 + u(2)*A(5)/SQ2 );
	v(1) += s*( u(0)*A(3)/SQ2 + u(1)*A(1)     + u(2)*A(4)/SQ2 );
	v(2) += s*( u(0)*A(5)/SQ2 + u(1)*A(4)/SQ2 + u(2)*A(2)     );
}

/** 12) Dot product: y += T*x  =>  y += T:x  =>  y(i,j)+=T(i,j,k,l)*x(k,l) */
inline void DotUp (STensor4 const & T, STensor2 const & x,  STensor2 & y)
{
	y(0) += T(0,0)*x(0)+T(0,1)*x(1)+T(0,2)*x(2)+T(0,3)*x(3)+T(0,4)*x(4)+T(0,5)*x(5);
	y(1) += T(1,0)*x(0)+T(1,1)*x(1)+T(1,2)*x(2)+T(1,3)*x(3)+T(1,4)*x(4)+T(1,5)*x(5);
	y(2) += T(2,0)*x(0)+T(2,1)*x(1)+T(2,2)*x(2)+T(2,3)*x(3)+T(2,4)*x(4)+T(2,5)*x(5);
	y(3) += T(3,0)*x(0)+T(3,1)*x(1)+T(3,2)*x(2)+T(3,3)*x(3)+T(3,4)*x(4)+T(3,5)*x(5);
	y(4) += T(4,0)*x(0)+T(4,1)*x(1)+T(4,2)*x(2)+T(4,3)*x(3)+T(4,4)*x(4)+T(4,5)*x(5);
	y(5) += T(5,0)*x(0)+T(5,1)*x(1)+T(5,2)*x(2)+T(5,3)*x(3)+T(5,4)*x(4)+T(5,5)*x(5);
}

/** 13) Dot product with addition: y = r + T*x  =>  y += T:x  =>  y(i,j)+=T(i,j,k,l)*x(k,l) */
inline void DotAdd (STensor2 const & r, STensor4 const & T, STensor2 const & x,  STensor2 & y)
{
	y(0) = r(0) + T(0,0)*x(0)+T(0,1)*x(1)+T(0,2)*x(2)+T(0,3)*x(3)+T(0,4)*x(4)+T(0,5)*x(5);
	y(1) = r(1) + T(1,0)*x(0)+T(1,1)*x(1)+T(1,2)*x(2)+T(1,3)*x(3)+T(1,4)*x(4)+T(1,5)*x(5);
	y(2) = r(2) + T(2,0)*x(0)+T(2,1)*x(1)+T(2,2)*x(2)+T(2,3)*x(3)+T(2,4)*x(4)+T(2,5)*x(5);
	y(3) = r(3) + T(3,0)*x(0)+T(3,1)*x(1)+T(3,2)*x(2)+T(3,3)*x(3)+T(3,4)*x(4)+T(3,5)*x(5);
	y(4) = r(4) + T(4,0)*x(0)+T(4,1)*x(1)+T(4,2)*x(2)+T(4,3)*x(3)+T(4,4)*x(4)+T(4,5)*x(5);
	y(5) = r(5) + T(5,0)*x(0)+T(5,1)*x(1)+T(5,2)*x(2)+T(5,3)*x(3)+T(5,4)*x(4)+T(5,5)*x(5);
}


// 14) s = xt * A * y
inline double Reduce (STensor2 const & x, STensor4 const & A, STensor2 const & y)
{
	return x(0)*(A(0,0)*y(0) + A(0,1)*y(1) + A(0,2)*y(2) + A(0,3)*y(3) + A(0,4)*y(4) + A(0,5)*y(5)) +
	       x(1)*(A(1,0)*y(0) + A(1,1)*y(1) + A(1,2)*y(2) + A(1,3)*y(3) + A(1,4)*y(4) + A(1,5)*y(5)) +
	       x(2)*(A(2,0)*y(0) + A(2,1)*y(1) + A(2,2)*y(2) + A(2,3)*y(3) + A(2,4)*y(4) + A(2,5)*y(5)) +
	       x(3)*(A(3,0)*y(0) + A(3,1)*y(1) + A(3,2)*y(2) + A(3,3)*y(3) + A(3,4)*y(4) + A(3,5)*y(5)) +
	       x(4)*(A(4,0)*y(0) + A(4,1)*y(1) + A(4,2)*y(2) + A(4,3)*y(3) + A(4,4)*y(4) + A(4,5)*y(5)) +
	       x(5)*(A(5,0)*y(0) + A(5,1)*y(1) + A(5,2)*y(2) + A(5,3)*y(3) + A(5,4)*y(4) + A(5,5)*y(5));
}

// 15) D = a * (A*x) dyadic (y*B) + C
inline void Generate (double const & a, STensor4 const & A, STensor2 const & x, STensor2 const & y, STensor4 const & B, STensor4 const & C, STensor4 & D)
{
	//C(i,j)=a* (A(i,k)*x(k)) * (y(l)*B(l,j)) +D(i,j);

	D(0,0)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(0,0);
	D(0,1)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(0,1);
	D(0,2)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(0,2);
	D(0,3)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(0,3);
	D(0,4)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(0,4);
	D(0,5)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(0,5);
	
	D(1,0)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(1,0);
	D(1,1)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(1,1);
	D(1,2)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(1,2);
	D(1,3)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(1,3);
	D(1,4)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(1,4);
	D(1,5)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(1,5);
	
	D(2,0)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(2,0);
	D(2,1)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(2,1);
	D(2,2)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(2,2);
	D(2,3)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(2,3);
	D(2,4)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(2,4);
	D(2,5)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(2,5);
	
	D(3,0)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(3,0);
	D(3,1)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(3,1);
	D(3,2)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(3,2);
	D(3,3)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(3,3);
	D(3,4)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(3,4);
	D(3,5)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(3,5);
	
	D(4,0)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(4,0);
	D(4,1)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(4,1);
	D(4,2)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(4,2);
	D(4,3)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(4,3);
	D(4,4)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(4,4);
	D(4,5)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(4,5);
	
	D(5,0)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(5,0);
	D(5,1)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(5,1);
	D(5,2)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(5,2);
	D(5,3)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(5,3);
	D(5,4)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(5,4);
	D(5,5)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(5,5);
}


///////////////////////////////////////////////////////////////////////////////////// Initialization /////

/** Initialize (global) tensors. */
int __initialize_tensors ()
{
	SymI    =  1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
	SymII   =  1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	           0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
	           0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
	           0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
	           0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
	           0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
	SymIdyI =  1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
	           1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
	           1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
	           0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	           0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	           0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	SymPsd  =  2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
	          -1.0/3.0,  2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
	          -1.0/3.0, -1.0/3.0,  2.0/3.0, 0.0, 0.0, 0.0,
	               0.0,      0.0,      0.0, 1.0, 0.0, 0.0,
	               0.0,      0.0,      0.0, 0.0, 1.0, 0.0,
	               0.0,      0.0,      0.0, 0.0, 0.0, 1.0;
	SymPiso =  1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
	           1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
	           1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
	               0.0,     0.0,     0.0, 0.0, 0.0, 0.0,
	               0.0,     0.0,     0.0, 0.0, 0.0, 0.0,
	               0.0,     0.0,     0.0, 0.0, 0.0, 0.0;
	return 0;
}

/** Dummy variable used only for tensors initialization. */
int __dummy_initialize_tensors = __initialize_tensors ();

}; // namespace MPM

#endif // MPM_TENSORS_H
