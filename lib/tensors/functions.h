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

#ifndef MECHSYS_TENSORS_FUNCTIONS_H
#define MECHSYS_TENSORS_FUNCTIONS_H

// STL
#include <cmath>  // for sqrt()
#include <cfloat> // for DBL_EPSILON
#include <cstring> // for strcmp

// MechSys
#include "tensors/tensors.h"
#include "tensors/operators.h"
#include "tensors/jacobirot.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "util/util.h"

using Util::SQ3;
using Util::SQ6;
using Util::SQ2BY3;

namespace Tensors
{

inline double Tr(Tensor2 const & A) { return A(0)+A(1)+A(2); }

inline double Val(Tensor2 const & Sig, Tensor2 const & Eps, char const * Name)
{
	     if (strcmp(Name,"Sx" )==0)                          return Sig(0);
	else if (strcmp(Name,"Sy" )==0)                          return Sig(1);
	else if (strcmp(Name,"Sz" )==0)                          return Sig(2);
	else if (strcmp(Name,"Sxy")==0 || strcmp(Name,"Syx")==0) return Sig(3)/SQ2;
	else if (strcmp(Name,"Syz")==0 || strcmp(Name,"Szy")==0) return Sig(4)/SQ2;
	else if (strcmp(Name,"Szx")==0 || strcmp(Name,"Sxz")==0) return Sig(5)/SQ2;
	else if (strcmp(Name,"p"  )==0)                          return (Sig(0)+Sig(1)+Sig(2))/3.0;
	else if (strcmp(Name,"q"  )==0)                          return sqrt(((Sig(0)-Sig(1))*(Sig(0)-Sig(1)) + (Sig(1)-Sig(2))*(Sig(1)-Sig(2)) + (Sig(2)-Sig(0))*(Sig(2)-Sig(0)) + 3.0*(Sig(3)*Sig(3) + Sig(4)*Sig(4) + Sig(5)*Sig(5)))/2.0);
	else if (strcmp(Name,"Ex" )==0)                          return Eps(0);
	else if (strcmp(Name,"Ey" )==0)                          return Eps(1);
	else if (strcmp(Name,"Ez" )==0)                          return Eps(2);
	else if (strcmp(Name,"Exy")==0 || strcmp(Name,"Eyx")==0) return Eps(3)/SQ2;
	else if (strcmp(Name,"Eyz")==0 || strcmp(Name,"Ezy")==0) return Eps(4)/SQ2;
	else if (strcmp(Name,"Ezx")==0 || strcmp(Name,"Exz")==0) return Eps(5)/SQ2;
	else if (strcmp(Name,"Ev" )==0)                          return Eps(0)+Eps(1)+Eps(2); 
	else if (strcmp(Name,"Ed" )==0)                          return sqrt(2.0*((Eps(0)-Eps(1))*(Eps(0)-Eps(1)) + (Eps(1)-Eps(2))*(Eps(1)-Eps(2)) + (Eps(2)-Eps(0))*(Eps(2)-Eps(0)) + 3.0*(Eps(3)*Eps(3) + Eps(4)*Eps(4) + Eps(5)*Eps(5))))/3.0;
	else throw new Fatal("Tensors::Val: Name==%s is invalid",Name);
}

inline void SetVal(char const * Name, double Val, Tensor2 & Sig, bool WithError=true)
{
	     if (strcmp(Name,"ZERO")==0)                          return;
	else if (strcmp(Name,"Sx"  )==0)                          Sig(0) = Val;
	else if (strcmp(Name,"Sy"  )==0)                          Sig(1) = Val;
	else if (strcmp(Name,"Sz"  )==0)                          Sig(2) = Val;
	else if (strcmp(Name,"Sxy" )==0 || strcmp(Name,"Syx")==0) Sig(3) = Val*SQ2;
	else if (strcmp(Name,"Syz" )==0 || strcmp(Name,"Szy")==0) Sig(4) = Val*SQ2;
	else if (strcmp(Name,"Szx" )==0 || strcmp(Name,"Sxz")==0) Sig(5) = Val*SQ2;
	else if (WithError) throw new Fatal("Tensors::SetVal: Name==%s is invalid",Name);
}

/** Tensor multiplication  \f$ \TeSe{R} \gets \TeSe{S}\bullet\TeSe{T} \f$.
 * \em IMPORTANT
 *
 *  This function is valid only if the result R=S*T is known to be symmetric,
 *  however this is \em NOT always true, even if S and T are both symmetric !!!
 */
inline void Mult(Tensor2 const & S, Tensor2 const & T, Tensor2 & R)
{
	const double ZERO = sqrt(DBL_EPSILON);

	R(0) = S(0)*T(0)     + S(3)*T(3)/2.0 + S(5)*T(5)/2.0;
	R(1) = S(3)*T(3)/2.0 + S(1)*T(1)     + S(4)*T(4)/2.0;
	R(2) = S(5)*T(5)/2.0 + S(4)*T(4)/2.0 + S(2)*T(2);

	R(3) = S(0)*T(3)     + S(3)*T(1)     + S(5)*T(4)/SQ2;
	R(4) = S(3)*T(5)/SQ2 + S(1)*T(4)     + S(4)*T(2);
	R(5) = S(0)*T(5)     + S(3)*T(4)/SQ2 + S(5)*T(2);

#ifndef NDEBUG
	double A = S(3)*T(0)     + S(1)*T(3)     + S(4)*T(5)/SQ2;
	double B = S(5)*T(3)/SQ2 + S(4)*T(1)     + S(2)*T(4);
	double C = S(5)*T(0)     + S(4)*T(3)/SQ2 + S(2)*T(5);
	//std::cout << _15_8 R(3)-A << " | " << _15_8 R(4)-B << " | " << _15_8 R(5)-C << std::endl;
	if (!(fabs(R(3)-A)<ZERO && fabs(R(4)-B)<ZERO && fabs(R(5)-C)<ZERO))
		throw new Fatal(_("Tensors::Mult: Result of matrix multiplication is NOT symmetric.\n (fabs(R(3)-A)<0 & fabs(R(4)-B)<0 & fabs(R(5)-C)<0) failed.\n S=<%.6f %.6f %.6f %.6f %.6f %.6f>,\n R=<%.6f %.6f %.6f %.6f %.6f %.6f>"),
		                S(0), S(1), S(2), S(3)/SQ2, S(4)/SQ2, S(5)/SQ2,
		                R(0), R(1), R(2), R(3)/SQ2, R(4)/SQ2, R(5)/SQ2);
#endif

}

/** Euclidian norm of a symmetric 2nd order tensor. */
inline double Norm(Tensor2 const & T)
{

	return sqrt(T(0)*T(0) + T(1)*T(1) + T(2)*T(2) + T(3)*T(3) + T(4)*T(4) + T(5)*T(5));

}

/** Determinant of T. */
inline double Det(Tensor2 const & T)
{

	return T(0)*T(1)*T(2) + T(3)*T(4)*T(5)/SQ2 - T(0)*T(4)*T(4)/2.0 - T(1)*T(5)*T(5)/2.0 - T(2)*T(3)*T(3)/2.0;

}

/** Inverse of a tensor T  \f$ \TeSe{R} \gets \TeSe{T}^{-1} \f$. */
inline void Inv(Tensor2 const & T, Tensor2 & R)
{
	const double ZERO = sqrt(DBL_EPSILON);

	double det = Det(T);
	
	if (fabs(det)<ZERO)
	   throw new Fatal(_("Tensors::Inv: For the computation of the inverse of a tensor, the determinant must be non null. T=<%.6f %.6f %.6f %.6f %.6f %.6f>"),
		                T(0), T(1), T(2), T(3)/SQ2, T(4)/SQ2, T(5)/SQ2);

	R(0) =     ( T(1)*T(2)     - T(4)*T(4)/2.0 )/det;
	R(1) =     ( T(2)*T(0)     - T(5)*T(5)/2.0 )/det;
	R(2) =     ( T(0)*T(1)     - T(3)*T(3)/2.0 )/det;
	R(3) = SQ2*( T(4)*T(5)/2.0 - T(2)*T(3)/SQ2 )/det;
	R(4) = SQ2*( T(5)*T(3)/2.0 - T(0)*T(4)/SQ2 )/det;
	R(5) = SQ2*( T(3)*T(4)/2.0 - T(1)*T(5)/SQ2 )/det;

}

/** Eigenvalues (L) of a tensor T  \f$ L_{i=1\cdots3} \gets eigenvalues(\TeSe{T}) \f$. */
inline void Eigenvals(Tensor2 const & T, Tensor1 & L)
{
	if (JacobiRot(T,L)==-1) throw new Fatal("Tensors::Eigenvals: Jacobi rotation did not converge. T=<%.6f %.6f %.6f %.6f %.6f %.6f>", T(0), T(1), T(2), T(3)/SQ2, T(4)/SQ2, T(5)/SQ2);
}

// Sort descending
inline void Sort(Tensor1 & L)
{
	for (int i=0; i<2; ++i)
	{
		// Find the index of the maximum element
		int max_index = i;
		for (int j=i+1; j<3; ++j) if (L(j)>L(max_index)) max_index = j;

		// Swap if i-th element is not the biggest yet
		if (max_index>i) Util::Swap(L(i), L(max_index));
	}
}

/** Eigenvalues (L) and Eigenprojectors (P) of a tensor T \f$ L_{i=1\cdots3} \gets eigenvalues(\TeSe{T}) \quad \TeSe{P}_{i=1\cdots3} \gets eigenprojectors(\TeSe{T}) \f$.
 *  Any symmetric second order tensor \f$\TeSe{T}\f$ may be calculated according to its \em Spectral \em Decomposition:
 * \f[
 *     \TeSe{T} = L_0\TeFi{V}_0\otimes\TeFi{V}_0 + L_1\TeFi{V}_1\otimes\TeFi{V}_1 + L_2\TeFi{V}_2\otimes\TeFi{V}_2
 * \f]
 * in which \f$L_k\f$ and \f$\TeFi{V}_k\f$ are the eigenvalues and eigenvectors of tensor \f$\TeSe{T}\f$, respectively
 *
 * The eigenprojector is defined according to:
 * \f[
 *     \TeSe{P}_k = \TeFi{V}_{(k)} \otimes \TeFi{V}_{(k)} \qquad no \quad sum \quad on \quad (k)
 * \f]
 *
 * Therefore, the spectral representation can be re-written as:
 * \f[
 *     \TeSe{T} = L_0\TeSe{P}_0 + L_1\TeSe{P}_1 + L_2\TeSe{P}_2
 * \f]
 *
 * Or, considering Einstein's summation rule (over index k, from 1 to 3)
 * \f[
 *     \TeSe{T} = L_k \TeSe{P}_k
 * \f]
 *
 * The following properties hold for the eigenprojectors, in fact, for any \em projector
 * (<a href="http://www.mech.utah.edu/~brannon/gobag.html">Brannon, 2000</a>),
 * \f[ \TeSe{P}_i \bullet \TeSe{P}_j = \TeSe{P}_i  \quad if \quad i = j   \f]
 * \f[ \TeSe{P}_i \bullet \TeSe{P}_j = \TeSe{0}    \quad if \quad i \ne j \f]
 * \f[ \TeSe{P}_1 + \TeSe{P}_2 + \TeSe{P}_3 = \TeSe{I}                    \f]
 * 
 * \return true if success, otherwise return false
 */
inline void Eigenvp(Tensor2 const & T, double L[3], Tensor2 P[3])
{

	double V0[3]; // eigenvector
	double V1[3]; // eigenvector
	double V2[3]; // eigenvector

	// Calculate eigenvalues and eigenvectors
	if (JacobiRot(T,V0,V1,V2,L)==-1)
	   throw new Fatal(_("Tensors::Eigenvp: Jacobi rotation did not converge. T=<%.6f %.6f %.6f %.6f %.6f %.6f>"),
		                T(0), T(1), T(2), T(3)/SQ2, T(4)/SQ2, T(5)/SQ2);

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

/** Square root of a tensor T  \f$ \TeSe{R} \gets \sqrt{\TeSe{T}} \f$. */
inline void Sqrt(Tensor2 const & T, Tensor2 & R)
{

	double  L[3]; // Eigenvalues
	Tensor2 P[3]; // Eigenprojectors

	// Calculate eigenvalues and eigenprojectors
	Eigenvp(T,L,P);
	
	// Calculate R = sqrt(T)
	R(0)=0.0;
	R(1)=0.0;
	R(2)=0.0;
	R(3)=0.0;
	R(4)=0.0;
	R(5)=0.0;
	for (int i=0; i<3; ++i)
		R = sqrt(L[i])*P[i] + R; // R <- sqrt(L[i])*P[i] + R

}

/** Characteristic invariants of a symmetric second order tensor. */
inline void CharInvs(Tensor2 const & T, double I[3])
{

	I[0] = T(0)+T(1)+T(2);
	I[1] = T(0)*T(1) + T(1)*T(2) + T(2)*T(0) - (T(3)*T(3) + T(4)*T(4) + T(5)*T(5))/2.0;
	I[2] = T(0)*T(1)*T(2) + T(3)*T(4)*T(5)/SQ2 - (T(0)*T(4)*T(4) + T(1)*T(5)*T(5) + T(2)*T(3)*T(3))/2.0;

}

/** Strain Invariants \f$ (\varepsilon_v, \varepsilon_d) \f$.
 *  \f[ Eps=\TeSe{\varepsilon} \f]
 *  \f[ Ev=\varepsilon_v=tr(\TeSe{\varepsilon}) \f]
 *  \f[ Ed=\varepsilon_d=\sqrt{\nicefrac{2}{3}\TeSe{e}:\TeSe{e}} \quad \TeSe{e}=dev(\TeSe{\varepsilon})=\TeSe{\varepsilon}-\varepsilon_V\TeSe{I}/3 \f]
 *  \param Eps Strain tensor
 *  \param Ev  Volumetric strain invariant
 *  \param Ed  Deviatoric strain invariant
 */
inline void Strain_Ev_Ed(Tensor2 const & Eps, double & Ev, double & Ed)
{

	Ev = Eps(0)+Eps(1)+Eps(2);
	Ed = sqrt(2.0*(   (Eps(0)-Eps(1))*(Eps(0)-Eps(1))
				    + (Eps(1)-Eps(2))*(Eps(1)-Eps(2))
					+ (Eps(2)-Eps(0))*(Eps(2)-Eps(0)) + 3.0*(Eps(3)*Eps(3) + Eps(4)*Eps(4) + Eps(5)*Eps(5)) ))/3.0;

}

/** Stress Invariants (Cambridge) \f$ (p,q) \f$.
 *  \f[ Sig=\TeSe{\sigma} \f]
 *  \f[ p=tr(\TeSe{\sigma})/3 \f]
 *  \f[ q=\sqrt{\nicefrac{3}{2}\TeSe{S}:\TeSe{S}} \quad \TeSe{S}=dev(\TeSe{\sigma})=\TeSe{\sigma}-p\TeSe{I} \f]
 *  \param Sig Stress tensor
 *  \param p   Cambridge mean stress
 *  \param q   Cambridge deviatoric stress
 */
inline void Stress_p_q(Tensor2 const & Sig, double & p, double & q)
{

	p = (Sig(0)+Sig(1)+Sig(2))/3.0;
	q = sqrt((   (Sig(0)-Sig(1))*(Sig(0)-Sig(1))
			   + (Sig(1)-Sig(2))*(Sig(1)-Sig(2))
			   + (Sig(2)-Sig(0))*(Sig(2)-Sig(0)) + 3.0*(Sig(3)*Sig(3) + Sig(4)*Sig(4) + Sig(5)*Sig(5)) )/2.0);

}

/** Cambridge's q deviatoric stress invariant \f$ (q) \f$.
 * \param Sig Stress tensor
 * \return Cambridge's q deviatoric stress invariant
 */
inline double Stress_q(Tensor2 const & Sig)
{

	return sqrt((   (Sig(0)-Sig(1))*(Sig(0)-Sig(1))
	              + (Sig(1)-Sig(2))*(Sig(1)-Sig(2))
	              + (Sig(2)-Sig(0))*(Sig(2)-Sig(0)) + 3.0*(Sig(3)*Sig(3) + Sig(4)*Sig(4) + Sig(5)*Sig(5)) )/2.0);

}

/** Sin3Th given deviator (S) of the stress tensor.
 * \param S Deviator of the stress tensor
 * \return sin(3*th)
 *
 * \verbatim
 *              SI
 *        alp:  0
 *         th: +30
 *              |         A: alp: alpha = 30 - theta
 *        alp <-|         T: th : theta = 30 - alpha
 *              |
 *       30     |         Sa = -r*sin(alp)
 *   60   0-> th|         Sb =  r*cos(alp)
 *   -30   \    |         Sc =  p*sqrt(3)
 *      ,   \   ^ Sb
 *       ',  \  |
 *         '. \ |
 *           '.(o).-----> Sa
 *          .' Sc  '.
 *        .'         '.
 *   SII.'             '. SIII
 * \endverbatim
 */
inline double Sin3ThDev(Tensor2 const & S)
{
	const double ZERO = sqrt(DBL_EPSILON);

	double normS = Norm(S);
	if (normS<=ZERO) return -1.0;
	else
	{
		double   L = 9.0*SQ2*Det(S)/(SQ3*pow(normS,3.0));
		if      (L>= 1.0) return  1.0;
		else if (L<=-1.0) return -1.0;
		else              return  L;
	}

}

/** Stress Invariants (Cambridge) + deviator of Sigma.
 * \param Sig    Stress tensor
 * \param p      Cambridge mean stress
 * \param q      Cambridge deviatoric stress
 * \param S      Deviator of stress tensor
 * \param sin3th Sine of 3*Lode angle
 */
inline void Stress_p_q_S_sin3th(Tensor2 const & Sig, double & p, double & q, Tensor2 & S, double & sin3th)
{
	const double ZERO = sqrt(DBL_EPSILON);

	Stress_p_q(Sig,p,q);
	S = Sig - I * p;
	if (q<=ZERO)  sin3th = -1.0; // 3th==-90 => th=30
	else          sin3th = 27.0*Det(S)/(2.0*pow(q,3.0));
	if (sin3th>= 1.0) sin3th =  1.0;
	if (sin3th<=-1.0) sin3th = -1.0;

}

/** Stress Invariants (Cambridge) + deviator of Sigma.
 * \param Sig    Stress tensor
 * \param p      Cambridge mean stress
 * \param q      Cambridge deviatoric stress
 * \param S      Deviator of stress tensor
 * \param t      27*det(S)/(2*q^3)
 */
inline void Stress_p_q_S_t(Tensor2 const & Sig, double & p, double & q, Tensor2 & S, double & t)
{
	const double MINQ = 1.0e-5;
	Stress_p_q (Sig,p,q);
	S = Sig - I * p;
	t = (q>MINQ ? 27.0*Det(S)/(2.0*pow(q,3.0)) : 0.0);
}

/** Stress Invariants (Cambridge) + deviator of Sigma and gradients
 * \param Sig    Stress tensor
 * \param p      Cambridge mean stress
 * \param q      Cambridge deviatoric stress
 * \param S      Deviator of stress tensor
 * \param t      27*det(S)/(2*q^3)
 * \param dqds   diff(q,Sig)
 * \param dtds   diff(t,Sig)
 */
inline void Stress_p_q_S_t(Tensor2 const & Sig, double & p, double & q, Tensor2 & S, double & t, Tensor2 & dqds, Tensor2 & dtds)
{
	const double MINQ = 1.0e-5;
	Stress_p_q (Sig,p,q);
	S    = Sig - I * p;
	t    = 0.0;
	dqds = 0.0;
	dtds = 0.0;
	if (q>MINQ)
	{
		Tensor2 SS;    Mult (S,S, SS);
		Tensor2 devSS; Dot  (Psd,SS, devSS);
		t    = 27.0*Det(S)/(2.0*pow(q,3.0));
		dqds = (3.0/(2.0*q))*S;
		dtds = (9.0/(2.0*q*q))*((3.0/q)*devSS-t*S);
	}
}

/** Stress Invariants (Cambridge) + deviator of Sigma and gradients + high-order gradients
 * \param Sig     Stress tensor
 * \param p       Cambridge mean stress
 * \param q       Cambridge deviatoric stress
 * \param S       Deviator of stress tensor
 * \param t       27*det(S)/(2*q^3)
 * \param dpds    diff(p,Sig)
 * \param dqds    diff(q,Sig)
 * \param dtds    diff(t,Sig)
 * \param d2qds2  diff(q,Sig, 2)
 * \param d2tds2  diff(t,Sig, 2)
 */
inline void Stress_p_q_S_t(Tensor2 const & Sig, double & p, double & q, Tensor2 & S, double & t, Tensor2 & dqds, Tensor2 & dtds, Tensor4 & d2qds2, Tensor4 & d2tds2)
{
	const double MINQ = 1.0e-5;
	Stress_p_q (Sig,p,q);
	S      = Sig - I * p;
	t      = 0.0;
	dqds   = 0.0;
	dtds   = 0.0;
	d2qds2 = 0.0;
	d2tds2 = 0.0;
	if (q>MINQ)
	{
		double q2 = q*q;
		double q3 = q2*q;
		double q4 = q3*q;
		double c  = 9.0/(2.0*q2);
		double d  = 3.0/(2.0*q);
		Tensor2 SS;    Mult (S,S, SS);
		Tensor2 devSS; Dot  (Psd,SS, devSS);
		Tensor4 SdyS;  Dyad (S,S, SdyS);
		t    = 27.0*Det(S)/(2.0*q3);
		dqds = d*S;
		dtds = c*((3.0/q)*devSS-t*S);
		AddScaled (d,Psd, -(9.0/(4.0*q3)),SdyS, d2qds2);
		Tensor2 tmp1; tmp1 = (9.0*t/q3)*dqds - c*dtds;
		Tensor4 T,tmp2,tmp3,tmp4;
		Dyad        ((2.0/3.0)*I,S, tmp2);
		LeafLeafPsd (I,S, T);
		AddScaledUp (-1.0,tmp2, T);
		Dyad        (devSS,dqds, tmp3);
		Dyad        (S,tmp1, tmp4);
		AddScaled   (27.0/(2.0*q3),T ,-81.0/(2.0*q4),tmp3, -t*c,Psd, 1.0,tmp4, d2tds2);
	}
}

/** Stress Invariants (Professor Nakai's invariants) \f$ (t_N,t_S) \f$.
 *  \f[ Sig=\TeSe{\sigma} \f]
 *  \f[ t_N=3\frac{I_3}{I2} \f]
 *  \f[ t_S=\frac{\sqrt{I_1I_2I_3-9I_3^2}}{I_2} \f]
 * \param Sig Stress tensor
 * \param tn  Mean stress
 * \param ts  Deviatoric stress
 */
inline void Stress_tn_ts(Tensor2 const & Sig, double & tn, double & ts)
{
	const double ZERO = sqrt(DBL_EPSILON);

	double I[3];
	CharInvs(Sig,I);
	if (fabs(I[1])<=ZERO)
	   throw new Fatal(_("Tensors::Stress_tn_ts: Second stress characteristic invariant must be non zero. Sig=<%.6f %.6f %.6f %.6f %.6f %.6f>"),
		                Sig(0), Sig(1), Sig(2), Sig(3)/SQ2, Sig(4)/SQ2, Sig(5)/SQ2);

	double s = I[0]*I[1]*I[2]-9.0*I[2]*I[2];
	if (s<0.0)
		throw new Fatal(_("Tensors::Stress_tn_ts: I1*I2*I3-9*I3^2 must be positive. Sig=<%.6f %.6f %.6f %.6f %.6f %.6f>"),
		                Sig(0), Sig(1), Sig(2), Sig(3)/SQ2, Sig(4)/SQ2, Sig(5)/SQ2);

	tn = 3.0*I[2]/I[1];
	if (tn<ZERO)
		throw new Fatal(_("Tensors::Stress_tn_ts: tn must be greater than/equal to 0 (tn must be >= 0). Sig=<%.6f %.6f %.6f %.6f %.6f %.6f>"),
		                Sig(0), Sig(1), Sig(2), Sig(3)/SQ2, Sig(4)/SQ2, Sig(5)/SQ2);

	ts = sqrt(s)/I[1];
	if (ts!=ts)
		throw new Fatal(_("Tensors::Stress_tn_ts: tn is NaN. Sig=<%.6f %.6f %.6f %.6f %.6f %.6f>"),
		                Sig(0), Sig(1), Sig(2), Sig(3)/SQ2, Sig(4)/SQ2, Sig(5)/SQ2);

}

/** Stress Invariants (Professor Brannon's isomorphic invariants) \f$ (P,Q) \f$.
 *  \f[ Sig=\TeSe{\sigma} \f]
 *  \f[ P=\sqrt{3}p \f]
 *  \f[ Q=\sqrt{\TeSe{S}:\TeSe{S}}  \quad \TeSe{S}=dev(\TeSe{\sigma})=\TeSe{\sigma}-p\TeSe{I} \f]
 * \param Sig Stress tensor
 * \param P   Mean stress invariant
 * \param Q   Deviatoric stress invariant
 */
inline void Stress_P_Q(Tensor2 const & Sig, double & P, double & Q)
{

	P = (Sig(0)+Sig(1)+Sig(2))/SQ3;
	Q = sqrt( ((Sig(0)-Sig(1))*(Sig(0)-Sig(1))
	         + (Sig(1)-Sig(2))*(Sig(1)-Sig(2))
	         + (Sig(2)-Sig(0))*(Sig(2)-Sig(0)))/3.0 + Sig(3)*Sig(3) + Sig(4)*Sig(4) + Sig(5)*Sig(5) );

}

/** Stress Invariants (Professor Brannon's isomorphic invariants) \f$ (P,Q) \f$. */
inline void Stress_P_Q_S_sin3th(Tensor2 const & Sig, double & P, double & Q, Tensor2 & S, double & sin3th)
{
	const double ZERO = sqrt(DBL_EPSILON);

	Stress_P_Q(Sig,P,Q);
	S = Sig - I * (P/SQ3);
	if (Q<=ZERO)  sin3th = -1.0; // 3th==-90 => th=30
	else          sin3th = 27.0*Det(S)/(2.0*pow(SQ3*Q/SQ2,3.0));
	if (sin3th>= 1.0) sin3th =  1.0;
	if (sin3th<=-1.0) sin3th = -1.0;

}

/** Returns Sin3Th, given three principal values, which are not necessary sorted. */
inline double Sin3Th(double SI, double SII, double SIII)
{

	double p       = (SI+SII+SIII)/3.0;
	double devS[3] = {SI-p, SII-p, SIII-p};
	double det     = devS[0]*devS[1]*devS[2];
	double norm    = sqrt(devS[0]*devS[0] + devS[1]*devS[1] + devS[2]*devS[2]);
	if (norm<=1.0e-5) return -1.0;
	else
	{
		double   L = 9.0*SQ2*det/(SQ3*pow(norm,3.0));
		if      (L>= 1.0) return  1.0;
		else if (L<=-1.0) return -1.0;
		else              return  L;
	}

}

/** Converts hydrostatic coordinates to sigma (principal) coord (A in radians). */
inline void Hid2Sig(double const &    p, double const &    q, double const &  alp,
		            double       & Sig1, double       & Sig2, double       & Sig3)
{

	double r    = q*SQ2BY3;
	double Sx   = -r*sin(alp);
	double Sy   =  r*cos(alp);
	       Sig1 = p + 2.0*Sy/SQ6;
	       Sig2 = p -     Sy/SQ6 - Sx/SQ2;
	       Sig3 = p -     Sy/SQ6 + Sx/SQ2;

}

/** Converts hydrostatic coordinates to sigma (principal) coord (A in radians).
 * NOTE:
 *    SI, SII and SIII must be previously allocated !!!
 */
inline void Hid2Sig(double const *  P, double const *   Q, double const *    A,
		            double       * SI, double       * SII, double       * SIII, int Size)
{

	for (int i=0; i<Size; ++i)
	{
		double r  = Q[i]*SQ2BY3;
		double Sx = -r*sin(A[i]);
		double Sy =  r*cos(A[i]);
		SI  [i] = P[i] + 2.0*Sy/SQ6;
		SII [i] = P[i] -     Sy/SQ6 - Sx/SQ2;
		SIII[i] = P[i] -     Sy/SQ6 + Sx/SQ2;
	}

}

/** Converts hydrostatic coordinates (Sx, Sy, Sz) to sigma (principal). */
inline void Hid2Sig_(double const & SX, double const &  SY, double const &    p,
		             double       & SI, double       & SII, double       & SIII)
{

	SI   = p + 2.0*SY/SQ6;
	SII  = p -     SY/SQ6 - SX/SQ2;
	SIII = p -     SY/SQ6 + SX/SQ2;

}

/** Professor Yao's Stress Invariants.
 * \param Sig   Stress tensor
 * \param p_yao Yao's modified Cambridge mean stress
 * \param q_yao Yao's modified Cambridge deviatoric stress
 * \param S_yao Yao's modified deviator of Sigma
 */
inline void Stress_p_q_S_yao(Tensor2 const & Sig, double & p_yao, double & q_yao, Tensor2 & S_yao)
{
	const double ZERO = sqrt(DBL_EPSILON);

	double Is[3];
	Tensors::CharInvs(Sig,Is); // I1=Is[0]; I2=Is[1]; I3=Is[2];
	           p_yao = (Sig(0)+Sig(1)+Sig(2))/3.0;
	Tensor2 S; S     = Sig - p_yao*I;
	double     s     = sqrt(Tensors::Dot(S,S));
	if (s>ZERO)
	{
		double a     = (Is[0]*Is[1]-Is[2])/(Is[0]*Is[1]-9.0*Is[2]);
		double l0    = 2.0*SQ2*Is[0]/(SQ3*(3.0*sqrt(a)-1.0));
		       S_yao = (l0/s)*S;
			   q_yao = sqrt(1.5)*sqrt(Tensors::Dot(S_yao,S_yao));
	}
	else q_yao = 0.0;

}

/** Output stress and strain values and invariants for use with plotn in plot.R - HEADER. */
inline void OutputHeader (std::ostream & Os)
{
	Os << Util::_8s<< "Sx" << Util::_8s<< "Sy" << Util::_8s<< "Sz" << Util::_8s<< "Sxy" << Util::_8s<< "Syz" << Util::_8s<< "Sxz";
	Os << Util::_8s<< "Ex" << Util::_8s<< "Ey" << Util::_8s<< "Ez" << Util::_8s<< "Exy" << Util::_8s<< "Eyz" << Util::_8s<< "Exz";
	Os << Util::_8s<< "p"  << Util::_8s<< "q"  << Util::_8s<< "Ev" << Util::_8s<< "Ed";
}

/** Output stress and strain values and invariants for use with plotn in plot.R. */
inline void Output (Tensor2 const & Sig, Tensor2 const & Eps,  std::ostream & Os)
{
	double p,q,ev,ed;
	Tensors::Stress_p_q   (Sig,p,q);
	Tensors::Strain_Ev_Ed (Eps,ev,ed);
	Os << Util::_8s<< Sig(0) << Util::_8s<< Sig(1) << Util::_8s<< Sig(2) << Util::_8s<< Sig(3)/SQ2 << Util::_8s<< Sig(4)/SQ2 << Util::_8s<< Sig(5)/SQ2;
	Os << Util::_8s<< Eps(0) << Util::_8s<< Eps(1) << Util::_8s<< Eps(2) << Util::_8s<< Eps(3)/SQ2 << Util::_8s<< Eps(4)/SQ2 << Util::_8s<< Eps(5)/SQ2;
	Os << Util::_8s<< p      << Util::_8s<< q      << Util::_8s<< ev*100.0<< Util::_8s<< ed*100.0;
}


}; // namespace Tensors

#endif // MECHSYS_TENSORS_FUNCTIONS_H
