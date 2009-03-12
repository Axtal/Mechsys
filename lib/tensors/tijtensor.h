/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo, Fernando Alonso                   *
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

#ifndef MECHSYS_TIJTENSOR_H
#define MECHSYS_TIJTENSOR_H

// STL
#include <cmath>

// MechSys
#include "tensors/tensors.h"
#include "tensors/functions.h"
#include "util/exception.h"

using Tensors::Tensor2;
using Tensors::CharInvs;
using Tensors::Stress_tn_ts;
using Tensors::Sqrt;
using Tensors::Inv;
using Tensors::Mult;
using Tensors::I;

namespace Tensors
{

/// Prof. Nakai's modified stress tensor tij
class TijTensor
{
public:
	/** tij tensor.
	 * \param TsMin Minimum allowed value for ts (derivatives related)
	 */
	TijTensor(double TsMin) : _ts_min(TsMin) {}

	/** Derivatives and tij tensor related quantities.
	 * \param Sig       In: Stress tensor \f$ \TeSe{\sigma} \f$
	 * \param tn        Out: Mean stress invariant \f$ t_N \gets 3\frac{I_3}{I2}\f$
	 * \param ts        Out: Deviatoric stress invariant \f$ t_S \gets \frac{\sqrt{I_1I_2I_3-9I_3^2}}{I_2}\f$
	 * \param Is        Out: Characterist invariants of stress (Sig) \f$ I_1,I_2,I_3 \f$
	 * \param sqrt_rz   Out: Sqrt of I3/I2 \f$ \sqrt{\frac{I_3}{I_2}} \f$
	 * \param Tau       Out: Tensor Sqrt of Sig \f$ \TeSe{\tau} \gets \sqrt{\TeSe{\sigma}} \f$
	 * \param inv_Tau   Out: Inverse of tensor Tau \f$ \TeSe{\tau}^{-1} \f$
	 * \param a         Out: "normal-to-SMP" tensor == Derivative of tn w.r.t tij tensor \f$ \TeSe{a} \gets \sqrt{\frac{I_3}{I_2}}\TeSe{\tau}^{-1} \equiv \pderiv{t_N}{\TeSe{t}} \f$
	 * \param t         Out: Modified tij tensor \f$ \TeSe{t} \gets \TeSe{a}\bullet\TeSe{\sigma} \f$
	 * \param dtn_ds    Out: Derivative of tn w.r.t Sig \f$ \pderiv{t_N}{\TeSe{\sigma}} \gets \frac{-\tn}{I_2}(I_1\TeSe{I}-\TeSe{\sigma}) + \frac{3}{I_2}(\TeSe{\sigma}\bullet\TeSe{\sigma}-I_1\TeSe{\sigma}+I_2\TeSe{I}) \f$
	 * \param dts_ds    Out: Derivative of ts w.r.t Sig \f$ \pderiv{t_S}{\TeSe{\sigma}} \gets  \frac{t_N}{6t_S}\TeSe{I}+\frac{t_N(6t_N-I_1)}{6t_S I_2}(I_1\TeSe{I}-\TeSe{\sigma}) + \frac{I_1-6t_N}{2t_S I_2} (\TeSe{\sigma}\bullet\TeSe{\sigma}-I_1\TeSe{\sigma}+I_2\TeSe{I}) \f$
	 * \param dts_dt    Out: Derivative of ts w.r.t tij tensor \f$ \pderiv{t_S}{\TeSe{t}} \gets \frac{\TeSe{t}-t_N\TeSe{a}}{t_S} \f$
	 */
	void Derivs(Tensor2 const & Sig,
	            double        & tn,
	            double        & ts,
		        double          Is[3],
	            double        & sqrt_rz,
	            Tensor2       & Tau,
	            Tensor2       & inv_Tau,
	            Tensor2       & a,
	            Tensor2       & t,
	            Tensor2       & dtn_ds,
	            Tensor2       & dts_ds,
	            Tensor2       & dts_dt);
private:
	double _ts_min; ///< Mininum allowed value for ts (derivatives related)
}; // class TijTensor


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void TijTensor::Derivs(Tensor2 const & Sig,
	                          double        & tn, 
	                          double        & ts, 
							  double          Is[3],
	                          double        & sqrt_rz, 
	                          Tensor2       & Tau,
	                          Tensor2       & inv_Tau,
	                          Tensor2       & a,
	                          Tensor2       & t,
	                          Tensor2       & dtn_ds,
	                          Tensor2       & dts_ds,
	                          Tensor2       & dts_dt)
{
	// ---------------------------------------------------------------------- tij calculation
	
	// Characteristic invariants:  I1=Is[0], I2=Is[1], I3=Is[2]
	CharInvs(Sig, Is);
	if (fabs(Is[1])<ZERO)
		throw new Fatal(_("TijTensor::Derivs: (I2>0.0) failed (I2=%.6f). Sig=<%.6f %.6f %.6f %.6f %.6f %.6f>"), Is[1],
				Sig(0), Sig(1), Sig(2), Sig(3)/sqrt(2.0), Sig(4)/sqrt(2.0), Sig(5)/sqrt(2.0));

	// Square root of I3/I2
	sqrt_rz = Is[2]/Is[1];
	if (fabs(sqrt_rz)<ZERO)
		throw new Fatal(_("TijTensor::Derivs: (I3/I2>0) failed (I3/I2=%.6f). Sig=<%.6f %.6f %.6f %.6f %.6f %.6f>"), sqrt_rz,
				Sig(0), Sig(1), Sig(2), Sig(3)/sqrt(2.0), Sig(4)/sqrt(2.0), Sig(5)/sqrt(2.0));
	sqrt_rz = sqrt(sqrt_rz); // sqrt_rz = sqrt(I3/I2)

	// Professr Nakai's stress invariants
	Stress_tn_ts(Sig, tn,ts);

	// Tau = square root of Sig
	Sqrt(Sig, Tau);

	// inverse of Tau
	Inv(Tau, inv_Tau);

	// a tensor
	a = inv_Tau * sqrt_rz;

	// tij tensor
	Mult(a,Sig, t);
	
	// -------------------------------------------------------------- derivatives calculation

	Tensor2 SigSig;
	Mult(Sig,Sig, SigSig);
	
	Tensor2 dI1_ds; dI1_ds = I*Is[0] - Sig;
	Tensor2 dI2_ds; dI2_ds = SigSig  - Sig*Is[0] + I*Is[1];
	
	dtn_ds = dI2_ds*(3.0/Is[1]) - dI1_ds*(tn/Is[1]);

	// Check ts
	if (fabs(ts)<_ts_min)
	{
		dts_ds = 0.0; // In fact, it is undetermined
		dts_dt = 0.0; // In fact, it is undetermined
	}
	else
	{
		double alp = tn/(6.0*ts);
		double bet = tn*(6.0*tn-Is[0])/(6.0*ts*Is[1]);
		double gam = (Is[0]-6.0*tn)/(2.0*ts*Is[1]);
		    dts_ds = I*alp + dI1_ds*bet + dI2_ds*gam;
		    dts_dt = t*(1.0/ts) - a*(tn/ts);
	}
		
}


}; // namespace Tensors

#endif // MECHSYS_TIJTENSOR_H
