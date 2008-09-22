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

#ifndef MECHSYS_FEM_QUADRATURE_H
#define MECHSYS_FEM_QUADRATURE_H

// MechSys
#include "util/array.h"
#include "util/exception.h"

namespace FEM
{

/** Auxiliar structure (integration point). */
struct IntegPoint
{
	double r; ///< xsi coordinate
	double s; ///< eta coordinate
	double t; ///< zeta coordinate
	double w; ///< weight (coeficient) W 
};

// Weights after Akin -- Rice University -- http://www.owlnet.rice.edu/~mech517/

// NumGP = 2  =>  precision = 3
const double GAUSSP_2[] = { -0.577350269189625764509149e+00,
                             0.577350269189625764509149e+00 };
const double WEIGHT_2[] = {  1.000000000000000000000000e+00,
                             1.000000000000000000000000e+00 };

// NumGP = 3  =>  precision = 5
const double GAUSSP_3[] = { -0.774596669241483377035835e+00,
                             0.000000000000000000000000e+00,
                             0.774596669241483377035835e+00 };
const double WEIGHT_3[] = {  0.555555555555555555555556e+00,
                             0.888888888888888888888889e+00,
                             0.555555555555555555555556e+00 };

// NumGP = 4  =>  precision = 7
const double GAUSSP_4[] = { -0.861136311594052575223946e+00,
                            -0.339981043584856264802666e+00,
                             0.339981043584856264802666e+00,
                             0.861136311594052575223946e+00 };
const double WEIGHT_4[] = {  0.34785484513745385737306e+00,
                             0.65214515486254614262694e+00,
                             0.65214515486254614262694e+00,
                             0.34785484513745385737306e+00 };

// NumGP = 5  =>  precision = 9
const double GAUSSP_5[] = { -0.906179845938663992797627e+00,
                            -0.538469310105683091036314e+00,
                             0.000000000000000000000000e+00,
                             0.538469310105683091036314e+00,
                             0.906179845938663992797627e+00 };
const double WEIGHT_5[] = {  0.23692688505618908751426e+00,
                             0.47862867049936646804129e+00,
                             0.56888888888888888888889e+00,
                             0.47862867049936646804129e+00,
                             0.23692688505618908751426e+00 };

// NumGP = 6  =>  precision = 11
const double GAUSSP_6[] = { -0.932469514203152027812302e+00,
                            -0.661209386466264513661400e+00,
                            -0.238619186083196908630502e+00,
                             0.238619186083196908630502e+00,
                             0.661209386466264513661400e+00,
                             0.932469514203152027812302e+00 };
const double WEIGHT_6[] = {  0.17132449237917034504030e+00,
                             0.36076157304813860756983e+00,
                             0.46791393457269104738987e+00,
                             0.46791393457269104738987e+00,
                             0.36076157304813860756983e+00,
                             0.17132449237917034504030e+00 };

// NumGP = 7  =>  precision = 13
const double GAUSSP_7[] = { -0.949107912342758524526190e+00,
                            -0.741531185599394439863865e+00,
                            -0.405845151377397166906607e+00,
                             0.000000000000000000000000e+00,
                             0.405845151377397166906607e+00,
                             0.741531185599394439863865e+00,
                             0.949107912342758524526190e+00 };
const double WEIGHT_7[] = {  0.12948496616886969327061e+00,
                             0.27970539148927666790147e+00,
                             0.38183005050511894495037e+00,
                             0.41795918367346938775510e+00,
                             0.38183005050511894495037e+00,
                             0.27970539148927666790147e+00,
                             0.12948496616886969327061e+00 };

// NumGP = 8  =>  precision = 15
const double GAUSSP_8[] = { -0.960289856497536231683561e+00,
                            -0.796666477413626739591554e+00,
                            -0.525532409916328985817739e+00,
                            -0.183434642495649804939476e+00,
                             0.183434642495649804939476e+00,
                             0.525532409916328985817739e+00,
                             0.796666477413626739591554e+00,
                             0.960289856497536231683561e+00 };
const double WEIGHT_8[] = {  0.10122853629037625915253e+00,
                             0.22238103445337447054436e+00,
                             0.31370664587788728733796e+00,
                             0.36268378337836198296515e+00,
                             0.36268378337836198296515e+00,
                             0.31370664587788728733796e+00,
                             0.22238103445337447054436e+00,
                             0.10122853629037625915253e+00 };

// NumGP = 9  =>  precision = 17
const double GAUSSP_9[] = { -0.968160239507626089835576e+00,
                            -0.836031107326635794299430e+00,
                            -0.613371432700590397308702e+00,
                            -0.324253423403808929038538e+00,
                             0.000000000000000000000000e+00,
                             0.324253423403808929038538e+00,
                             0.613371432700590397308702e+00,
                             0.836031107326635794299430e+00,
                             0.968160239507626089835576e+00 };
const double WEIGHT_9[] = {  0.081274388361574411971890e+00,
                             0.18064816069485740405847e+00,
                             0.26061069640293546231874e+00,
                             0.31234707704000284006863e+00,
                             0.33023935500125976316453e+00,
                             0.31234707704000284006863e+00,
                             0.26061069640293546231874e+00,
                             0.18064816069485740405847e+00,
                             0.08127438836157441197189e+00 };

// NumGP = 10  =>  precision = 19
const double GAUSSP_10[] = { -0.973906528517171720077964e+00,
                             -0.865063366688984510732097e+00,
                             -0.679409568299024406234327e+00,
                             -0.433395394129247190799266e+00,
                             -0.148874338981631210884826e+00,
                              0.148874338981631210884826e+00,
                              0.433395394129247190799266e+00,
                              0.865063366688984510732097e+00,
                              0.679409568299024406234327e+00,
                              0.973906528517171720077964e+00 };
const double WEIGHT_10[] = {  0.066671344308688137593570e+00,
                              0.14945134915058059314578e+00,
                              0.21908636251598204399554e+00,
                              0.26926671930999635509123e+00,
                              0.29552422471475287017389e+00,
                              0.29552422471475287017389e+00,
                              0.26926671930999635509123e+00,
                              0.14945134915058059314578e+00,
                              0.21908636251598204399554e+00,
                              0.06667134430868813759357e+00 };

inline void SetGaussIP (size_t NDim, size_t NumGP1D, Array<FEM::IntegPoint> & IPs, bool IsTriOrTet=false)
{
	if (IsTriOrTet)
	{
		if (NDim==2) // triangle
		{
			switch (NumGP1D) // NumGP1D==NumGPTotal for triangle/tetrahedron
			{
				case 3:
				{
					IPs.Resize(3);
					IPs[0].r = 1.0/6.0;   IPs[0].s = 1.0/6.0;   IPs[0].t = 0.0;   IPs[0].w = 1.0/6.0;
					IPs[1].r = 2.0/3.0;   IPs[1].s = 1.0/6.0;   IPs[1].t = 0.0;   IPs[1].w = 1.0/6.0;
					IPs[2].r = 1.0/6.0;   IPs[2].s = 2.0/3.0;   IPs[2].t = 0.0;   IPs[2].w = 1.0/6.0;
					return;
				}
				case 7:
				{
					IPs.Resize(7);
					IPs[0].r = 0.1012865073235;   IPs[0].s = 0.1012865073235;   IPs[0].t = 0.0;   IPs[0].w = 0.0629695902724;
					IPs[1].r = 0.7974269853531;   IPs[1].s = 0.1012865073235;   IPs[1].t = 0.0;   IPs[1].w = 0.0629695902724;
					IPs[2].r = 0.1012865073235;   IPs[2].s = 0.7974269853531;   IPs[2].t = 0.0;   IPs[2].w = 0.0629695902724;
					IPs[3].r = 0.4701420641051;   IPs[3].s = 0.0597158717898;   IPs[3].t = 0.0;   IPs[3].w = 0.0661970763942;
					IPs[4].r = 0.4701420641051;   IPs[4].s = 0.4701420641051;   IPs[4].t = 0.0;   IPs[4].w = 0.0661970763942;
					IPs[5].r = 0.0597158717898;   IPs[5].s = 0.4701420641051;   IPs[5].t = 0.0;   IPs[5].w = 0.0661970763942;
					IPs[6].r = 0.3333333333333;   IPs[6].s = 0.3333333333333;   IPs[6].t = 0.0;   IPs[6].w = 0.1125000000000;
					return;
				}
				case 13:
				{
					IPs.Resize(13);
					IPs[ 0].r = 0.0651301029022;   IPs[ 0].s = 0.0651301029022;   IPs[ 0].t = 0.0;   IPs[ 0].w =  0.0266736178044;
					IPs[ 1].r = 0.8697397941956;   IPs[ 1].s = 0.0651301029022;   IPs[ 1].t = 0.0;   IPs[ 1].w =  0.0266736178044;
					IPs[ 2].r = 0.0651301029022;   IPs[ 2].s = 0.8697397941956;   IPs[ 2].t = 0.0;   IPs[ 2].w =  0.0266736178044;
					IPs[ 3].r = 0.3128654960049;   IPs[ 3].s = 0.0486903154253;   IPs[ 3].t = 0.0;   IPs[ 3].w =  0.0385568804452;
					IPs[ 4].r = 0.6384441885698;   IPs[ 4].s = 0.3128654960049;   IPs[ 4].t = 0.0;   IPs[ 4].w =  0.0385568804452;
					IPs[ 5].r = 0.0486903154253;   IPs[ 5].s = 0.6384441885698;   IPs[ 5].t = 0.0;   IPs[ 5].w =  0.0385568804452;
					IPs[ 6].r = 0.6384441885698;   IPs[ 6].s = 0.0486903154253;   IPs[ 6].t = 0.0;   IPs[ 6].w =  0.0385568804452;
					IPs[ 7].r = 0.3128654960049;   IPs[ 7].s = 0.6384441885698;   IPs[ 7].t = 0.0;   IPs[ 7].w =  0.0385568804452;
					IPs[ 8].r = 0.0486903154253;   IPs[ 8].s = 0.0486903154253;   IPs[ 8].t = 0.0;   IPs[ 8].w =  0.0385568804452;
					IPs[ 9].r = 0.2603459660790;   IPs[ 9].s = 0.2603459660790;   IPs[ 9].t = 0.0;   IPs[ 9].w =  0.0878076287166;
					IPs[10].r = 0.4793080678419;   IPs[10].s = 0.2603459660790;   IPs[10].t = 0.0;   IPs[10].w =  0.0878076287166;
					IPs[11].r = 0.2603459660790;   IPs[11].s = 0.4793080678419;   IPs[11].t = 0.0;   IPs[11].w =  0.0878076287166;
					IPs[12].r = 0.3333333333333;   IPs[12].s = 0.3333333333333;   IPs[12].t = 0.0;   IPs[12].w = -0.0747850222338;
					return;
				}
				default: throw new Fatal("FEM::SetGaussIP: Total number of Gauss-Points for triangle (NumGPTotal==%d) is invalid.",NumGP1D);
			}
		}
		else if (NDim==3) // tetrahedron
		{
		}
		else throw new Fatal("FEM::SetGaussIP: NDim==%d for triangle/tetrahedron is invalid.",NDim);
	}
	else
	{
		double const * gps = NULL;
		double const * wgt = NULL;
		switch (NumGP1D)
		{
			case  2: { gps=GAUSSP_2;   wgt=WEIGHT_2;   break; }
			case  3: { gps=GAUSSP_3;   wgt=WEIGHT_3;   break; }
			case  4: { gps=GAUSSP_4;   wgt=WEIGHT_4;   break; }
			case  5: { gps=GAUSSP_5;   wgt=WEIGHT_5;   break; }
			case  6: { gps=GAUSSP_6;   wgt=WEIGHT_6;   break; }
			case  7: { gps=GAUSSP_7;   wgt=WEIGHT_7;   break; }
			case  8: { gps=GAUSSP_8;   wgt=WEIGHT_8;   break; }
			case  9: { gps=GAUSSP_9;   wgt=WEIGHT_9;   break; }
			case 10: { gps=GAUSSP_10;  wgt=WEIGHT_10;  break; }
			default: throw new Fatal("FEM::SetGaussIP: Number of Gauss-Points 1D (NumGP1D==%d) is invalid.",NumGP1D);
		}
		if (NDim==1)
		{
			IPs.Resize (NumGP1D);
			for (size_t i=0; i<NumGP1D; ++i)
			{
				IPs[i].r = gps[i];
				IPs[i].s = 0.0;
				IPs[i].t = 0.0;
				IPs[i].w = wgt[i];
			}
		}
		else if (NDim==2)
		{
			size_t m = 0;
			IPs.Resize (NumGP1D*NumGP1D);
			for (size_t i=0; i<NumGP1D; ++i)
			for (size_t j=0; j<NumGP1D; ++j)
			{
				IPs[m].r = gps[j];
				IPs[m].s = gps[i];
				IPs[m].t = 0.0;
				IPs[m].w = wgt[j]*wgt[i];
				m++;
			}
		}
		else if (NDim==3)
		{
			size_t m = 0;
			IPs.Resize (NumGP1D*NumGP1D*NumGP1D);
			for (size_t i=0; i<NumGP1D; ++i)
			for (size_t j=0; j<NumGP1D; ++j)
			for (size_t k=0; k<NumGP1D; ++k)
			{
				IPs[m].r = gps[k];
				IPs[m].s = gps[j];
				IPs[m].t = gps[i];
				IPs[m].w = wgt[k]*wgt[j]*wgt[i];
				m++;
			}
		}
		else throw new Fatal("FEM::SetGaussIP: NDim==%d is invalid.",NDim);
	}
}

}; // namespace FEM

#endif // MECHSYS_FEM_QUADRATURE_H
