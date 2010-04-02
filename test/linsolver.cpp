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
#include <iomanip>

// MechSys
#ifdef HAVE_UMFPACK
  #include <mechsys/linalg/umfpack.h>
#endif
#ifdef HAVE_SUPERLU
  #include <mechsys/linalg/superlu.h>
#endif
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/sparse_matrix.h>
#include <mechsys/linalg/sparse_triplet.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;

void Dense2Sparse(Mat_t const & D, Sparse::Matrix<double,int> & S, double Zero=1.0e-9)
{
	// Size
	int m = D.num_rows();
	int n = D.num_cols();

	// Number of nonzero elements
	int nz = 0;
	for (int i=0; i<m; ++i)
	for (int j=0; j<n; ++j)
		if (fabs(D(i,j))>Zero) nz++;

	// Sparse triplet
	Sparse::Triplet<double,int> T;
	T.AllocSpace(/*nRows*/m, /*nCols*/n, /*Size*/nz);

	// Fill arrays
	for (int i=0; i<m; ++i)
	for (int j=0; j<n; ++j)
		if (fabs(D(i,j))>Zero) T.PushEntry(i,j,D(i,j));

	// Sparse matrix
	S.Set(T);
}

int SolveAndOutput(int SolverType, Mat_t & A, Vec_t & Y, Vec_t const & Ycorr)
{
	// SolverType:
	//             1: => Dense
	//             2: => Sparse with UMFPACK
	//             3: => Sparse with SuperLU
	
	// Tolerance
	const double MAXERR = 1.0e-10;

	// Size
	int sz = A.num_rows();

	// Solve:  Y <- inv(A)*Y
	if (SolverType==1)
	{
		cout << "[1;33m... dense solver ...[0m" << endl;
        Sol (A,Y);
	}
	else if (SolverType==2)
	{
#ifdef HAVE_UMFPACK
		cout << "[1;33m... UMFPACK ...[0m" << endl;
		Vec_t             B(Y);
		Sparse::Matrix<double,int> S;
		Dense2Sparse   (A, S);
		UMFPACK::Solve (S,B, Y); // Y <- inv(S)*B
#else
		throw new Fatal("UMFPACK is not available");
#endif
	}
	else if (SolverType==3)
	{
#ifdef HAVE_SUPERLU
		cout << "[1;33m... SuperLU ...[0m" << endl;
		Sparse::Matrix<double,int> S;
		Dense2Sparse(A,S);
		//SuperLU::Solve(S, Y); // Y <- inv(S)*Y
#else
		throw new Fatal("SuperLU is not available");
#endif
	}
	else throw new Fatal("SolverType (%d) is invalid",SolverType);

	// Output
	cout << sz << " x " << sz << endl;
	cout << "Y     = " << Y;
	cout << "Ycorr = " << Ycorr << endl; // Ycorrect

	// Calculate error
	double norm_corr = 0.0;
	double     error = 0.0;
	for (int i=0; i<sz; ++i)
	{
		norm_corr += Ycorr(i)*Ycorr(i);
			error += (Y(i)-Ycorr(i))*(Y(i)-Ycorr(i));
	}
	error = sqrt(error)/sqrt(norm_corr);

	// Check
	bool ok = false;
	if (error>MAXERR) { ok=false; cout << "error="<<error<< "   .. FAILED .... FAILED .... FAILED .... FAILED .... FAILED .... FAILED .." << endl; }
	else              { ok=true;  cout << "error="<<error<< "   ..ok.." << endl; }
	cout << "\n\n" << endl;

	// Return a count to sum the number of errors
	return (ok ? 0 : 1);

}

int main(int argc, char **argv) try
{
	// SolverType:
	//             1: => Dense
	//             2: => Sparse with UMFPACK
	//             3: => Sparse with SuperLU
	
	// Input
	int solvertype = 1;
	if (argc>1) solvertype = atoi(argv[1]);
	if (!(solvertype==1 || solvertype==2 || solvertype==3))
		throw new Fatal("Please, call this program as follows: %s SolverType, where:\nSolverType:\n \t1: => Dense\n \t2: => Sparse with UMFPACK\n \t3: => Sparse with SuperLU\n",argv[0]);

	int num_errors=0;

	{
		Mat_t A(9,9); A = 19,   3,   1,   12,   1,   16,   1,   3,  11,
		                 -19,   3,   1,   12,   1,   16,   1,   3,  11,
		                 -19,  -3,   1,   12,   1,   16,   1,   3,  11,
		                 -19,  -3,  -1,   12,   1,   16,   1,   3,  11,
		                 -19,  -3,  -1,  -12,   1,   16,   1,   3,  11,
		                 -19,  -3,  -1,  -12,  -1,   16,   1,   3,  11,
		                 -19,  -3,  -1,  -12,  -1,  -16,   1,   3,  11,
		                 -19,   3,  -1,  -12,  -1,  -16,  -1,   3,  11,
		                 -19,  -3,  -1,  -12,  -1,  -16,  -1,  -3,  11;
		Vec_t Y    (9); Y = 0, 0, 1,  0, 0, 0,  0, 0, 0;
		Vec_t Ycorr(9); Ycorr = 0,  -0.16666666666666666, 0.5, 0,  0,  0,  -0.5, 0.16666666666666666, 0;
		//num_errors += SolveAndOutput(solvertype, A,Y,Ycorr);
	}

	{
		Mat_t A(5,5); A = 2,  3,  0,  0,  0,
		                  3,  0,  4,  0,  6,
		                  0, -1, -3,  2,  0,
		                  0,  0,  1,  0,  0,
		                  0,  4,  2,  0,  1;
		Vec_t Y    (5); Y = 8, 45, -3, 3, 19;
		Vec_t Ycorr(5); Ycorr = 1, 2, 3, 4, 5;
		num_errors += SolveAndOutput(solvertype, A,Y,Ycorr);
	}

	{
		Mat_t A(5,5); A = 5,  4,  3,  2,  1,
		                  4,  4,  3,  2,  1,
		                  0,  3,  3,  2,  1,
		                  0,  0,  2,  2,  1,
		                  0,  0,  0,  1,  1;
		Vec_t Y    (5); Y = 1, 2, 3, 4, 5;
		Vec_t Ycorr(5); Ycorr = -1, 3, -10,  19, -14;
		num_errors += SolveAndOutput(solvertype, A,Y,Ycorr);
	}

	{
		Mat_t A(6,6); A = 6,  5,  4,  3,  2,  1,
		                  5,  5,  4,  3,  2,  1,
		                  0,  4,  4,  3,  2,  1,
		                  0,  0,  3,  3,  2,  1,
		                  0,  0,  0,  2,  2,  1,
		                  0,  0,  0,  0,  1,  1;
		Vec_t Y    (6); Y = 1, 2, 3, 4, 5, 6;
		Vec_t Ycorr(6); Ycorr = -1,  4,  -17,  50, -101,  107;
		num_errors += SolveAndOutput(solvertype, A,Y,Ycorr);
	}

	{
		Mat_t A(7,7); A = 7,  6,  5,  4,  3,  2,  1,
		                  6,  6,  5,  4,  3,  2,  1,
		                  0,  5,  5,  4,  3,  2,  1,
		                  0,  0,  4,  4,  3,  2,  1,
		                  0,  0,  0,  3,  3,  2,  1,
		                  0,  0,  0,  0,  2,  2,  1,
		                  0,  0,  0,  0,  0,  1,  1;
		Vec_t Y    (7); Y = 1, 2, 3, 4, 5, 6, 7;
		Vec_t Ycorr(7); Ycorr = -1,  5,  -26,  103, -310,  619, -612;
		num_errors += SolveAndOutput(solvertype, A,Y,Ycorr);
	}

	{
		Mat_t A(10,10); A = 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     9,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,  8,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,  0,  7,  7,  6,  5,  4,  3,  2,  1,
		                     0,  0,  0,  6,  6,  5,  4,  3,  2,  1,
		                     0,  0,  0,  0,  5,  5,  4,  3,  2,  1,
		                     0,  0,  0,  0,  0,  4,  4,  3,  2,  1,
		                     0,  0,  0,  0,  0,  0,  3,  3,  2,  1,
		                     0,  0,  0,  0,  0,  0,  0,  2,  2,  1,
		                     0,  0,  0,  0,  0,  0,  0,  0,  1,  1;
		Vec_t Y    (10); Y = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
		Vec_t Ycorr(10); Ycorr = -1,  8,  -65,  454,  -2725,  13624,  -54497,  163490, -326981, 326991;
		num_errors += SolveAndOutput(solvertype, A,Y,Ycorr);
	}

	{
		Mat_t A(20,20); A = 0                  ,  0                  ,  0                 ,  3.5675854203801425,  5.586999952339394  ,  0                 ,  9.124672011124929  ,  9.145908826954795 ,  1.2416508676904114,  6.38300935148842  ,  8.372597363621663 ,  0                  ,  4.8937055932552855 ,  2.122155506388934 ,  1.9249910701053075,  8.10075353261826  ,  0                  ,  0                 ,  0                 ,  0                 ,
		                    0                  ,  8.034018687286288  ,  9.138640099458767 ,  0                 ,  7.167512496608924  ,  0                 ,  0                  ,  0.5224156775260869,  9.151075423706791 ,  0.5527127803943133,  0                 ,  0                  ,  0                  ,  0                 ,  2.272425659398155 ,  0                 ,  0.09821762809327228,  0                 ,  1.0346169405287153,  0                 ,
		                    3.933323734304129  ,  0                  ,  6.045976642762936 ,  0.6280320819039942,  0                  ,  0                 ,  0.6140203425194934 ,  6.39310615769223  ,  0                 ,  4.037854158235598 ,  1.9017865519581822,  3.602733303541066  ,  0                  ,  2.686808823172111 ,  5.9502664007032156,  0                 ,  4.002550260487107  ,  4.627875905089286 ,  1.4056849404831984,  0                 ,
		                    0                  ,  9.660110167917344  ,  0                 ,  5.075437326736201 ,  8.229344625037534  ,  0                 ,  8.583627552501572  ,  1.7807198014928283,  0                 ,  1.6191894089316983,  3.180176627991914 ,  0.5028035436711342 ,  0                  ,  0                 ,  0                 ,  0                 ,  0                  ,  1.1694508123497538,  0                 ,  4.973393581488593 ,
		                    0                  ,  0                  ,  0                 ,  0                 ,  8.955387440376533  ,  0                 ,  0                  ,  0                 ,  0                 ,  6.737033231793585 ,  0                 ,  0                  ,  0.8406995803220718 ,  0                 ,  0                 ,  4.3216109960992855,  0                  ,  3.4671405554185952,  0                 ,  0                 ,
		                    1.385736985109698  ,  0                  ,  0                 ,  0                 ,  0                  ,  0                 ,  5.756620211148549  ,  8.43296076919799  ,  5.689564996693641 ,  0                 ,  0                 ,  0                  ,  0                  ,  0                 ,  0                 ,  0                 ,  0                  ,  0                 ,  0                 ,  3.4900697970019823,
		                    2.0789969292331802 ,  5.974867554471115  ,  4.338542348294354 ,  8.046348353368327 ,  0.08751824121803531,  0                 ,  0.43480370689359615,  0                 ,  0                 ,  0                 ,  0                 ,  0.21197661565830694,  5.878441983814603  ,  0                 ,  0                 ,  0                 ,  0                  ,  6.788854452000153 ,  6.062662267233091 ,  6.728162823511978 ,
		                    0.5374390217351777 ,  8.771369977500436  ,  0                 ,  1.8018373020167933,  0                  ,  5.246951774894724 ,  0                  ,  0                 ,  4.200189050235592 ,  0                 ,  5.824521414137147 ,  7.705584595537237  ,  8.316324494690875  ,  0                 ,  0                 ,  0                 ,  0                  ,  0                 ,  1.277254092094513 ,  7.574874413062603 ,
		                    0                  ,  3.9584190446893897 ,  0                 ,  5.013657729990476 ,  0                  ,  6.385530108566956 ,  0                  ,  0                 ,  0                 ,  0                 ,  8.277569130450104 ,  6.899156203121713  ,  3.0655986254352285 ,  3.302314983284859 ,  0                 ,  5.656808001497019 ,  0                  ,  3.562520664465325 ,  0.920620182277091 ,  0                 ,
		                    2.9652891618821853 ,  9.605098735961636  ,  5.969015261347524 ,  8.486840233026586 ,  6.570461155812925  ,  2.7689233965970916,  4.891827503661292  ,  6.103225568342411 ,  0                 ,  1.192249182768621 ,  9.755511421595996 ,  0                  ,  1.315292257739018  ,  0                 ,  0                 ,  5.035028083029985 ,  0                  ,  0                 ,  7.71746739264559  ,  0                 ,
		                    3.236407921267915  ,  0                  ,  0                 ,  0                 ,  0.1393753667108022 ,  7.857758401045775 ,  0.11713860216390648,  1.4631471861831646,  8.147193514176983 ,  0                 ,  0                 ,  8.32428646975454   ,  0                  ,  0                 ,  0                 ,  9.360562742016489 ,  0                  ,  0                 ,  3.0181093147037608,  0.1813351619034509,
		                    0                  ,  0                  ,  7.7372331399768415,  0                 ,  4.300505563006509  ,  0                 ,  0                  ,  0                 ,  0                 ,  8.967138515878984 ,  0                 ,  0                  ,  0                  ,  0                 ,  5.516559493051853 ,  0                 ,  0                  ,  0                 ,  0                 ,  7.5308398601016755,
		                    0                  ,  0                  ,  0                 ,  0                 ,  0.7260314156986958 ,  9.641906392874098 ,  8.647227787313017  ,  0                 ,  4.843826767243984 ,  7.162499107286214 ,  0.7161359833354453,  0                  ,  5.055131381762084  ,  0.8666467595989003,  0                 ,  1.0712736864402794,  0                  ,  4.6832902662542   ,  6.741092859798396 ,  0                 ,
		                    6.723683426451984  ,  1.7347989745713865 ,  2.81269142489325  ,  3.161818033414278 ,  3.789589775557017  ,  0.2973896849528046,  9.636707282777255  ,  0                 ,  7.243793402664997 ,  0                 ,  9.141605456332742 ,  0                  ,  0                  ,  1.7563111778436968,  0                 ,  0                 ,  0                  ,  0                 ,  0                 ,  9.107769992791393 ,
		                    7.660177348481911  ,  9.757878945492967  ,  0                 ,  0                 ,  5.454970391258724  ,  0                 ,  1.480182691130142  ,  0                 ,  3.980512003770281 ,  5.912202052303711 ,  5.275390362376423 ,  9.402134883591792  ,  1.5325609264828544 ,  0                 ,  0                 ,  0.9688270177141955,  8.961998004688787  ,  9.682075344044918 ,  2.3331578962435398,  0                 ,
		                    0                  ,  0.7420651996198635 ,  5.499413067550481 ,  0                 ,  0                  ,  0                 ,  8.985616883119826  ,  0                 ,  0                 ,  0                 ,  1.632453790716828 ,  6.233979849609247  ,  0                  ,  6.179287061888416 ,  2.8706643104880083,  0                 ,  0                  ,  0                 ,  0                 ,  0                 ,
		                    0                  ,  0                  ,  8.522563946980364 ,  6.962833350821894 ,  6.390638362572334  ,  7.7503519153029155,  6.6588096153608    ,  0                 ,  3.148540643401213 ,  0                 ,  0                 ,  0                  ,  2.7875828171548225 ,  0                 ,  0                 ,  0                 ,  6.355005714351777  ,  4.733037768787803 ,  0.683004053735865 ,  0                 ,
		                    0                  ,  2.654866147252987  ,  0                 ,  0.2777425437946157,  0                  ,  4.871086319563457 ,  0                  ,  0                 ,  0.2486323604990659,  6.157788721206641 ,  9.262154654959115 ,  0                  ,  0.02717699005374219,  0                 ,  0                 ,  6.033192161053535 ,  0                  ,  4.740589637450165 ,  6.601074976942364 ,  0.6917756252120999,
		                    0.22916959815091897,  3.9489130728041766 ,  0                 ,  3.816313079140624 ,  7.486380981073855  ,  0                 ,  2.5056964723817443 ,  0                 ,  4.839742316423194 ,  0.309707325646863 ,  0                 ,  1.3273283407419922 ,  0                  ,  7.707075561098562 ,  8.049161871827437 ,  0                 ,  7.586628280796931  ,  3.544465739650361 ,  3.133452745731952 ,  9.449910658959729 ,
		                    0                  ,  0.44951693816252636,  0                 ,  2.696042434612286 ,  9.721936217783874  ,  0                 ,  4.732090227606945  ,  9.77128819033402  ,  0                 ,  0                 ,  3.7389230532001907,  0                  ,  7.470756543006561  ,  0                 ,  0                 ,  0.3021467296140301,  4.713176433137849  ,  7.67146049488812  ,  0                 ,  3.053629908626001 ;
		Vec_t Y    (20); Y = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;
		Vec_t Ycorr(20); Ycorr = -5.83298035975012485466, -5.02528090315008135747, -0.65339477675371748777,
		                          3.28342871409988967812,  2.17233316234833395697, -3.02582636229451429344,
		                          0.49178772642741530596, -0.61316205636004328383,  2.81474223668284828648,
		                          0.26549722703766331922,  3.03791749623312412609,  6.32611698147693246597,
		                         -2.48386676081233570557, -4.50044353516062134446,  0.80704125950832916736,
		                         -3.97526809879191134200, -0.57278944457015512626,  0.87247668887712592767,
		                          5.23101362560551041980,  0.11691861137762919742;
		num_errors += SolveAndOutput(solvertype, A,Y,Ycorr);
	}
	
	{
		Mat_t A(30,30); A = 30,  29,  28,  27,  26,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                    29,  29,  28,  27,  26,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,  28,  28,  27,  26,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,  27,  27,  26,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,  26,  26,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,  25,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,  24,  24,  23,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,  23,  23,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,  22,  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,  21,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,  20,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  19,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  18,  18,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  17,  17,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  16,  16,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  15,  15,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  14,  14,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  13,  13,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  12,  12,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  11,  11,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  10,  10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  8,  8,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  7,  7,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  6,  6,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  5,  5,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  4,  4,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  3,  3,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,  2,  2,  1,
		                     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,  0,  1,  1;
		Vec_t Y(30); Y = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30;
		Vec_t Ycorr(30);
		Ycorr =  0.033254755522863,  -0.068694428519376,  -0.035873196189045,
		        -0.036919475284371,  -0.038407556802543,  -0.039933381619058,
		        -0.041591631634103,  -0.043392870280541,  -0.045356832247099,
		        -0.047506523988869,  -0.049869520156956,  -0.052479117021596,
		        -0.055375893610983,  -0.058609808613377,  -0.062243062185938,
		        -0.066354067210954,  -0.071043059046661,  -0.076440232393419,
		        -0.082717211278872,  -0.090110675932401,  -0.098893240676007,
		        -0.109960833915966,  -0.120313328672264,  -0.157806699294130,
		        -0.053159804235202,  -0.734200978823981,   1.936803915295993,
		        -6.810411745888430,  12.620823491777838,  16.379176508221523;
		//num_errors += SolveAndOutput(solvertype, A,Y,Ycorr);
	}

	cout << "Number of errors = " << num_errors << endl;
	return (num_errors>0? 1 : 0);
}
MECHSYS_CATCH