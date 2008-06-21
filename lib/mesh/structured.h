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

#ifndef MECHSYS_MESH_STRUCTURED_H
#define MECHSYS_MESH_STRUCTURED_H

// STL
#include <iostream>
#include <fstream>

// MechSys
#include "util/array.h"
#include "util/exception.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/laexpr.h"

#ifndef VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

using LinAlg::Vector;
using LinAlg::Matrix;

namespace Mesh
{

struct Vertex
{
	long           MyID;   ///< ID
	bool           OnBry;  ///< On boundary?
	int            BryMrk; ///< Boundary mark
	int            Count;  ///< How many times this vertex was defined (repeated)
	Vector<double> C;      ///< X, Y, and Z coordinates
};

struct Elem
{
	long           MyID;  ///< ID
	bool           OnBry; ///< On boundary?
	Array<Vertex*> C;     ///< Connectivity
};

class Block
{
public:
	// Methods
	void Set (Matrix<double> const & C, Matrix<double> const & W, int nDivX, int nDivY, int nDivZ=0); ///< C=coordinates, W=weights
	/* 
	 * 2D: 8 nodes => C.Resize(2, 8)  :  W.Resize(2, max(nDivX,nDivY) +1 )
	 *         _                         _
	 *    C = |  x0 x1 x2 x3 x4 x5 x6 x7  |
	 *        |_ y0 y1 y2 y3 y4 y5 y6 y7 _|
	 *         _                                          _
	 *    W = |  wx0 wx1 wx2 wx3 ..... wx(nDivX-1)  NULL   |
	 *        |_ wy0 wy1 ... wy(nDivY-1) 0 0 0 0 0  NULL  _|
	 *
	 *
	 * 3D: 20 nodes => C.Resize(3, 20)  :  W.Resize(3, max(nDivX,nDivY,nDivZ) +1 )
	 *         _                                         _
	 *        |  x0 x1 x2 x3 x4 x5 x6 x7 ... x17 x18 x19  |
	 *    C = |  y0 y1 y2 y3 y4 y5 y6 y7 ... y17 y18 y19  |
	 *        |_ z0 z1 z2 z3 z4 z5 z6 z7 ... z17 z18 z19 _|
	 *         _                                          _
	 *        |  wx0 wx1 wx2 wx3 ..... wx(nDivX-1)  NULL   |
	 *    W = |  wy0 wy1 ... wy(nDivY-1) 0 0 0 0 0  NULL   |
	 *        |_ wz0 wz1 wz2 ... wz(nDivZ-1) 0 0 0  NULL  _|
	 */

	// Access methods
	int    nDivX      ()                   const { return _n_div_x;      }
	int    nDivY      ()                   const { return _n_div_y;      }
	int    nDivZ      ()                   const { return _n_div_z;      }
	double SumWeightX ()                   const { return _sum_weight_x; }
	double SumWeightY ()                   const { return _sum_weight_y; }
	double SumWeightZ ()                   const { return _sum_weight_z; }
	double W          (int iXYZ, int iDiv) const { return _W(iXYZ,iDiv); }

	// Access the coordinates of all 8 or 20 nodes
	Matrix<double> const & C() const { return _C; }

private:
	Matrix<double> _C;            ///< X, Y, and Z coordinates of all 8 or 20 nodes of this block
	Matrix<double> _W;            ///< X, Y, and Z weights
	int            _n_div_x;      ///< number of divisions along X
	int            _n_div_y;      ///< number of divisions along Y
	int            _n_div_z;      ///< number of divisions along Z
	double         _sum_weight_x; ///< sum of weights along X
	double         _sum_weight_y; ///< sum of weights along Y
	double         _sum_weight_z; ///< sum of weights along Z
}; // class Block

class Structured
{
public:
	// Constructor
	Structured (bool Is3D=false) : _is_3d(Is3D) { _s.Resize((Is3D?20:8)); }

	// Methods
	void Generate (Array<Block> const & Blocks);
	void WriteVTU (char const * FileName) const;

private:
	// Data
	bool           _is_3d;  ///< Is 3D mesh?
	Array<Vertex*> _verts;  ///< Vertices
	Array<Elem*>   _elems;  ///< Elements
	Vector<double> _s;      ///< Current shape (interpolation) values, computed just after _shape(r,s,t)

	// Private methods
	void _shape_2d (double r, double s);
	void _shape_3d (double r, double s, double t);
	void _vtk_con  (Elem const * E, String & Connect) const;

}; // class Structured


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Methods -- Block

inline void Block::Set(Matrix<double> const & C, Matrix<double> const & W, int nDivX, int nDivY, int nDivZ)
{
	_C            = C;
	_W            = W;
	_n_div_x      = nDivX;
	_n_div_y      = nDivY;
	_n_div_z      = nDivZ;
	_sum_weight_x = 0.0;
	_sum_weight_y = 0.0;
	_sum_weight_z = 0.0;
	for (int i=0; i<_n_div_x; ++i) _sum_weight_x += _W(0,i);
	for (int i=0; i<_n_div_y; ++i) _sum_weight_y += _W(1,i);
	for (int i=0; i<_n_div_z; ++i) _sum_weight_z += _W(2,i);
}

// Methods -- Structured

inline void Structured::Generate(Array<Block> const & Blocks)
{
	// Generate vertices and elements (with duplicates)
	for (size_t b=0; b<Blocks.Size(); ++b)
	{
		double t = -1.0; // initial Z natural coordinate
		for (int k=0; k<(_is_3d ? Blocks[b].nDivZ()+1 : 1); ++k)
		{
			double s = -1.0; // initial Y natural coordinate
			for (int j=0; j<Blocks[b].nDivY()+1; ++j)
			{
				double r = -1.0; // initial X natural coordinate
				for (int i=0; i<Blocks[b].nDivX()+1; ++i)
				{
					// Compute shape (interpolation) functions
					if (_is_3d) _shape_3d (r,s,t);
					else        _shape_2d (r,s);

					// New vertex
					Vertex * v = new Vertex;
					v->MyID = _verts.Size();
					v->C    = Blocks[b].C() * _s;
					if (_is_3d) v->OnBry = (i==0 || j==0 || k==0 || i==Blocks[b].nDivX() || j==Blocks[b].nDivY() || k==Blocks[b].nDivZ());
					else        v->OnBry = (i==0 || j==0 ||         i==Blocks[b].nDivX() || j==Blocks[b].nDivY());
					_verts.Push(v);

					// New element
					if (i!=0 && j!=0 && (_is_3d ? k!=0 : true))
					{
						Elem * e = new Elem;
						e->MyID = _elems.Size();
						if (_is_3d)
						{
							e->C.Resize(8);
							e->C[0] = _verts[v->MyID - 1 - (Blocks[b].nDivX()+1) - (Blocks[b].nDivX()+1)*(Blocks[b].nDivY()+1)];
							e->C[1] = _verts[v->MyID     - (Blocks[b].nDivX()+1) - (Blocks[b].nDivX()+1)*(Blocks[b].nDivY()+1)];
							e->C[2] = _verts[v->MyID                             - (Blocks[b].nDivX()+1)*(Blocks[b].nDivY()+1)];
							e->C[3] = _verts[v->MyID - 1                         - (Blocks[b].nDivX()+1)*(Blocks[b].nDivY()+1)];
							e->C[4] = _verts[v->MyID - 1 - (Blocks[b].nDivX()+1)];
							e->C[5] = _verts[v->MyID -     (Blocks[b].nDivX()+1)];
							e->C[6] = _verts[v->MyID];
							e->C[7] = _verts[v->MyID - 1];
							e->OnBry = (e->C[0]->OnBry || e->C[1]->OnBry || e->C[2]->OnBry || e->C[3]->OnBry || e->C[4]->OnBry || e->C[5]->OnBry || e->C[6]->OnBry || e->C[7]->OnBry);
						}
						else
						{
							e->C.Resize(4);
							e->C[0] = _verts[v->MyID - 1 - (Blocks[b].nDivX()+1)];
							e->C[1] = _verts[v->MyID     - (Blocks[b].nDivX()+1)];
							e->C[2] = _verts[v->MyID];
							e->C[3] = _verts[v->MyID - 1];
							e->OnBry = (e->C[0]->OnBry || e->C[1]->OnBry || e->C[2]->OnBry || e->C[3]->OnBry);
						}
						_elems.Push(e);
					}
					// Next r
					r += (2.0/Blocks[b].SumWeightX()) * Blocks[b].W(0,i);
				}
				// Next s
				s += (2.0/Blocks[b].SumWeightY()) * Blocks[b].W(1,j);
			}
			// Next t
			t += (_is_3d ? (2.0/Blocks[b].SumWeightZ()) * Blocks[b].W(2,k) : 0.0);
		}
	}

	// Remove duplicates
	for (size_t i=0; i<_verts.Size(); ++i)
	{
	}

	if (!_is_3d)
	{
	}
}

inline void Structured::WriteVTU(char const * FileName) const
{
	// Open File
	std::ofstream      of(FileName, std::ios::out);
	std::ostringstream oss;

	// Data
	size_t nn = _verts.Size(); // Number of Nodes
	size_t ne = _elems.Size(); // Number of Elements

	// Constants
	size_t          nimax = 40;        // number of integers in a line
	size_t          nfmax = 12;        // number of floats in a line
	Util::NumStream nsflo = Util::_6_3; // number format for floats

	// Header
	oss << "<?xml version=\"1.0\"?>\n";
	oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	oss << "  <UnstructuredGrid>\n";
	oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << ne << "\">\n";

	// Nodes: coordinates
	oss << "      <Points>\n";
	oss << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << "  " << nsflo <<         _verts[i]->C(0) << " ";
		oss <<         nsflo <<         _verts[i]->C(1) << " ";
		oss <<         nsflo << (_is_3d?_verts[i]->C(2):0.0);
		k++;
		VTU_NEWLINE (i,k,nn,nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Points>\n";

	// Elements: connectivity, offsets, types
	oss << "      <Cells>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		String con; _vtk_con (_elems[i], con);
		oss << "  " << con;
		k++;
		//VTU_NEWLINE (i,k,ne,nimax/(_is_3d?20:8),oss);
		VTU_NEWLINE (i,k,ne,nimax/(_is_3d?8:4),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t ossfset = 0;
	for (size_t i=0; i<ne; ++i)
	{
		//ossfset += (_is_3d?20:8);
		ossfset += (_is_3d?8:4);
		oss << (k==0?"  ":" ") << ossfset;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		//oss << (k==0?"  ":" ") << (_is_3d?25:23); // VTK_QUADRATIC_HEXAHEDRON or VTK_QUADRATIC_QUAD
		oss << (k==0?"  ":" ") << (_is_3d?12:9); // VTK_HEXAHEDRON or VTK_QUAD
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Cells>\n";

	// Data -- nodes
	oss << "      <PointData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "onbry" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << (k==0?"  ":" ") << _verts[i]->OnBry;
		k++;
		VTU_NEWLINE (i,k,nn,nfmax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </PointData>\n";

	// Data -- elements
	oss << "      <CellData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "onbry" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << _elems[i]->OnBry;
		k++;
		VTU_NEWLINE (i,k,ne,nfmax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </CellData>\n";

	// Bottom
	oss << "    </Piece>\n";
	oss << "  </UnstructuredGrid>\n";
	oss << "</VTKFile>" << std::endl;

	// Write to file
	of << oss.str();
	of.close();
}

// Private methods

inline void Structured::_shape_2d(double r, double s)
{
	/*      3           6            2
	 *        @---------@----------@
	 *        |               (1,1)|
	 *        |       s ^          |
	 *        |         |          |
	 *        |         |          |
	 *      7 @         +----> r   @ 5
	 *        |       (0,0)        |
	 *        |                    |
	 *        |                    |
	 *        |(-1,-1)             |
	 *        @---------@----------@
	 *      0           4            1
	 */
	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	_s(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
	_s(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
	_s(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
	_s(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
	_s(4) = 0.50*sm1*(1.0-r*r);
	_s(5) = 0.50*rp1*(1.0-s*s);
	_s(6) = 0.50*sp1*(1.0-r*r);
	_s(7) = 0.50*rm1*(1.0-s*s);
}

inline void Structured::_shape_3d(double r, double s, double t)
{
	/*                   t                         
	 *                   ^                         
	 *                   |                         
	 *                 4         15         7      
	 *                   o_______o________o        
	 *                 ,'|              ,'|        
	 *            12 o'  |         14 ,'  |        
	 *             ,'    |16        ,o    |19      
	 *       5   ,'      o      6 ,'      o        
	 *         o'=======o=======o'        |        
	 *         |      13 |      |         |        
	 *         |         |      |  11     |        
	 *      17 |       0 o______|_o_______o  ---> s
	 *         o       ,'       o       ,'  3      
	 *         |   8 o'      18 |     ,'           
	 *         |   ,'           |   ,o 10          
	 *         | ,'             | ,'               
	 *      1  o_______o________o'                 
	 *       ,'        9         2                 
	 *     ,'                                      
	 *    r                                        
	 */
	_s( 0) = 0.125*(1.0-r)  *(1.0-s)  *(1.0-t)  *(-r-s-t-2.0);
	_s( 1) = 0.125*(1.0+r)  *(1.0-s)  *(1.0-t)  *( r-s-t-2.0);
	_s( 2) = 0.125*(1.0+r)  *(1.0+s)  *(1.0-t)  *( r+s-t-2.0);
	_s( 3) = 0.125*(1.0-r)  *(1.0+s)  *(1.0-t)  *(-r+s-t-2.0);
	_s( 4) = 0.125*(1.0-r)  *(1.0-s)  *(1.0+t)  *(-r-s+t-2.0);
	_s( 5) = 0.125*(1.0+r)  *(1.0-s)  *(1.0+t)  *( r-s+t-2.0);
	_s( 6) = 0.125*(1.0+r)  *(1.0+s)  *(1.0+t)  *( r+s+t-2.0);
	_s( 7) = 0.125*(1.0-r)  *(1.0+s)  *(1.0+t)  *(-r+s+t-2.0);

	_s( 8) = 0.25 *(1.0-r*r)*(1.0-s)  *(1.0-t);
	_s( 9) = 0.25 *(1.0+r)  *(1.0-s*s)*(1.0-t);
	_s(10) = 0.25 *(1.0-r*r)*(1.0+s)  *(1.0-t);
	_s(11) = 0.25 *(1.0-r)  *(1.0-s*s)*(1.0-t);

	_s(12) = 0.25 *(1.0-r*r)*(1.0-s)  *(1.0+t);
	_s(13) = 0.25 *(1.0+r)  *(1.0-s*s)*(1.0+t);
	_s(14) = 0.25 *(1.0-r*r)*(1.0+s)  *(1.0+t);
	_s(15) = 0.25 *(1.0-r)  *(1.0-s*s)*(1.0+t);

	_s(16) = 0.25 *(1.0-r)  *(1.0-s)  *(1.0-t*t);
	_s(17) = 0.25 *(1.0+r)  *(1.0-s)  *(1.0-t*t);
	_s(18) = 0.25 *(1.0+r)  *(1.0+s)  *(1.0-t*t);
	_s(19) = 0.25 *(1.0-r)  *(1.0+s)  *(1.0-t*t);
}

inline void Structured::_vtk_con(Elem const * E, String & Connect) const
{
	if (_is_3d)
	{
		Connect.Printf("%d %d %d %d %d %d %d %d",E->C[1]->MyID,
		                                         E->C[2]->MyID,
		                                         E->C[3]->MyID,
		                                         E->C[0]->MyID,
		                                         E->C[5]->MyID,
		                                         E->C[6]->MyID,
		                                         E->C[7]->MyID,
		                                         E->C[4]->MyID);
	}
	else
	{
		Connect.Printf("%d %d %d %d",E->C[0]->MyID,
		                             E->C[1]->MyID,
		                             E->C[2]->MyID,
		                             E->C[3]->MyID);
	}
}


}; // namespace Mesh

#endif // MECHSYS_MESH_STRUCTURED_H
