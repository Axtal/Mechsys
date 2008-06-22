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
#include <cfloat>

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

	/* TODO:
	 *        1) Add additional checks, especially for 3D meshes
	 *        2) Remove boundary information for nodes duplicated not on boundary
	 *        3) Add boundary marks evaluation
	 *        4) Extend to second-order (o2) elements
	 *        5) Add quality check and improvement
	 */

struct Elem;

struct Share
{
	Elem* E; ///< The element
	int   N; ///< Local node index: 2D=>0,1,2,3, 3D=>0,1,2,3,4,5,6,7
};

struct Vertex
{
	long           MyID;   ///< ID
	bool           OnBry;  ///< On boundary?
	int            BryMrk; ///< Boundary mark
	bool           Dupl;   ///< Is this a duplicated node?
	Vector<double> C;      ///< X, Y, and Z coordinates
	Array<Share>   Shares; ///< Shared elements
};

struct Elem
{
	long           MyID;  ///< ID
	bool           OnBry; ///< On boundary?
	Array<Vertex*> V;     ///< Connectivity
};

class Block
{
public:
	// Constructor
	Block () : _n_div_x(0), _n_div_y(0), _n_div_z(0) {}

	// Methods
	void Set (Matrix<double> * C, Array<double> * Wx, Array<double> * Wy, Array<double> * Wz=NULL); ///< C=coordinates, W=weights
	/* 
	 * 2D: 8 nodes => C. Resize(2, 8)
	 *                Wx.Resize(nDivX)
	 *                Wy.Resize(nDivY)
	 *          _                         _
	 *    C  = |  x0 x1 x2 x3 x4 x5 x6 x7  |
	 *         |_ y0 y1 y2 y3 y4 y5 y6 y7 _|
	 *
	 *
	 * 3D: 20 nodes => C. Resize(3, 20)
	 *                 Wx.Resize(nDivX)
	 *                 Wy.Resize(nDivY)
	 *                 Wz.Resize(nDivZ)
	 *         _                                         _
	 *        |  x0 x1 x2 x3 x4 x5 x6 x7 ... x17 x18 x19  |
	 *    C = |  y0 y1 y2 y3 y4 y5 y6 y7 ... y17 y18 y19  |
	 *        |_ z0 z1 z2 z3 z4 z5 z6 z7 ... z17 z18 z19 _|
	 */

	// Access methods
	bool   Is3D       ()         const { return _is_3d;        }
	int    nDivX      ()         const { return _n_div_x;      }
	int    nDivY      ()         const { return _n_div_y;      }
	int    nDivZ      ()         const { return _n_div_z;      }
	double SumWeightX ()         const { return _sum_weight_x; }
	double SumWeightY ()         const { return _sum_weight_y; }
	double SumWeightZ ()         const { return _sum_weight_z; }
	double Wx         (int iDiv) const { return (*_wx)[iDiv];  }
	double Wy         (int iDiv) const { return (*_wy)[iDiv];  }
	double Wz         (int iDiv) const { return (*_wz)[iDiv];  }

	// Access the coordinates of all 8 or 20 nodes
	Matrix<double> const & C() const { return (*_c); }

private:
	Matrix<double> * _c;            ///< X, Y, and Z coordinates of all 8 or 20 nodes of this block
	Array <double> * _wx;           ///< X weights
	Array <double> * _wy;           ///< Y weights
	Array <double> * _wz;           ///< Z weights
	int              _n_div_x;      ///< number of divisions along X
	int              _n_div_y;      ///< number of divisions along Y
	int              _n_div_z;      ///< number of divisions along Z
	double           _sum_weight_x; ///< sum of weights along X
	double           _sum_weight_y; ///< sum of weights along Y
	double           _sum_weight_z; ///< sum of weights along Z
	bool             _is_3d;        ///< Is 3D block?
}; // class Block

class Structured
{
public:
	// Constructor
	Structured (double Tol=sqrt(DBL_EPSILON)) : _tol(Tol) {}

	// Destructor
	~Structured ();

	// Methods
	size_t Generate (Array<Block*> const & Blocks); ///< Returns the number of elements
	void   WriteVTU (char const * FileName) const;

private:
	// Data
	double         _tol;     ///< Tolerance to remove duplicate nodes
	bool           _is_3d;   ///< Is 3D mesh?
	Array<Vertex*> _verts_d; ///< Vertices (with duplicates)
	Array<Vertex*> _vbry_d;  ///< Vertices on boundary (with duplicates)
	Array<Vertex*> _verts;   ///< Vertices
	Array<Elem*>   _elems;   ///< Elements
	Array<Elem*>   _ebry;    ///< Elements on boundary
	Vector<double> _s;       ///< Current shape (interpolation) values, computed just after _shape(r,s,t)

	// Private methods
	void _shape_2d (double r, double s);
	void _shape_3d (double r, double s, double t);
	void _vtk_con  (Elem const * E, String & Connect) const;

}; // class Structured


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Methods -- Block

inline void Block::Set(Matrix<double> * C, Array<double> * Wx, Array<double> * Wy, Array<double> * Wz)
{
	// Coordinates
	_c     = C;
	_is_3d = (C->Rows()>2 ? true : false);

	// 2D data
	_wx           = Wx;
	_wy           = Wy;
	_n_div_x      = _wx->Size();
	_n_div_y      = _wy->Size();
	_sum_weight_x = 0.0;
	_sum_weight_y = 0.0;
	for (int i=0; i<_n_div_x; ++i) _sum_weight_x += (*_wx)[i];
	for (int i=0; i<_n_div_y; ++i) _sum_weight_y += (*_wy)[i];
	_wx->Push(0.0); // extra value just to allow looping over weights
	_wy->Push(0.0);

	// 3D data
	if (_is_3d)
	{
		if (Wz==NULL) throw new Fatal("Block::Set: For 3D blocks, Wz must be given (not NULL)");
		_wz           = Wz;
		_n_div_z      = _wz->Size();
		_sum_weight_z = 0.0;
		for (int i=0; i<_n_div_z; ++i) _sum_weight_z += (*_wz)[i];
		_wz->Push(0.0);
	}
}

// Destructor -- Structured

inline Structured::~Structured()
{
	for (size_t i=0; i<_verts_d.Size(); ++i) if (_verts_d[i]!=NULL) delete _verts_d[i]; // it is only necessary to delete nodes in _verts_d array
	for (size_t i=0; i<_elems.  Size(); ++i) if (_elems  [i]!=NULL) delete _elems  [i];
}

// Methods -- Structured

inline size_t Structured::Generate(Array<Block*> const & Blocks)
{
	// Check
	if (Blocks.Size()<1) throw new Fatal("Structured::Generate: Number of blocks must be greater than 0 (%d is invalid)",Blocks.Size());

	// Check if the first block is 3D
	_is_3d = Blocks[0]->Is3D();
	_s.Resize((_is_3d ? 20 : 8)); // resize the shape values vector

	// Generate vertices and elements (with duplicates)
	for (size_t b=0; b<Blocks.Size(); ++b)
	{
		// Check if all blocks have the same space dimension
		if (Blocks[b]->Is3D()!=_is_3d) throw new Fatal("Structured::Generate: All blocks must have the same space dimension");
		
		// Generate
		double t = -1.0; // initial Z natural coordinate
		for (int k=0; k<(_is_3d ? Blocks[b]->nDivZ()+1 : 1); ++k)
		{
			double s = -1.0; // initial Y natural coordinate
			for (int j=0; j<Blocks[b]->nDivY()+1; ++j)
			{
				double r = -1.0; // initial X natural coordinate
				for (int i=0; i<Blocks[b]->nDivX()+1; ++i)
				{
					// Compute shape (interpolation) functions
					if (_is_3d) _shape_3d (r,s,t);
					else        _shape_2d (r,s);

					// New vertex
					Vertex * v = new Vertex;
					v->MyID  = _verts_d.Size();     // id
					v->C     = Blocks[b]->C() * _s; // new x-y-z coordinates
					v->Dupl  = false;               // is this a duplicated node?
					if (_is_3d) v->OnBry = (i==0 || j==0 || k==0 || i==Blocks[b]->nDivX() || j==Blocks[b]->nDivY() || k==Blocks[b]->nDivZ());
					else        v->OnBry = (i==0 || j==0 ||         i==Blocks[b]->nDivX() || j==Blocks[b]->nDivY());
					_verts_d.Push(v);
					if (v->OnBry) _vbry_d.Push(v); // on boundary?

					// New element
					if (i!=0 && j!=0 && (_is_3d ? k!=0 : true))
					{
						Elem * e = new Elem;
						e->MyID = _elems.Size(); // id
						if (_is_3d)
						{
							e->V.Resize(8);
							e->V[0] = _verts_d[v->MyID - 1 - (Blocks[b]->nDivX()+1) - (Blocks[b]->nDivX()+1)*(Blocks[b]->nDivY()+1)];
							e->V[1] = _verts_d[v->MyID     - (Blocks[b]->nDivX()+1) - (Blocks[b]->nDivX()+1)*(Blocks[b]->nDivY()+1)];
							e->V[2] = _verts_d[v->MyID                              - (Blocks[b]->nDivX()+1)*(Blocks[b]->nDivY()+1)];
							e->V[3] = _verts_d[v->MyID - 1                          - (Blocks[b]->nDivX()+1)*(Blocks[b]->nDivY()+1)];
							e->V[4] = _verts_d[v->MyID - 1 - (Blocks[b]->nDivX()+1)];
							e->V[5] = _verts_d[v->MyID -     (Blocks[b]->nDivX()+1)];
							e->V[6] = _verts_d[v->MyID];
							e->V[7] = _verts_d[v->MyID - 1];
							e->OnBry = (e->V[0]->OnBry || e->V[1]->OnBry || e->V[2]->OnBry || e->V[3]->OnBry || e->V[4]->OnBry || e->V[5]->OnBry || e->V[6]->OnBry || e->V[7]->OnBry); // if any node is on Bry, yes
							for (size_t m=0; m<8; ++m)
							{
								Share s = {e,m}; // The node shares information will have this element and the local node index
								e->V[m]->Shares.Push(s);
							}
						}
						else
						{
							e->V.Resize(4);
							e->V[0] = _verts_d[v->MyID - 1 - (Blocks[b]->nDivX()+1)];
							e->V[1] = _verts_d[v->MyID     - (Blocks[b]->nDivX()+1)];
							e->V[2] = _verts_d[v->MyID];
							e->V[3] = _verts_d[v->MyID - 1];
							e->OnBry = (e->V[0]->OnBry || e->V[1]->OnBry || e->V[2]->OnBry || e->V[3]->OnBry); // if any node is on Bry, yes
							for (size_t m=0; m<4; ++m)
							{
								Share s = {e,m};
								e->V[m]->Shares.Push(s);
							}
						}
						_elems.Push(e);
						if (e->OnBry) _ebry.Push(e);
					}
					// Next r
					r += (2.0/Blocks[b]->SumWeightX()) * Blocks[b]->Wx(i);
				}
				// Next s
				s += (2.0/Blocks[b]->SumWeightY()) * Blocks[b]->Wy(j);
			}
			// Next t
			t += (_is_3d ? (2.0/Blocks[b]->SumWeightZ()) * Blocks[b]->Wz(k) : 0.0);
		}
	}

	// Remove duplicates
	long ncomp = 0; // number of comparisons
	long ndupl = 0; // number of duplicates
	for (size_t i=0; i<_vbry_d.Size(); ++i)
	{
		for (size_t j=i+1; j<_vbry_d.Size(); ++j)
		{
			// check distance
			double dist = sqrt(          pow(_vbry_d[i]->C(0)-_vbry_d[j]->C(0),2.0)+
										 pow(_vbry_d[i]->C(1)-_vbry_d[j]->C(1),2.0)+
							   (_is_3d ? pow(_vbry_d[i]->C(2)-_vbry_d[j]->C(2),2.0) : 0.0));
			if (dist<_tol)
			{
				/* TODO: this is wrong, since corner nodes can be duplicated and are still on boundary
				// If this node is duplicated, than it is not on-boundary any longer
				_vbry_d[i]->OnBry = false;
				*/
				// Mark duplicated
				_vbry_d[j]->Dupl = true;
				// Chage elements' connectivities
				for (size_t k=0; k<_vbry_d[j]->Shares.Size(); ++k)
				{
					Elem * e = _vbry_d[j]->Shares[k].E;
					int    n = _vbry_d[j]->Shares[k].N;
					e->V[n] = _vbry_d[i];
				}
				ndupl++;
			}
			ncomp++;
		}
	}

	// Set new array with non-duplicated vertices
	size_t k = 0;
	_verts.Resize(_verts_d.Size()-ndupl);
	for (size_t i=0; i<_verts_d.Size(); ++i)
	{
		if (_verts_d[i]->Dupl==false)
		{
			_verts[k]       = _verts_d[i]; // copy pointer
			_verts[k]->MyID = k;           // new ID
			k++;
		}
	}

	std::cout << "number of comparisons = " << ncomp << ", number of duplicates = " << ndupl << std::endl;

	return _elems.Size();
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
	Util::NumStream nsflo = Util::_8s; // number format for floats

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

	/* TODO:
	 * remove this, because, after the removal of duplicates,
	 * the shares information is not valid any longer
	 * I kept this here because it may be useful for debugging purposes
	 * (before the removal of duplicates)*/
	/*
	oss << "        <DataArray type=\"Float32\" Name=\"" << "shares" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << (k==0?"  ":" ") << _verts[i]->Shares.Size();
		k++;
		VTU_NEWLINE (i,k,nn,nfmax,oss);
	}
	oss << "        </DataArray>\n";
	*/

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
		Connect.Printf("%d %d %d %d %d %d %d %d",E->V[1]->MyID,
		                                         E->V[2]->MyID,
		                                         E->V[3]->MyID,
		                                         E->V[0]->MyID,
		                                         E->V[5]->MyID,
		                                         E->V[6]->MyID,
		                                         E->V[7]->MyID,
		                                         E->V[4]->MyID);
	}
	else
	{
		Connect.Printf("%d %d %d %d",E->V[0]->MyID,
		                             E->V[1]->MyID,
		                             E->V[2]->MyID,
		                             E->V[3]->MyID);
	}
}


}; // namespace Mesh


#ifdef USE_BOOST_PYTHON

namespace boopy = boost::python;

class PyMeshBlock
{
public:
	void Set (boopy::list const & C, boopy::list const & Wx, boopy::list const & Wy)
	{
		// Read C
		int nrow = boopy::len(C); if (nrow<1) throw new Fatal("PyMeshBlock: Number of rows of C matrix must be greater than 0 (%d is invalid)",nrow);
		int ncol = boopy::len(C[0]);
		_c.Resize (nrow,ncol);
		for (int i=0; i<nrow; ++i)
		{
			boopy::list row = boopy::extract<boopy::list>(C[i])();
			if (boopy::len(row)!=ncol) throw new Fatal("PyMeshBlock: All rows of C matrix must have the same number of columns (%d is invalid)",boopy::len(row));
			for (int j=0; j<ncol; ++j) _c(i,j) = boopy::extract<double>(row[j])();
		}

		// Read Wx
		int sz_wx = boopy::len(Wx); if (sz_wx<1) throw new Fatal("PyMeshBlock: Number of elements in Wx list must be greater than 0 (%d is invalid)",sz_wx);
		_wx.Resize (sz_wx);
		for (int i=0; i<sz_wx; ++i) _wx[i] = boopy::extract<double>(Wx[i])();

		// Read Wy
		int sz_wy = boopy::len(Wy); if (sz_wy<1) throw new Fatal("PyMeshBlock: Number of elements in Wy list must be greater than 0 (%d is invalid)",sz_wy);
		_wy.Resize (sz_wy);
		for (int i=0; i<sz_wy; ++i) _wy[i] = boopy::extract<double>(Wy[i])();

		// Set _block
		_block.Set (&_c, &_wx, &_wy);
	}
	void Set (boopy::list const & C, boopy::list const & Wx, boopy::list const & Wy, boopy::list const & Wz)
	{
		// Read C
		int nrow = boopy::len(C); if (nrow<1) throw new Fatal("PyMeshBlock: Number of rows of C matrix must be greater than 0 (%d is invalid)",nrow);
		int ncol = boopy::len(C[0]);
		_c.Resize (nrow,ncol);
		for (int i=0; i<nrow; ++i)
		{
			boopy::list row = boopy::extract<boopy::list>(C[i])();
			if (boopy::len(row)!=ncol) throw new Fatal("PyMeshBlock: All rows of C matrix must have the same number of columns (%d is invalid)",boopy::len(row));
			for (int j=0; j<ncol; ++j) _c(i,j) = boopy::extract<double>(row[j])();
		}

		// Read Wx
		int sz_wx = boopy::len(Wx); if (sz_wx<1) throw new Fatal("PyMeshBlock: Number of elements in Wx list must be greater than 0 (%d is invalid)",sz_wx);
		_wx.Resize (sz_wx);
		for (int i=0; i<sz_wx; ++i) _wx[i] = boopy::extract<double>(Wx[i])();

		// Read Wy
		int sz_wy = boopy::len(Wy); if (sz_wy<1) throw new Fatal("PyMeshBlock: Number of elements in Wy list must be greater than 0 (%d is invalid)",sz_wy);
		_wy.Resize (sz_wy);
		for (int i=0; i<sz_wy; ++i) _wy[i] = boopy::extract<double>(Wy[i])();

		// Read Wz
		int sz_wz = boopy::len(Wz); if (sz_wz<1) throw new Fatal("PyMeshBlock: Number of elements in Wz list must be greater than 0 (%d is invalid)",sz_wz);
		_wz.Resize (sz_wz);
		for (int i=0; i<sz_wz; ++i) _wz[i] = boopy::extract<double>(Wz[i])();

		// Set _block
		_block.Set (&_c, &_wx, &_wy, &_wz);
	}
	Mesh::Block * GetBlock () { return &_block; }
private:
	LinAlg::Matrix<double> _c;
	Array<double>          _wx;
	Array<double>          _wy;
	Array<double>          _wz;
	Mesh::Block            _block;
}; // class PyMeshBlock

void (PyMeshBlock::*PMBSet1)(boopy::list const & C, boopy::list const & Wx, boopy::list const & Wy)                         = &PyMeshBlock::Set;
void (PyMeshBlock::*PMBSet2)(boopy::list const & C, boopy::list const & Wx, boopy::list const & Wy, boopy::list const & Wz) = &PyMeshBlock::Set;

class PyMeshStruct
{
public:
	PyMeshStruct ()           : _ms(sqrt(DBL_EPSILON)) {}
	PyMeshStruct (double Tol) : _ms(Tol)               {}
	size_t Generate (boopy::list & ListOfPyMeshBlock)
	{
		int nb = boopy::len(ListOfPyMeshBlock); if (nb<1) throw new Fatal("PyMeshStruct: Number of blocks must be greater than 0 (%d is invalid)",nb);
		Array<Mesh::Block*> blocks;
		blocks.Resize (nb);
		for (int i=0; i<nb; ++i)
		{
			PyMeshBlock * blk = boopy::extract<PyMeshBlock*>(ListOfPyMeshBlock[i])();
			blocks[i] = blk->GetBlock();
		}
		return _ms.Generate (blocks);
	}
	void WriteVTU (boopy::str FileName) { _ms.WriteVTU(boopy::extract<char const *>(FileName)()); }
private:
	Mesh::Structured _ms;
}; // class PyMeshStruct

#endif


#endif // MECHSYS_MESH_STRUCTURED_H
