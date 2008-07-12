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

#ifndef MECHSYS_FEM_FUNCTIONS_H
#define MECHSYS_FEM_FUNCTIONS_H

// STL
#include <iostream>
#include <fstream>
#include <cfloat>  // for DBL_EPSILON

// Boost
#ifdef USE_BOOST
  #include <boost/tuple/tuple_io.hpp>
#endif

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/geometry.h"
#include "util/array.h"
#include "util/numstreams.h"
#include "util/exception.h"
#include "mesh/structured.h"

#ifndef VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

namespace FEM
{


#if defined(USE_BOOST) || defined(USE_BOOST_PYTHON)

typedef Array< boost::tuple<double,double,double, char const *,double> >               NBrys_T; // x,y,z, key, val
typedef Array< boost::tuple<                 int, char const *,double> >               FBrys_T; //   tag, key, val
typedef Array< boost::tuple<int, char const*, char const*, char const*, char const*> > EAtts_T; // tag, type, model, prms, inis

inline void SetGeom (Mesh::Structured const * M, NBrys_T const & NodesBrys, FBrys_T const & FacesBrys, EAtts_T const & ElemsAtts, FEM::Geom * G)
{
	/* Example:
	
		// 1) Nodes brys
		FEM::NBrys_T nbrys;
		nbrys.Push (make_tuple(L/2., 0.0, 0.0, "ux", 0.0)); // x,y,z, key, val

		// 2) Faces brys
		FEM::FBrys_T fbrys;
		fbrys.Push (make_tuple(-10, "uy", 0.0)); // tag, key, val
		fbrys.Push (make_tuple(-20, "fy",  -q)); // tag, key, val

		// 3) Elements attributes
		FEM::EAtts_T eatts;
		eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", "E=207 nu=0.3", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0")); // tag, type, model, prms, inis

	*/

	// 3D mesh?
	bool is3d = M->Is3D();

	// Set nodes
	size_t nn = M->Verts().Size();
	G->SetNNodes (nn);
	for (size_t i=0; i<nn; ++i) // loop over all vertices
	{
		// Current vertex
		Mesh::Vertex const * v = M->Verts()[i];

		// New node
		G->SetNode (v->MyID, v->C(0), v->C(1), (is3d ? v->C(2) : 0.0));
	}

	// Set elements
	size_t ne = M->Elems().Size();
	G->SetNElems (ne);
	for (size_t i=0; i<ne; ++i)
	{
		// Current mesh element
		Mesh::Elem const * me = M->Elems()[i];

		// Set element
		bool found = false;
		for (size_t j=0; j<ElemsAtts.Size(); ++j)
		{
			if (me->Tag==ElemsAtts[j].get<0>())
			{
				// New finite element
				found = true;
				FEM::Element * fe = G->SetElem (me->MyID, ElemsAtts[j].get<1>());

				// Set connectivity
				for (size_t k=0; k<me->V.Size(); ++k)
					fe->Connect (k, G->Nod(me->V[k]->MyID));

				// Set parameters and initial values
				fe->SetModel (ElemsAtts[j].get<2>(), ElemsAtts[j].get<3>(), ElemsAtts[j].get<4>());
				break;
			}
		}
		if (found==false) throw new Fatal("SetGeom: Could not find Tag==%d for Element %d in the ElemsAtts list",me->Tag,i);
	}

	// Set faces boundaries
	if (FacesBrys.Size()>0)
	{
		for (size_t i=0; i<M->ElemsBry().Size(); ++i) // loop over all elements on boundary
		{
			// Current mesh element on boundary
			Mesh::Elem const * me = M->ElemsBry()[i];

			for (int j=0; j<me->ETags.Size(); ++j) // j is the local face id
			{
				int tag = me->ETags(j);
				if (tag<0) // this element has a face tag
				{
					bool found = false;
					for (size_t k=0; k<FacesBrys.Size(); ++k)
					{
						if (tag==FacesBrys[k].get<0>())
						{
							found = true;
							G->Ele(me->MyID)->Bry (FacesBrys[k].get<1>(), FacesBrys[k].get<2>(), j);
							break;
						}
					}
					if (found==false) throw new Fatal("SetGeom: Could not find Tag==%d for Face %d of Element %d in the FacesBrys list",tag,j,me->MyID);
				}
			}
		}
	}

	// Set nodes boundaries
	if (NodesBrys.Size()>0)
	{
		for (size_t i=0; i<M->VertsBry().Size(); ++i) // loop over all vertices on boundary
		{
			// Current vertex
			Mesh::Vertex const * v = M->VertsBry()[i];

			for (size_t j=0; j<NodesBrys.Size(); ++j)
			{
				double x =         NodesBrys[j].get<0>();
				double y =         NodesBrys[j].get<1>();
				double z = (is3d ? NodesBrys[j].get<2>() : 0.0);
				double d = sqrt(pow(x - v->C(0),2.0) + pow(y - v->C(1),2.0) + (is3d ? pow(z - v->C(2),2.0) : 0.0));
				if (d<sqrt(DBL_EPSILON))
					G->Nod(v->MyID)->Bry (NodesBrys[j].get<3>(), NodesBrys[j].get<4>());
			}
		}
	}
}

#endif // USE_BOOST || USE_BOOST_PYTHON


inline void _write_elem_val (FEM::Geom const & G, size_t ne, size_t nfmax, char const * Key, std::ostringstream & oss)
{
	oss << "        <DataArray type=\"Float32\" Name=\""<< Key <<"\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		double val = 0.0;
		try { val = G.Ele(i)->Val(Key); } catch (Exception * e) { delete e; }
		oss << (k==0?"  ":" ") << val;
		k++;
		VTU_NEWLINE (i,k,ne,nfmax,oss);
	}
	oss << "        </DataArray>\n";
}

inline void WriteVTUEquilib (FEM::Geom const & G, char const * FileName)
{
	// Open File
	std::ofstream      of(FileName, std::ios::out);
	std::ostringstream oss;

	// Data
	size_t nn = G.nNodes(); // Number of Nodes
	size_t ne = G.nElems(); // Number of Elements

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
		oss << "  " << nsflo << G.Nod(i)->X() << " ";
		oss <<         nsflo << G.Nod(i)->Y() << " ";
		oss <<         nsflo << G.Nod(i)->Z();
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
		String con;  G.Ele(i)->VTKConnect(con);
		oss << "  " << con;
		k++;
		VTU_NEWLINE (i,k,ne,nimax/G.Ele(i)->nNodes(),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t ossfset = 0;
	for (size_t i=0; i<ne; ++i)
	{
		ossfset += G.Ele(i)->nNodes();
		oss << (k==0?"  ":" ") << ossfset;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << G.Ele(i)->VTKCellType();
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Cells>\n";

	// Data -- nodes
	oss << "      <PointData Vectors=\"TheVectors\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "displ" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << " " << nsflo << (G.Nod(i)->HasVar("ux") ? G.Nod(i)->Val("ux"): 0.0) << " ";
		oss <<        nsflo << (G.Nod(i)->HasVar("uy") ? G.Nod(i)->Val("uy"): 0.0) << " ";
		oss <<        nsflo << (G.Nod(i)->HasVar("uz") ? G.Nod(i)->Val("uz"): 0.0) << " ";
		k++;
		VTU_NEWLINE (i,k,nn,nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "force" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << " " << nsflo << (G.Nod(i)->HasVar("fx") ? G.Nod(i)->Val("fx"): 0.0) << " ";
		oss <<        nsflo << (G.Nod(i)->HasVar("fy") ? G.Nod(i)->Val("fy"): 0.0) << " ";
		oss <<        nsflo << (G.Nod(i)->HasVar("fz") ? G.Nod(i)->Val("fz"): 0.0) << " ";
		k++;
		VTU_NEWLINE (i,k,nn,nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </PointData>\n";

	// Data -- elements
	oss << "      <CellData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"active\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << G.Ele(i)->IsActive();
		k++;
		VTU_NEWLINE (i,k,ne,nfmax,oss);
	}
	oss << "        </DataArray>\n";
	_write_elem_val(G, ne, nfmax, "Sx", oss);
	_write_elem_val(G, ne, nfmax, "Sy", oss);
	_write_elem_val(G, ne, nfmax, "Sz", oss);
	_write_elem_val(G, ne, nfmax, "Sxy",oss);
	_write_elem_val(G, ne, nfmax, "Syz",oss);
	_write_elem_val(G, ne, nfmax, "Szx",oss);
	_write_elem_val(G, ne, nfmax, "p",  oss);
	_write_elem_val(G, ne, nfmax, "q",  oss);
	_write_elem_val(G, ne, nfmax, "Ex", oss);
	_write_elem_val(G, ne, nfmax, "Ey", oss);
	_write_elem_val(G, ne, nfmax, "Ez", oss);
	_write_elem_val(G, ne, nfmax, "Exy",oss);
	_write_elem_val(G, ne, nfmax, "Eyz",oss);
	_write_elem_val(G, ne, nfmax, "Ezx",oss);
	_write_elem_val(G, ne, nfmax, "Ev", oss);
	_write_elem_val(G, ne, nfmax, "Ed", oss);
	oss << "      </CellData>\n";

	// Bottom
	oss << "    </Piece>\n";
	oss << "  </UnstructuredGrid>\n";
	oss << "</VTKFile>" << std::endl;

	// Write to file
	of << oss.str();
	of.close();
}

inline void WriteVTK (FEM::Geom const & G, char const * FileName)
{
	// Filter elements
	Array<Element const*> act_elems; // Array for active elements
	for (size_t i=0; i<G.nElems(); ++i)
	{
		if (G.Ele(i)->IsActive()) // Only active elements are considered
			act_elems.Push(G.Ele(i));
	}

	// Data
	size_t n_nodes = G.nNodes();       // Number of Nodes
	size_t n_elems = act_elems.Size(); // Number of Elements
	std::map<String, int>  index_map;  // Map to associate labels with indexes

	// Get all possible labels from elements
	for (size_t i_elem=0; i_elem<n_elems; ++i_elem)
	{
		Array<String>  elem_labels;
		act_elems[i_elem]->GetLabels(elem_labels);
		int n_labels = elem_labels.Size();
		for (int j_label=0; j_label<n_labels; ++j_label)
		{
			String & current_label = elem_labels[j_label];
			if (index_map.find(current_label)==index_map.end())
				index_map[current_label] = index_map.size()-1; // add a new entry
		}
	}

	// Collect nodal values
	size_t n_comps = index_map.size();
	LinAlg::Matrix<double> values(n_nodes, n_comps); values.SetValues(0.0);  // Matrix for nodal values collected from elements
	LinAlg::Matrix<size_t> refs  (n_nodes, n_comps); refs  .SetValues(0);    // Matrix for nodal references of variablea
	for (size_t i_elem=0; i_elem<n_elems; ++i_elem)
	{
		LinAlg::Matrix<double> elem_values;
		Array<String>          elem_labels;
		act_elems[i_elem]->OutNodes(elem_values, elem_labels);
		int                    n_labels     = elem_labels.Size();
		int                    n_elem_nodes = act_elems[i_elem]->nNodes();
		for (int j_label=0; j_label<n_labels; ++j_label)
		{
			String & current_label = elem_labels[j_label];
			int index = index_map[current_label];
			for (int j_node=0; j_node<n_elem_nodes; ++j_node)
			{
				int node_number            = act_elems[i_elem]->Nod(j_node)->GetID();
				values(node_number,index) += elem_values(j_node, j_label);  // accumulate values
				refs  (node_number,index) ++;                               // update references number
			}
		}
	}
	
	// Compute average values
	for (size_t i=0; i<n_nodes; i++)
	for (size_t j=0; j<n_comps; j++)
	{
		if   (refs(i,j)!=0) values(i,j) /= refs(i,j);
		else                values(i,j)  = 0.0;
	}
	
	// Total number of CELLS data = Sum (1 + nNodes); 1 for the numPts label
	int n_data = 0;
	for (size_t i=0; i<n_elems; ++i) n_data += 1 + act_elems[i]->nNodes();

	// Define variables for displacements
	const String UX = "ux";
	const String UY = "uy";
	const String UZ = "uz";
	
	// Define variables for velocity
	const String VX = "vx";
	const String VY = "vy";
	const String VZ = "vz";

	// Check if exists data about displacements
	bool has_disp = false;
	if (index_map.find(UX)!=index_map.end()) has_disp = true;
	if (index_map.find(UY)!=index_map.end()) has_disp = true;
	if (index_map.find(UZ)!=index_map.end()) has_disp = true;

	// Check if exists data about velocities
	bool has_vel = false;
	if (index_map.find(VX)!=index_map.end()) has_vel  = true;
	if (index_map.find(VY)!=index_map.end()) has_vel  = true;
	if (index_map.find(VZ)!=index_map.end()) has_vel  = true;

	// Structure for output and float number output format
	std::ostringstream oss;
	Util::NumStream nsflo = Util::_8s; // number format for floats

	// Write Legacy VTK file header
	oss << "# vtk DataFile Version 3.0"     << std::endl;
	oss << "MechSys/FEM - "                 << std::endl;
	oss << "ASCII"                          << std::endl;
	oss << "DATASET UNSTRUCTURED_GRID"      << std::endl;
	oss <<                                     std::endl;

	// Node coordinates
	oss << "POINTS " << n_nodes << " float" << std::endl;
	for (size_t i=0; i<n_nodes; ++i)
		oss << nsflo << G.Nod(i)->X() << nsflo << G.Nod(i)->Y() << nsflo << G.Nod(i)->Z() << std::endl;
	oss << std::endl;

	// Elements connectivities
	oss << "CELLS "<< n_elems << " " << n_data << std::endl;
	for (size_t i=0; i<n_elems; ++i)
	{
		String connect; act_elems[i]->VTKConnect(connect);
		oss << act_elems[i]->nNodes() << " " << connect << std::endl;
	}
	oss << std::endl;

	// Cell types
	oss << "CELL_TYPES " << n_elems << std::endl;
	for (size_t i=0; i<n_elems; ++i)
		oss << act_elems[i]->VTKCellType() << std::endl;
	oss << std::endl;

	oss << "POINT_DATA " << n_nodes << std::endl;

	// Vectors
	if (has_disp)
	{
		oss << "VECTORS " << "Disp float" << std::endl;
		for (size_t j=0; j<n_nodes; ++j)
		{
			oss << nsflo << values(j, index_map[UX]) << nsflo << values(j, index_map[UY])
			    << nsflo << ((index_map.find(UZ)==index_map.end())?0.0:values(j, index_map[UZ])) << std::endl;
		}
		oss << std::endl;
	}
	if (has_vel)
	{
		oss << "VECTORS " << "Vel float" << std::endl;
		for (size_t j=0; j<n_nodes; ++j)
		{
			oss << nsflo << values(j, index_map[VX]) << nsflo << values(j, index_map[VY])
			    << nsflo << ((index_map.find(VZ)==index_map.end())?0.0:values(j, index_map[VZ])) << std::endl; 
		}
		oss << std::endl;
	}

	// Scalars
	std::map<String,int>::iterator iter;
	for (iter=index_map.begin(); iter!=index_map.end(); iter++)
	{
		oss << "SCALARS " << iter->first << " float 1" << std::endl;
		oss << "LOOKUP_TABLE default" << std::endl;
		for (size_t j=0; j<n_nodes; ++j)
			oss << nsflo << values(j,iter->second) << std::endl;
		oss << std::endl;
	}

	// Create file and copy contents of 'oss' into it
	std::ofstream ofile;
	ofile.open (FileName, std::ios::out);
	ofile << oss.str();
	ofile.close();
}

}; // namespace FEM


#ifdef USE_BOOST_PYTHON
// {

namespace boopy = boost::python;

void PyWriteVTUEquilib (PyGeom const & G, boopy::str const & FileName) { FEM::WriteVTUEquilib((*G.GetGeom()), boopy::extract<char const *>(FileName)()); }
void PyWriteVTK        (PyGeom const & G, boopy::str const & FileName) { FEM::WriteVTK       ((*G.GetGeom()), boopy::extract<char const *>(FileName)()); }

void PySetGeom (PyMeshStruct const & M, boopy::list const & NodesBrys, boopy::list const & FacesBrys, boopy::list const & ElemsAtts, PyGeom & G)
{
	/* Example:
	 *           # 1) Nodes brys
	 *           nbrys = [[L/2., 0.0, 0.0, 'ux', 0.0]] # x,y,z, key, val
	 *           
	 *           # 2) Faces brys
	 *           fbrys = [[-10, 'uy', 0.0], # [tag], [key], [val]
	 *                    [-20, 'fy',  -q]] # [tag], [key], [val]
	 *           
	 *           # 3) Elements attributes
	 *           eatts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=%f nu=%f'%(E,nu), 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0']] # [tag], [type], [model], [prms], [inis]
	 */

	// Extract list with nodes boundaries
	FEM::NBrys_T nbrys;
	int nbrys_size = len(NodesBrys);
	if (nbrys_size>0) nbrys.Resize(nbrys_size);
	for (int i=0; i<nbrys_size; ++i)
	{
		if (len(NodesBrys[i])==5)
		{
			boopy::list lst = boopy::extract<boopy::list>(NodesBrys[i])();
			nbrys[i] = boost::make_tuple(boopy::extract<double     >(lst[0])(),
			                             boopy::extract<double     >(lst[1])(),
			                             boopy::extract<double     >(lst[2])(),
			                             boopy::extract<char const*>(lst[3])(),
			                             boopy::extract<double     >(lst[4])());
		}
		else throw new Fatal("PySetGeom: Each sublist in NodesBrys must have 5 items: x,y,z, key, val\n\tExample: NodesBrys = [[1.0, 0.0, 0.0, 'ux', 0.0]]");
	}

	// Extract list with faces boundaries
	FEM::FBrys_T fbrys;
	int fbrys_size = len(FacesBrys);
	if (fbrys_size>0) fbrys.Resize(fbrys_size);
	for (int i=0; i<fbrys_size; ++i)
	{
		if (len(FacesBrys[i])==3)
		{
			boopy::list lst = boopy::extract<boopy::list>(FacesBrys[i])();
			fbrys[i] = boost::make_tuple(boopy::extract<int        >(lst[0])(),
			                             boopy::extract<char const*>(lst[1])(),
			                             boopy::extract<double     >(lst[2])());
		}
		else throw new Fatal("PySetGeom: Each sublist in FacesBrys must have 3 items: tag, key, val\n\tExample: FacesBrys = [[-10, 'uy', 0.0], [-20, 'fy', -1]]");
	}

	// Extract list with elements attributes
	FEM::EAtts_T eatts;
	int eatts_size = len(ElemsAtts);
	if (eatts_size>0) eatts.Resize(eatts_size);
	for (int i=0; i<eatts_size; ++i)
	{
		if (len(ElemsAtts[i])==5)
		{
			boopy::list lst = boopy::extract<boopy::list>(ElemsAtts[i])();
			eatts[i] = boost::make_tuple(boopy::extract<int        >(lst[0])(),
			                             boopy::extract<char const*>(lst[1])(),
			                             boopy::extract<char const*>(lst[2])(),
			                             boopy::extract<char const*>(lst[3])(),
			                             boopy::extract<char const*>(lst[4])());
		}
		else throw new Fatal("PySetGeom: Each sublist in ElemsAtts must have 5 items: tag, type, model, prms, inis\n\tExample: ElemsAtts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=207.0 nu=0.3', 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0']]");
	}

	// Set geometry
	FEM::SetGeom (M.GetMesh(), nbrys, fbrys, eatts, G.GetGeom());
}

// }
#endif // USE_BOOST_PYTHON


#endif // MECHSYS_FEM_FUNCTIONS_H
