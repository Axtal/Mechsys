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

#ifndef MECHSYS_FEM_EMBEDDED_H
#define MECHSYS_FEM_EMBEDDED_H

// STL
#include <iostream>
#include <fstream>
#include <cfloat> // for DBL_EPSILON
#include <cstring>
#include <map>

// Boost::Python
#ifdef USE_BOOST_PYTHON
  //#include <boost/python.hpp> // this includes everything
  namespace BPy = boost::python;
#endif // USE_BOOST_PYTHON

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/geometry.h"
#include "util/array.h"
#include "util/numstreams.h"
#include "util/exception.h"
#include "mesh/mesh.h"
#include "mesh/structured.h"

#include "fem/elems/rod3.h"
#include "fem/elems/embspring.h"

using std::cout;
using std::endl;

namespace FEM
{

typedef Array< boost::tuple<double,double,double, char const *,double> > NBrys_T; // Node: x,y,z, key, val
typedef Array< boost::tuple<                 int, char const *,double> > EBrys_T; // Edge:   tag, key, val
typedef Array< boost::tuple<                 int, char const *,double> > FBrys_T; // Face:   tag, key, val
typedef Array< boost::tuple<int, char const*, char const*, char const*,
                            char const*, char const*, bool> > EAtts_T; // Elem: tag, type, model, prms, inis, props, active

inline void GetSegments( double x1, double y1, double z1,   ///< In: Initial point coordinates
                         double x2, double y2, double z2,   ///< In: Final   point coordinates
                         int nNodes,                        ///< In: Number of nodes for the new rod elements
                         FEM::Geom           * G,           ///< Out: The FE geometry
                         Array<int> & aElems,               ///< Out: Array of tresspassed elements
                         Array<Array<double> > & aSegments) ///< Out: Array of arrays containig nodal coordinates of the new nodes
{
	// Segments disposition Array<Array<double> >
	// { {s1x1 s1y1 s1z1 s1x2 s1y2 s1z2 } {s2x1 s2y1 s2z1 s2x2 s2y2 s2z2 } ... }

	// Determination of embedded portions
	//int    n_elems = DA.Elements.Size();
	int    n_elems = G->NElems();
	double tiny    = 1E-4;

	double length      = sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));
	double tinylength  = 0.01*length;
	double l           = (x2-x1)/length;
	double m           = (y2-y1)/length;
	double n           = (z2-z1)/length;
	Vector<double>       Direction(3);
	Direction          = l, m, n;

	int    init_elem   = 0;
	int    prev_elem   = 0;
	int    actual_elem = 0;
	int    next_elem   = 0;
	int    final_elem  = 0;
	//find the initial element
	for (int j=0; j<n_elems; j++)
		//if (DA.Elements[j]->IsInside(x1+tinylength*l, y1+tinylength*m, z1+tinylength*n))
		if (G->Ele(j)->IsInside(x1+tinylength*l, y1+tinylength*m, z1+tinylength*n)) { init_elem = j; break; }

	//find the final element
	for (int j=0; j<n_elems; j++)
		if (G->Ele(j)->IsInside(x2-tinylength*l, y2-tinylength*m, z2-tinylength*n)) { final_elem = j; break; }

	double xp = x1;  // previous x coordinate
	double yp = y1;  // previous y coordinate
	double zp = z1;  // previous z coordinate
	double x  = xp;  // current  x coordinate
	double y  = yp;  // current  y coordinate
	double z  = zp;  // current  z coordinate
	actual_elem  = prev_elem = init_elem;
	next_elem    = 0;
	double step  = length;
	bool initial = true;  // flag for the initial element
	bool   final = false; // flag for the final element

	if (init_elem == final_elem) { x=x2; y=y2; z=z2; final=true; }

	//Splitting bar:
	do
	{
		if (!final)
		{
			step *= 0.501;
			x    += step*l;
			y    += step*m;
			z    += step*n;
			//if (DA.Elements[prev_elem]->IsInside(x, y, z))
			if (G->Ele(prev_elem)->IsInside(x, y, z))
				actual_elem = prev_elem;
			else
				for (int j=0; j<n_elems; j++)
					if (G->Ele(j)->IsInside(x, y, z)) { actual_elem = j; break; }
		}
		if (prev_elem == actual_elem)
		{
			double bf;  // boundary function value
			if (!final)
			{
				double r, s, t;
				//DA.Elements[next_elem]->InverseMap(x,y,z,r,s,t);
				G->Ele(next_elem)->InverseMap(x,y,z,r,s,t);
				//bf = fabs(DA.Elements[next_elem]->BoundDistance(r,s,t));
				bf = fabs(G->Ele(next_elem)->BoundDistance(r,s,t));
			}
			if ( final || bf <= tiny )  // intersection is reached
			{
				Array<double> segment;

				// Initialize connectivities of element
				int n_div = nNodes - 1;
				for (int j=0; j<nNodes; ++j)
				{
					double _x = xp + j*(x-xp)/n_div;
					double _y = yp + j*(y-yp)/n_div;
					double _z = zp + j*(z-zp)/n_div;
					segment.Push(_x);
					segment.Push(_y);
					segment.Push(_z);
				}
				aElems   .Push(actual_elem);
				aSegments.Push(segment);

				if (final) break; //if this is the final element
				else{
					xp=x; yp=y; zp=z;
					if (next_elem==final_elem)
					{
						x=x2; y=y2; z=z2;
						actual_elem = prev_elem = next_elem;
						final = true;
					}
					else
					{
						prev_elem = next_elem;
						step = sqrt(pow(x-x2,2)+pow(y-y2,2)+pow(z-z2,2));
					}
				}
				initial = false;
			}
			else{ step=fabs(step); } // Advance: when intersection is forward
		}
		else{ step=-fabs(step); next_elem = actual_elem; } // Retrocession: when intersection is backward
	}while (true);
}

void AddReinf(double x1, double y1, double z1, // In:  Coordinates of the initial point
              double x2, double y2, double z2, // In:  Coordinates of the final point
              char const *          Prms,      // In:  Properties and parameters [E Ar At ks]
              bool                  IsActive,  // In:  Define if the entire reinforcement is active or not
              int                   Tag,       // In:  Tag for all generated elements
              FEM::Geom  *          G)         // Out: The FE geometry
{
	// E=2.0E8 Ar=0.02 ks=2.0E6
	LineParser    lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	double E = -1.0, Ar = -1.0, At = -1, ks = -1.0;
	double c = 0.0,  phi = 0.0;
	for (size_t i=0; i<names.Size(); ++i)
	{
		     if (names[i]=="E" )  E   = values[i];
		else if (names[i]=="Ar")  Ar  = values[i];
		else if (names[i]=="At")  At  = values[i];
		else if (names[i]=="ks")  ks  = values[i];
		else if (names[i]=="c")   c   = values[i];
		else if (names[i]=="phi") phi = values[i];
		else throw new Fatal("AddReinforcement: Parameter name (%s) is invalid",names[i].CStr());
	}

	if (At<0) At=Ar;

	// Check
	if (E  <=0.0) throw new Fatal("AddReinf: Young modulus (E) must be provided (and positive). E==%f is invalid",E);
	if (Ar <=0.0) throw new Fatal("AddReinf: Steel Section Area (Ar) must be provided (and positive). Ar==%f is invalid",Ar);
	if (At <=0.0) throw new Fatal("AddReinf: Total Section Area (At) must be provided (and positive). At==%f is invalid",At);
	if (ks <=0.0) throw new Fatal("AddReinf: Spring Stiffness (ks) must be provided (and positive). ks=%f is invalid",ks);

	// Determine the internal segments
	Array<Array<double> > a_segments;  // { {E1x1 E1y1 E1z1 E1x2 E1y2 E1z2 } {E2x1 E2y1 E2z1 E2x2 E2y2 E2z2 } ... }
	Array<int>            a_elems;     // { E1, E2, E3, ... }
	int n_seg_nodes = 3;               // Second order
	GetSegments(x1, y1, z1, x2, y2, z2, n_seg_nodes, G, a_elems, a_segments);

	// Direction
	double length    = sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));
	double l         = (x2-x1)/length;
	double m         = (y2-y1)/length;
	double n         = (z2-z1)/length;
	Vector<double>     Direction(3); Direction = l, m, n;

	// Assign weigths for interface distribution
	Vector<double> weights;
	assert(n_seg_nodes<=4);
	     if (n_seg_nodes==2) { weights.Resize(2); weights = 0.5, 0.5; }
	else if (n_seg_nodes==3) { weights.Resize(3); weights = 0.25, 0.5, 0.25; } //weights = 0.16666666, 0.67777777, 0.16666666;
	else if (n_seg_nodes==4) { weights.Resize(4); weights = 1./6., 1./3., 1./3., 1./6.; }

	// Loop along the segments
	Node * last_end_node = 0; // Pointer to the end node of the last segment
	for (size_t i_seg=0; i_seg<a_segments.Size(); ++i_seg)
	{
		// Flags
		bool initial_segment = (i_seg>0)?false:true;
		bool final_segment   = (i_seg == a_segments.Size()-1)?true:false;
		Array<double> & c_segm = a_segments[i_seg];

		// segment extreme nodes coordinates
		double x1 = c_segm[0];
		double y1 = c_segm[1];
		double z1 = c_segm[2];
		double x2 = c_segm[n_seg_nodes*3-3];
		double y2 = c_segm[n_seg_nodes*3-2];
		double z2 = c_segm[n_seg_nodes*3-1];

		double L  = sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2)); // lenght of the segment
		double LA = 2.0*L*sqrt(3.14159*At);                       // total interface area (lateral area)

		//Array for the new nodes
		Array<Node*> nodes;
		nodes.Resize(n_seg_nodes);

		// Adding the initial and final nodes
		Node * begin_node = new Node;
		Node *   end_node = new Node;
		G->PushNode(begin_node); begin_node->Initialize(G->NNodes()-1, x1, y1, z1);
		G->PushNode(  end_node);   end_node->Initialize(G->NNodes()-1, x2, y2, z2);

		// Allocate bar reinforcement
		Element * bar_elem = new Rod3; // Three nodes Rod
		G->PushElem(bar_elem);
		bar_elem->Initialize(G->NElems()-1, IsActive, G->NDim(), Tag);

		// Loop along the nodes from the segment
		for (int j_node=0; j_node<n_seg_nodes; j_node++)
		{
			// Generating the bar element:

			// Flags
			bool initial_node = (j_node>0)?false:true;
			bool final_node   = (j_node==n_seg_nodes-1)?true:false;

			// Coordinates of the current node
			double x = c_segm[j_node*3+0];
			double y = c_segm[j_node*3+1];
			double z = c_segm[j_node*3+2];

			if (initial_node) // first node of the current segment
			{
				if (initial_segment) nodes[j_node] = begin_node;
				else                 nodes[j_node] = last_end_node;
			}
			else if (final_node && final_segment)
				nodes[j_node] = end_node; // final node of the final segment
			else
			{
				nodes[j_node] = new Node;
				G->PushNode(nodes[j_node]); nodes[j_node]->Initialize(G->NNodes()-1, x, y, z);
				last_end_node = nodes[j_node];
			}
			bar_elem->Connect(j_node, nodes[j_node]);

			// Generating the interface element:
			Element * solid_elem = G->Ele(a_elems[i_seg]);
			Element * conn_elem = new EmbSpring(solid_elem, nodes[j_node], Direction);
			G->PushElem(conn_elem);
			conn_elem->Initialize(G->NElems()-1, IsActive, G->NDim(), Tag);

			// Connector (EmbSpring) connectivities
			for (size_t i=0; i<solid_elem->NNodes(); i++)
				conn_elem->Connect(i, solid_elem->Nod(i));           // Add the nodes from the solid element
			conn_elem->Connect(solid_elem->NNodes(), nodes[j_node]); // Add the node attached to the bar

			// Set parameters for the connector spring
			std::ostringstream os;
			os << "ks=" << ks << " Al=" << LA*weights(j_node);
			conn_elem->SetModel  ("", os.str().c_str(), "");
		}

        // Set parameters for the bar element
		std::ostringstream os;
		os << "E=" << E << " A=" << Ar;
		bar_elem->SetModel  ("", os.str().c_str(), "");
	}
}

} // Namespace FEM

#ifdef USE_BOOST_PYTHON

void PyAddReinf(double x1, double y1, double z1, // In:  Coordinates of the initial point
                double x2, double y2, double z2, // In:  Coordinates of the final point
                BPy::str const      & Prms,      // In:  Properties and parameters [E Ar At ks]
                bool                  IsActive,  // In:  Define if the entire reinforcement is active or not
                int                   Tag,       // In:  Tag for all generated elements
                FEM::Geom           & G)         // Out: The FE geometry
{
	AddReinf (x1,y1,z1, x2,y2,z2, BPy::extract<char const *>(Prms)(), IsActive, Tag, &G);
}

void PyAddReinfs (bool              Is3D,  ///< 3D ?
                  BPy::dict const & Verts, ///< {IDv1:(x1,y1,z1), IDv2:(x2,y2,z2), ... num vertices}
                  BPy::dict const & Edges, ///< {(v1,v2):tag1, (v3,v4):tag2, ... num edges} v# => vertex number
                  BPy::list const & EAtts, ///< Elements attributes
                  FEM::Geom       &     G) ///< Out: The FE geometry
{
	/* Example:
	 *           
	 *           # Elements attributes
	 *           eatts = [[-1, 'Reinforcement', '', 'E=%g Ar=%g At=%g ks=%g c=%g phi=%g', 'ZERO', 'gam=20', True]] # tag, type, model, prms, inis, props, active?
	 */

	// Map element tag to index in EAtts list
	int neatts = BPy::len(EAtts);
	if (neatts<1) throw new Fatal("embedded.h::PyAddReinfs: EAtts, Element(reinforcement) Attributes, must contain at least one element");
	std::map<int,int> tag2idx; 
	for (int i=0; i<neatts; ++i)
	{
		BPy::list const & lst = BPy::extract<BPy::list>(EAtts[i])();
		tag2idx[BPy::extract<int>(lst[0])()] = i;
		if (BPy::len(EAtts[i])!=7) throw new Fatal("embedded.h::PyAddReinfs: Each sublist in EAtts must have 7 items: tag, type, model, prms, inis, props, active?\n\tExample: eatts = [[-1, 'Reinforcement', '', 'E=200 Ar=1 At=1 ks=1e+12 c=0 phi=20', 'ZERO', 'gam=20', True]]\n\tlen(EAtts[i])==%d is invalid.",BPy::len(EAtts[i]));
	}

	// Read edges
	double x0=0.0; double y0=0.0; double z0=0.0;
	double x1=0.0; double y1=0.0; double z1=0.0;
	BPy::object const & e_keys = BPy::extract<BPy::dict>(Edges)().iterkeys();
	BPy::object const & e_vals = BPy::extract<BPy::dict>(Edges)().itervalues();
	for (int i=0; i<BPy::len(Edges); ++i)
	{
		// Extract reinforcement data
		BPy::tuple const & edge = BPy::extract<BPy::tuple> (e_keys.attr("next")())();
		int                tag  = BPy::extract<int>        (e_vals.attr("next")())();
		long               v0   = BPy::extract<long>       (edge[0])();
		long               v1   = BPy::extract<long>       (edge[1])();
		BPy::tuple const & xyz0 = BPy::extract<BPy::tuple> (Verts[v0])();
		BPy::tuple const & xyz1 = BPy::extract<BPy::tuple> (Verts[v1])();
		x0 = BPy::extract<double>(xyz0[0])();   y0 = BPy::extract<double>(xyz0[1])();   if (Is3D) z0 = BPy::extract<double>(xyz0[2])();
		x1 = BPy::extract<double>(xyz1[0])();   y1 = BPy::extract<double>(xyz1[1])();   if (Is3D) z1 = BPy::extract<double>(xyz1[2])();
		//std::cout << "Edge: tag="<<tag << ", v0="<<v0 << ", v1="<<v1 << ", x0="<<x0 << ", y0="<<y0 << ", z0="<<z0 << ", x1="<<x1 << ", y1="<<y1 << ", z1="<<z1 << std::endl;
		
		// Find element attributes
		std::map<int,int>::const_iterator iter = tag2idx.find(tag);
		if (iter==tag2idx.end()) throw new Fatal("embedded.h::PyAddReinfs: Could not find tag < %s > in the list of Element Attributes", tag);
		int idx_eatt = iter->second;

		// Add reinforcement to FE geometry
		BPy::str const & prms = BPy::extract<BPy::str>(EAtts[idx_eatt][3])();
		bool             act  = BPy::extract<bool>    (EAtts[idx_eatt][6])();
		AddReinf (x0,y0,z0, x1,y1,z1, BPy::extract<char const *>(prms)(), act, tag, &G);
	}
}

#endif // USE_BOOST_PYTHON

#endif // MECHSYS_FEM_EMBEDDED_H
