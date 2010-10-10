/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_HYDROMECH_H
#define MECHSYS_HYDROMECH_H

// MechSys
#include <mechsys/fem/equilibelem.h>

namespace FEM
{

class HydroMechElem : public EquilibElem
{
public:
    // Static
    static size_t NDp; ///< Number of DOFs of pressure = GEp->NN

    // Constructor
    HydroMechElem (int                  NDim,   ///< Space dimension
                   Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                   Model        const * Mdl,    ///< Model
                   SDPair       const & Prp,    ///< Properties
                   SDPair       const & Ini,    ///< Initial values
                   Array<Node*> const & Nodes); ///< Connectivity

    // Destructor
    ~HydroMechElem () { if (GEp!=NULL) delete GEp; }

    // Data
    GeomElem * GEp; ///< Element for the pressure DOFs
};

size_t HydroMechElem::NDp = 0;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline HydroMechElem::HydroMechElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : EquilibElem(NDim,Cell,Mdl,Prp,Ini,Nodes)
{
    // allocate GEp
    if (GE->Name=="Hex20") GEp = AllocGeomElem ("Hex8", NDim);
    else throw new Fatal("HydroMechElem::HydroMechElem: Geom element==%s is not implemented yet",GE->Name.CStr());

    // set constants of this class (just once)
    if (NDp==0)
    {
        NDp = GEp->NN;
    }

    // set initial values at Nodes
    if (Ini.HasKey("geostatic"))
    {
        if (!Ini.HasKey("gamW"))               throw new Fatal("HydroMechElem::HydroMechElem: For geostatic stresses, 'gamW' must be provided in 'Ini' dictionary");
        if (NDim==2 && !Ini.HasKey("y_water")) throw new Fatal("HydroMechElem::HydroMechElem: For geostatic stresses in 2D, 'y_water' must be provided in 'Ini' dictionary");
        if (NDim==3 && !Ini.HasKey("z_water")) throw new Fatal("HydroMechElem::HydroMechElem: For geostatic stresses in 3D, 'z_water' must be provided in 'Ini' dictionary");
        bool   pos_pw  = (Ini.HasKey("only_positive_pw") ? true : false);
        double gamW    = Ini("gamW");
        double z_water = (NDim==2 ? Ini("y_water") : Ini("z_water"));
        for (size_t i=0; i<GE->NN; ++i)
        {
            // elevation of node
            double z = (NDim==2 ? Con[i]->Vert.C(1) : Con[i]->Vert.C(2));

            // pore-water pressure
            double hw = z_water-z; // column of water
            double pw = (hw>0.0 ? gamW*hw : (pos_pw ? 0.0 : gamW*hw));

            // set U in node
            Con[i]->U("pw") = pw;
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * HydroMechElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new HydroMechElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int HydroMechElemRegister()
{
    ElementFactory["HydroMech"]   = HydroMechElemMaker;
    ElementVarKeys["HydroMech2D"] = std::make_pair ("ux uy pw",    "fx fy qw");
    ElementVarKeys["HydroMech3D"] = std::make_pair ("ux uy uz pw", "fx fy fz qw");
    PROB.Set ("HydroMech", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __HydroMechElem_dummy_int  = HydroMechElemRegister();

}; // namespace FEM

#endif // MECHSYS_HYDROMECH_H
