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

#ifndef MECHSYS_MECHSYS_H
#define MECHSYS_MECHSYS_H

// MechSys -- mesh
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "mesh/unstructured.h"
#include "mesh/alphashape.h"

// MechSys -- fem -- basic
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/output.h"
#include "models/model.h"
#include "util/exception.h"

// MechSys -- fem -- Elements
#include "fem/elems/lin2.h"
#include "fem/elems/tri3.h"
#include "fem/elems/tri6.h"
#include "fem/elems/quad4.h"
#include "fem/elems/quad8.h"
#include "fem/elems/hex8.h"
#include "fem/elems/hex20.h"
#include "fem/elems/rod.h"
#include "fem/elems/beam.h"
#include "fem/elems/spring.h"
#include "fem/diffusionelem.h"
#include "fem/equilibelem.h"
#include "fem/biotelem.h"

// MechSys -- fem -- Embedded
//#include "fem/embedded.h"
//#include "fem/elems/rod3.h"
//#include "fem/elems/embspring.h"

// MechSys -- Models
#include "models/equilibs/beamelastic.h"
#include "models/equilibs/biotelastic.h"
#include "models/equilibs/camclay.h"
#include "models/equilibs/linelastic.h"
#include "models/equilibs/rodelastic.h"
#include "models/equilibs/springelastic.h"
#include "models/equilibs/pyequilib.h"

#endif // MECHSYS_MECHSYS_H
