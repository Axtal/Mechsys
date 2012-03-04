/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso                                *
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

/* ColorMaps - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_COLORMAP_H
#define MPM_COLORMAP_H

// FLTK
#include <FL/Enumerations.H> // for Fl_Color
#include <FL/Fl_Choice.H>

// MechSys
#include <mechsys/util/colors.h>

/** Color Maps. */
namespace ClrMap
{

Fl_Menu_Item Items[] =
{
	{"BW/Paper" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Hot"      , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Jet"      , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{0,0,0,0,0,0,0,0,0}
};

/** Colour. */
struct Colour
{
	double R;
	double G;
	double B;
}; // struct Colour

// Color maps
const size_t BWP_NCOLORS = 9;
const Colour BWP_RGB[]=
{{ 0.0  , 0.0  , 0.0  },
 { 0.15 , 0.15 , 0.5  },
 { 0.3  , 0.15 , 0.75 },
 { 0.6  , 0.2  , 0.5  },
 { 1.0  , 0.25 , 0.15 },
 { 0.9  , 0.5  , 0.0  },
 { 0.9  , 0.75 , 0.1  },
 { 0.9  , 0.9  , 0.5  },
 { 1.0  , 1.0  , 1.0  }};
const size_t HOT_NCOLORS = 12;
const Colour HOT_RGB[]=
{{ 0.25 , 0.0  , 0.0  },
 { 0.5  , 0.0  , 0.0  },
 { 0.75 , 0.0  , 0.0  },
 { 1.0  , 0.0  , 0.0  },
 { 1.0  , 0.25 , 0.0  },
 { 1.0  , 0.5  , 0.0  },
 { 1.0  , 0.75 , 0.0  },
 { 1.0  , 1.0  , 0.0  },
 { 1.0  , 1.0  , 0.25 },
 { 1.0  , 1.0  , 0.5  },
 { 1.0  , 1.0  , 0.75 },
 { 1.0  , 1.0  , 1.0  }};
const size_t JET_NCOLORS = 12;
const Colour JET_RGB[]=
{{ 0.0     , 0.0     , 2.0/3.0 },
 { 0.0     , 0.0     , 1.0     },
 { 0.0     , 1.0/3.0 , 1.0     },
 { 0.0     , 2.0/3.0 , 1.0     },
 { 0.0     , 1.0     , 1.0     },
 { 1.0/3.0 , 1.0     , 2.0/3.0 },
 { 2.0/3.0 , 1.0     , 1.0/3.0 },
 { 1.0     , 1.0     , 0.0     },
 { 1.0     , 2.0/3.0 , 0.0     },
 { 1.0     , 1.0/3.0 , 0.0     },
 { 1.0     , 0.0     , 0.0     },
 { 2.0/3.0 , 0.0     , 0.0     }};

static inline Fl_Color Color2FlColor (Colour const & C)
{
	uchar r = static_cast<uchar>(C.R*255.0);
	uchar g = static_cast<uchar>(C.G*255.0);
	uchar b = static_cast<uchar>(C.B*255.0);
	return fl_rgb_color (r,g,b);
}

Fl_Color BWP [BWP_NCOLORS];
Fl_Color HOT [HOT_NCOLORS];
Fl_Color JET [JET_NCOLORS];

static inline int __initialize_colormaps ()
{
	for (size_t i=0; i<BWP_NCOLORS; ++i) BWP[i] = Color2FlColor (BWP_RGB[i]);
	for (size_t i=0; i<HOT_NCOLORS; ++i) HOT[i] = Color2FlColor (HOT_RGB[i]);
	for (size_t i=0; i<JET_NCOLORS; ++i) JET[i] = Color2FlColor (JET_RGB[i]);
	return 0;
}

int __dummy_initialize_colormaps = __initialize_colormaps ();



// Auxiliary functions
Array<Fl_Color> __colorset;

int __init_colorset()
{
    Vec3_t const * c = &Colors::Get("yellow_light");
	uchar r = static_cast<uchar>((*c)(0)*255.0);
	uchar g = static_cast<uchar>((*c)(1)*255.0);
	uchar b = static_cast<uchar>((*c)(2)*255.0);
	__colorset.Push (fl_rgb_color(r,g,b));

    c = &Colors::Get("melon");
	r = static_cast<uchar>((*c)(0)*255.0);
	g = static_cast<uchar>((*c)(1)*255.0);
	b = static_cast<uchar>((*c)(2)*255.0);
	__colorset.Push (fl_rgb_color(r,g,b));

    c = &Colors::Get("peacock");
	r = static_cast<uchar>((*c)(0)*255.0);
	g = static_cast<uchar>((*c)(1)*255.0);
	b = static_cast<uchar>((*c)(2)*255.0);
	__colorset.Push (fl_rgb_color(r,g,b));

    return 0;
}

int __dummy_init_colorset = __init_colorset();

Fl_Color GetColor (int Tag)
{
    if      (Tag==-1) return __colorset[0];
    else if (Tag==-2) return __colorset[1];
    else              return __colorset[2];
}

}; // namespace ClrMap

#endif // MPM_COLORMAP_H
