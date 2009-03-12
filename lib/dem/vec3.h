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

#ifndef __VEC3_H
#define __VEC3_H

/*
#define DO_INLINE_VEC3 1

#if DO_INLINE_VEC3 >= 1
#define  inline
#else
#define 
#endif
*/

#include <fstream>
#include <iostream>
#include <math.h>
#include <string>

// #include "Foundation/Error.h"


//=== included by Fernando Alonso=======
#include <vector>
using namespace std;
//class Vec3;

//=========================================

using std::ostream;
using std::istream;
using std::string;
using std::cout; //!added by F. Alonso
using std::endl; //!added by F. Alonso

class Matrix3;

// class VecErr:public MError
// {
//  public:
//   VecErr(const string&);
//   virtual ~VecErr(){};
// };

struct VDMulVadd;
struct VDMul;
/** Class Vec3
*/
class Vec3
{
protected:
  double data[3]; ///< vector components

public:
  static const Vec3 ZERO; //! The zero vector.
  static const Vec3 I;
  static const Vec3 J;
  static const Vec3 K;
  
  
  //! constructors
   Vec3();   //! create Vec3 V
   explicit Vec3(double s);  //! create Vec3 V(s)
   Vec3(double,double,double); //! create Vec3 V(a,b,c)
   Vec3(const Vec3&);  //! create Vec3 V(v)

  // vec-vec operators
   Vec3& operator=(const Vec3&);  //! V=0
   Vec3& operator=(double s);     //! V=s
   Vec3& operator-=(const Vec3&); //! V=U
   Vec3& operator+=(const Vec3&); //! V+=U
   Vec3 operator+(const Vec3&) const; //! V+U
   Vec3 operator-(const Vec3&) const; //! V-U

//   Vec3 operator*(const Matrix3 &m) const; //! V*M
   double operator*(const Vec3&) const;    //! s*V
   Vec3 operator-() const;                //!-V
  
  // vec-dbl ops
   Vec3 operator*(double) const;  //! V*s
   Vec3& operator*=(double);      //! V*=s
   Vec3 operator/(double) const;  //! V/s
   Vec3 operator-(double) const;  //! V-s
   Vec3 operator+(double) const;  //! V+s
   Vec3& operator+=(double);      //! V+=s
   Vec3& operator-=(double);      //! V-=s

   Vec3& operator/=(double);      //! V/=s
   double norm() const;           //! V.norm()=||V||
   double norm2() const;          //! V.norm2()=V*V
   Vec3 unit() const;             //! V.unit2()=V/||V||
   Vec3 unit_s() const; //!safe version of unit() (throw exceptions)
   double max() const;  //! maximal component of U
   double min() const;  //! minimal component of U

   Vec3 rotate(const Vec3 &axis, const Vec3 &axisPt) const;
   
   void show_on_screen(void) {cout <<data[0]<<"i+"<<data[1]<<"j+"<<data[2]<<"k "; };

   bool operator==(const Vec3&) const; //! V==U
   bool operator!=(const Vec3&) const; //! V!=U 

   friend Vec3 cmax(const Vec3&,const Vec3&);
   friend Vec3 cmin(const Vec3&,const Vec3&);

   friend Vec3 cross(const Vec3&,const Vec3&); //! cros(U,V)= U X V
   friend double dot(const Vec3&,const Vec3&); //! dot(U,V) = U*V
   friend Vec3 operator*(double,const Vec3&);  //! U*V

  //n+1-ary operators
   void mul_add_and_assign(const Vec3*,const Vec3*,const double&);
   void mul_and_assign(const Vec3*,const double&);

   Vec3(const VDMulVadd&);
   Vec3& operator=(const VDMulVadd&);

   Vec3(const VDMul&);
   Vec3& operator=(const VDMul&);


   void set_x(double x) {data[0] = x;}
   void set_y(double y) {data[1] = y;}
   void set_z(double z) {data[2] = z;}
//  void set_xyz(double x, double y, double z)
//  { data[0] = x; data[1] = y; data[2] = z;}
  
   double& X() {return data[0];};
   double& Y() {return data[1];};
   double& Z() {return data[2];};
   double X() const {return data[0];};
   double Y() const {return data[1];};
   double Z() const {return data[2];};

   const double &operator[](int i) const {return data[i];}
   double& operator[](int i) {return data[i];}

  //! in/output
   friend ostream& operator << (ostream&,const Vec3&);
   friend istream& operator >> (istream&,Vec3&);

  //! comparison -> enable to use of Vec3 as key in STL map and set
  bool operator<(const Vec3&) const; 

  friend class Matrix3;
};

 Vec3 comp_max(const Vec3&,const Vec3&); //!< per component maximum
 Vec3 comp_min(const Vec3&,const Vec3&); //!< per component minimum

// #if DO_INLINE_VEC3 >= 1
// #include "Foundation/vec3.hpp"
// #endif

#endif // __VEC3_H
