#include "vec3.h"



const Vec3 Vec3::ZERO = Vec3(0.0, 0.0, 0.0);
const Vec3 Vec3::I = Vec3(1, 0.0, 0.0);
const Vec3 Vec3::J = Vec3(0.0, 1, 0.0);
const Vec3 Vec3::K = Vec3(0.0, 0.0, 1);
// commented by F. Alonso
// #include "Foundation/Matrix3.h"

//the error...
//  VecErr::VecErr(const string& m):MError(m)
// {
//   message.insert(0,"Vec3 "); 
// }

// constructors
 Vec3::Vec3()
{
  data[0]=0;
  data[1]=0;
  data[2]=0;
}

 Vec3::Vec3(double s)
{
  data[0]=s;
  data[1]=s;
  data[2]=s;
}

// anexed by fernando Alonso
// Vec3::Vec3(double a,double b)
//{
//  data[0]=a;
//  data[1]=b;
//  data[2]=0;
//}

 Vec3::Vec3(double a,double b,double c)
{
  data[0]=a;
  data[1]=b;
  data[2]=c;
}

 Vec3::Vec3(const Vec3& rhs)
{
  data[0]=rhs.data[0];
  data[1]=rhs.data[1];
  data[2]=rhs.data[2];
}

// operators

 Vec3& Vec3::operator=(const Vec3& rhs)
{
  data[0]=rhs.data[0];
  data[1]=rhs.data[1];
  data[2]=rhs.data[2];
  return *this;
}

 Vec3& Vec3::operator=(double s)
{
  data[0]=s;
  data[1]=s;
  data[2]=s;
  return *this;
}

 Vec3& Vec3::operator-=(const Vec3& rhs)
{
  data[0]-=rhs.data[0];
  data[1]-=rhs.data[1];
  data[2]-=rhs.data[2];
  return *this;
}

 Vec3& Vec3::operator+=(const Vec3& rhs)
{
  data[0]+=rhs.data[0];
  data[1]+=rhs.data[1];
  data[2]+=rhs.data[2];
  return *this;
}

 Vec3 Vec3::operator+(const Vec3& rhs) const
{
  return Vec3(data[0]+rhs.data[0], data[1]+rhs.data[1], data[2]+rhs.data[2]);
}

 Vec3 Vec3::operator-(const Vec3& rhs) const
{
  return Vec3(data[0]-rhs.data[0], data[1]-rhs.data[1], data[2]-rhs.data[2]);
}

 Vec3 Vec3::operator-() const
{
  return Vec3( -data[0],-data[1],-data[2] );
}

/*
commented by F. Alonso
 Vec3 Vec3::operator*(const Matrix3 &m) const
{
  const double x = m(0,0)*data[0] + m(1,0)*data[1] + m(2,0)*data[2];
  const double y = m(0,1)*data[0] + m(1,1)*data[1] + m(2,1)*data[2];
  const double z = m(0,2)*data[0] + m(1,2)*data[1] + m(2,2)*data[2];

  return Vec3(x,y,z);
}*/

 double Vec3::operator*(const Vec3& rhs) const
{
  return data[0]*rhs.data[0]+data[1]*rhs.data[1]+data[2]*rhs.data[2];
}

 Vec3 Vec3::operator*(double s) const 
{ 
   return Vec3(data[0]*s,data[1]*s,data[2]*s) ; 
} 

 Vec3& Vec3::operator*=(double rhs)
{
  data[0]*=rhs;
  data[1]*=rhs;
  data[2]*=rhs;
  return *this;
}

 Vec3& Vec3::operator/=(double c)
{
  data[0] /= c;
  data[1] /= c;
  data[2] /= c;

  return *this;
}

 Vec3 Vec3::operator/(double s) const 
{ 
   return Vec3(data[0]/s,data[1]/s,data[2]/s) ; 
}

 Vec3 Vec3::operator+(double s) const 
{
   return Vec3(data[0]+s, data[1]+s, data[2]+s) ; 
}

 Vec3 Vec3::operator-(double s) const 
{
   return Vec3(data[0]-s, data[1]-s, data[2]-s); 
}

 Vec3 Vec3::rotate(const Vec3 &axis, const Vec3 &axisPt) const
{
  const double phi = axis.norm();
  if (phi > 0.0)
  {
    const Vec3 r = *this - axisPt;
    const Vec3 n = axis/phi;
    const double cosPhi = cos(phi);
    const Vec3 rotatedR =
      r*cosPhi + n*((dot(n, r))*(1-cosPhi)) + cross(r, n)*sin(phi);
    return rotatedR + axisPt;
  }
  return *this;
}

 Vec3 &Vec3::operator+=(double s)
{
   data[0] += s;
   data[1] += s;
   data[2] += s;
   return *this; 
}

 Vec3 &Vec3::operator-=(double s)
{
   data[0] -= s;
   data[1] -= s;
   data[2] -= s;
   return *this; 
}

// included by Fernando Alonso
// cross product s=cross(v1,v2)==v1 x v2
//  double cross2d(const Vec3& v1, const Vec3& v2)
//{
//  return v1.data[0] * v2.data[1] -  v1.data[1] * v2.data[0];
//}

// vector product
// 9 Flops ( 6 mult, 3 sub ) 
 Vec3 cross(const Vec3& lhs,const Vec3& rhs)
{
  return Vec3(lhs.data[1]*rhs.data[2]-lhs.data[2]*rhs.data[1],
              lhs.data[2]*rhs.data[0]-lhs.data[0]*rhs.data[2],
              lhs.data[0]*rhs.data[1]-lhs.data[1]*rhs.data[0]);
}

//  dot product  

 double dot(const Vec3& v1, const Vec3& v2)
{
  return v1.data[0] * v2.data[0] + 
         v1.data[1] * v2.data[1] + 
         v1.data[2] * v2.data[2];
}

 Vec3 operator*(double f,const Vec3& rhs)
{
  return Vec3(f*rhs.data[0], f*rhs.data[1], f*rhs.data[2]);
}


// euclidian norm
// 6 Flops ( 3 mult, 2 add, 1 sqrt )
 double Vec3::norm() const
{
  return sqrt(data[0]*data[0]+data[1]*data[1]+data[2]*data[2]);
}

// square of the euclidian norm
// 5 Flops ( 3 mult, 2 add)
 double Vec3::norm2() const
{
  return data[0]*data[0]+data[1]*data[1]+data[2]*data[2];
}

// returns unit vector in direction of the original vector
// 9 Flops ( 3 mult, 2 add, 3 div, 1 sqrt ) 
 Vec3 Vec3::unit() const
{
  return (*this)/norm();
}

// per element min/max
 Vec3 cmax(const Vec3& v1,const Vec3& v2)
{
  Vec3 res;
  res.data[0]=v1.data[0]>v2.data[0] ? v1.data[0] : v2.data[0];
  res.data[1]=v1.data[1]>v2.data[1] ? v1.data[1] : v2.data[1];
  res.data[2]=v1.data[2]>v2.data[2] ? v1.data[2] : v2.data[2];
  return res;
}

 Vec3 cmin(const Vec3& v1,const Vec3& v2)
{
  Vec3 res;
  res.data[0]=v1.data[0]<v2.data[0] ? v1.data[0] : v2.data[0];
  res.data[1]=v1.data[1]<v2.data[1] ? v1.data[1] : v2.data[1];
  res.data[2]=v1.data[2]<v2.data[2] ? v1.data[2] : v2.data[2];
  return res;
}

// save version, throws exception if norm()==0
 Vec3 Vec3::unit_s() const
{
  double n=norm();
  if(n==0) cout<<"norm() of data[2]ero-vector"; 
  Vec3 res(data[0],data[1],data[2]);
  return res/n;
}

 double Vec3::max() const
{
  double m = ( data[0]>data[1] ? data[0] : data[1] );

  return ( m>data[2] ? m : data[2] );

}

 double Vec3::min() const
{
  double m = ( data[0]<data[1] ? data[0] : data[1] );

  return ( m<data[2] ? m : data[2] );
}

 bool Vec3::operator==(const Vec3& V) const
{
  return((data[0]==V.data[0])&&(data[1]==V.data[1])&&(data[2]==V.data[2]));
}

 bool Vec3::operator!=(const Vec3& V) const
{
  return((data[0]!=V.data[0])||(data[1]!=V.data[1])||(data[2]!=V.data[2]));
}

// per component min/max

 Vec3 comp_max(const Vec3& V1,const Vec3& V2)
{
  double x=(V1.X() > V2.X()) ? V1.X() : V2.X();
  double y=(V1.Y() > V2.Y()) ? V1.Y() : V2.Y();
  double z=(V1.Z() > V2.Z()) ? V1.Z() : V2.Z();
 
  return Vec3(x,y,z);
}

 Vec3 comp_min(const Vec3& V1,const Vec3& V2)
{
  double x=(V1.X() < V2.X()) ? V1.X() : V2.X();
  double y=(V1.Y() < V2.Y()) ? V1.Y() : V2.Y();
  double z=(V1.Z() < V2.Z()) ? V1.Z() : V2.Z();
 
  return Vec3(x,y,z);
}

// in/output

 std::ostream& operator << (std::ostream& ostr,const Vec3& V)
{
  const char delimiter = ' ';
  ostr
    << V.data[0] << delimiter
    << V.data[1] << delimiter
    << V.data[2];

  return ostr;
}

 std::istream& operator >> (std::istream& istr,Vec3& V)
{
  istr
    >> V.data[0]
    >> V.data[1]
    >> V.data[2];

  return istr;
}


bool Vec3::operator<(const Vec3& rhs) const
{
  bool res;

  if(data[0]!=rhs.data[0]) {
    res=data[0]<rhs.data[0];
  } else if(data[1]!=rhs.data[1]){
    res=data[1]<rhs.data[1];
  } else if(data[2]!=rhs.data[2]){
    res=data[2]<rhs.data[2];
  } else {
    res=false;
  }

  return res;
}
