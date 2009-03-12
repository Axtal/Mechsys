#ifndef __QUATERNION_H
#define __QUATERNION_H

#include "vec3.h"
/** Class defining quaternion entities for tri-dimensional rotations. */
class Quaternion {
    protected:
        double 	q0; 	///< Scalar part of the quaternion representing the angle of rotation
        Vec3 	q; 	///< Vectorial part of the quaternion representing the axis of rotation
    public:

	// Constructors
        Quaternion() {}; 											///< Default constructor, does nothing
        Quaternion(double a,Vec3 v) {q0=a; q=v;}; 								///< Constructor in which you include the scalar and vectorial part
        Quaternion(double a,double b, double c,double d) {q0=a; q.set_x(b); q.set_y(c); q.set_z(d);};		///< Constructor requiring 4 scalars

	// Methods
        double norm(void);											///< Gives the euclidean norm of the quaternion
        void normalize(void);											///< Normalize the quaternion to unit norm
        friend Quaternion normalize_rotation(double a,Vec3 & v);						///< Gives a quaternion given by the angle of rotation in radians around the axis vector v
        Quaternion conj(void);											///< Conjugates the quaternion (conj(q0+q)=q0-q);
        Quaternion operator *(Quaternion Q);									///< Quaternion product
        Quaternion operator +(Quaternion Q);									///< Quaternion addition
        Quaternion operator -(Quaternion Q);									///< Quaternion substraction
        friend Vec3 rotate(Quaternion & Q,Vec3 & V);								///< Gives a vector as a rotated vector by a Quaternion
        Quaternion operator *(double temp);
        
	//Access methods
	Vec3 getvector(void) {return q;};									///< Return the vectorial part of the quaternion
        double getscalar(void) {return q0;};									///< Return the scalar part of the quaternion
        void show_on_screen(void) {cout <<q0<<" "; q.show_on_screen(); cout <<endl;};
};



#endif
