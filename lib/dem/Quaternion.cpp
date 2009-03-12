#include "Quaternion.h"
#include <math.h>

double Quaternion::norm() {
    return sqrt(q0*q0+dot(q,q));
}

void Quaternion::normalize(void) {
    double n=norm();
    q0=q0/n;
    q=q/n;
}

Quaternion normalize_rotation(double theta,Vec3 &axis) {
    Vec3 n=axis.unit_s();
    Quaternion r;
    r.q0=cos(theta/2.);
    r.q=n*(sin(theta/2.));
    return r;
}

Quaternion Quaternion::conj(void) {
    return Quaternion(q0,-q);
}

Quaternion Quaternion::operator * (Quaternion temp) {
    Quaternion r;
    r.q0=q0*temp.q0-dot(q,temp.q);
    r.q=q0*temp.q+temp.q0*q+cross(q,temp.q);
    return r;
}

Quaternion Quaternion::operator + (Quaternion temp) {
    Quaternion r;
    r.q0=q0+temp.q0;
    r.q=q+temp.q;
    return r;
}

Quaternion Quaternion::operator - (Quaternion temp) {
    Quaternion r;
    r.q0=q0-temp.q0;
    r.q=q-temp.q;
    return r;
}

Quaternion Quaternion::operator*(double temp) {
    Quaternion r;
    r.q0=q0*temp;
    r.q=q*temp;
    return r;
    
}

Vec3 rotate(Quaternion & Q,Vec3 & V) {
    return (Q*(Quaternion(0,V)*Q.conj())).getvector();
}

