#ifndef Vect_H
#define Vect_H

#include <math.h>

class Vect {
    double x, y, z;

    public:

    Vect();

    Vect(double, double, double);

    // method functions

    double getVX() { return x; }
    double getVY() { return y; }
    double getVZ() { return z; }

    double magnitude() { return sqrt(x*x + y*y +z*z); }

    double squaredMagnitude() { return (x*x + y*y +z*z); }

    Vect normalize() {
        double Vmag = this->magnitude();
        return Vect(x/ Vmag, y/ Vmag, z/ Vmag);
    }

    Vect negative () {
        return Vect (-x, -y, -z);
    }

    double dotProduct(Vect v) {
        return x*v.getVX() + y*v.getVY() + z*v.getVZ();
    }

    Vect crossProduct(Vect v) {
        return Vect (y*v.getVZ() - z*v.getVY(), z*v.getVX() - x*v.getVZ(), x*v.getVY() - y*v.getVX());
    }

    Vect vectAdd (Vect v) {
        return Vect (x + v.getVX(), y + v.getVY(), z + v.getVZ());
    }

    Vect vectMult (double scalar) {
        return Vect (x*scalar, y*scalar, z*scalar);
    }

    Vect operator+(Vect v_1){
        return Vect (x + v_1.getVX(), y + v_1.getVY(), z + v_1.getVZ());
    }

//    Vect operator=(Vect){
//        return Vect (-x, -y, -z);
//    }

    Vect operator-(Vect v_1){
        return Vect (x - v_1.getVX(), y - v_1.getVY(), z - v_1.getVZ());
    }

    Vect operator*(double scalar){
        return Vect (x*scalar, y*scalar, z*scalar);
    }

    Vect operator/(double scalar){
        return Vect (x/scalar, y/scalar, z/scalar);
    }

};

Vect::Vect () {
    x = 0;
    y = 0;
    z = 0;
}

Vect::Vect (double i, double j, double k) {
    x = i;
    y = j;
    z = k;
}
#endif // Vect_H
