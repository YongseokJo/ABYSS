//
// Created by mario on 07/11/18.
//

#ifndef SEVN_REVISED_TYPES_H
#define SEVN_REVISED_TYPES_H

#include <iostream>
#include <cmath>
#include <vector>

struct double3 {

    double x;
    double y;
    double z;

    inline bool isnan() const{
        if(std::isnan(x) || std::isnan(y) || std::isnan(z)) return true;
        else return false;
    }

    inline bool isinf() const{
        if(std::isinf(x) || std::isinf(y) || std::isinf(z)) return true;
        else return false;
    }

};


struct double4{

public:

    double4(){ x = y = z = w = 0.0;}
    explicit double4(const double3 *a) {x = a->x; y = a->y; z = a->z; w = mod();}


    inline void set_x(double a) {x = a; mod();}
    inline void set_y(double a) {y = a; mod();}
    inline void set_z(double a) {z = a; mod();}
    inline void set_w(double a) {w = a;}


private:
    double x;
    double y;
    double z;
    double w;

    inline double mod(){
        return std::sqrt(x*x + y*y + z*z);
    }

};

struct starprint{
    double m;
    double mhe;
    double rco;
};


#endif //SEVN_REVISED_TYPES_H
