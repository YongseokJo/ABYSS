//
// Created by mario on 03/05/18.
//

#ifndef SEVN_REVISED_STARPARAMETER_H
#define SEVN_REVISED_STARPARAMETER_H

#include <string>
#include <vector>
#include <map>

//TODO Clean and migrate everything to params.h
namespace starparameter{

    extern const double maximum_variation;
    extern const double wrtolerance;
    extern const size_t min_points_per_phase;

    extern const std::string semimajor;

    extern double MAX_ZAMS, MIN_ZAMS, MAX_Z, MIN_Z;


    enum PROC_ID {CE = 0, RL, TD, COUNT};
}

#endif //SEVN_REVISED_STARPARAMETER_H
