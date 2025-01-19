//
// Created by Giuliano Iorio on 14/12/2021.
//

#ifndef SEVN_LAMBDA_BASE_H
#define SEVN_LAMBDA_BASE_H

#include <utilities.h>

class Star;
class IO;

/**
 * Base abstract class for Lambda implementation.
 *
 * The derived classes have to define their own public operator(Star *star)
 *
 */
class Lambda_Base{

public:

    Lambda_Base(){};

    Lambda_Base(_UNUSED const Star *s) {};

    Lambda_Base(_UNUSED const IO *io) {};

    virtual ~Lambda_Base() = default;

    virtual double operator()(const Star *s) = 0;
};


#endif //SEVN_LAMBDA_BASE_H
