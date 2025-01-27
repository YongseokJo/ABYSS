/*
 Created by iorio on 7/26/22.
 it includes:
- Class Zeros: Kicks velocity always=0
*/

#ifndef SEVN_ZEROS_H
#define SEVN_ZEROS_H

#include <kicks.h>

/**
 * Kicks velocity always=0
 */
class Zeros : public Kicks{

public:
    Zeros(bool reg = true){
        if(reg) {
            Register(this, name());
        }
    }

    static Zeros _zeros;

    void _apply(Star *s) override;

    Zeros* instance() {
        return (new Zeros(false));
    }

    inline std::string name() override { return "zeros"; }

};

#endif //SEVN_ZEROS_H
