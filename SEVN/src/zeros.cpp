//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <zeros.h>
#include <star.h>

void Zeros::_apply(Star *s) {
    random_velocity_kick = s->vkick[0] = s->vkick[1] = s->vkick[2] = s->vkick[3] = 0.;
}
