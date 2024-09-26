#define ASSERT(x) assert(x)

#include "FB_defs.h"

#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "stdio.h"
#include <vector>
#include <algorithm>
#include "../global.h"
#include "../defs.h"


#include "Common/binary_tree.h"
#include "Common/Float.h"
#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "ar_interaction.hpp"
#include "ar_perturber.hpp"


void FBTermination2(Particle* ptclCM, REAL current_time, std::vector<Particle*> &particle);


// made 2024.08.12 by Eunwoo Chung

// Reference: SDAR/sample/AR/ar.cxx & PeTar/src/hard.hpp
// No debugging yet

// true: integrated normally, false: terminated by stellar merger, TDE, GW merger, etc.
// If Intererrupt_mode != none, then bin_termination = true;
bool Group::ARIntegration(REAL next_time, std::vector<Particle*> &particle){

    auto bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);
    // ASSERT((next_time - StartTime) > 0);
    // auto bin_interrupt = sym_int.integrateToTime((next_time - StartTime)*EnzoTimeStep);

    sym_int.particles.shiftToOriginFrame();
    sym_int.particles.template writeBackMemberAll<Particle>(); // Eunwoo: I'm not sure
    // for (Particle* member : Members) {
    //     fprintf(binout, "PID: %d. posx: %e, posy: %e, posz: %e\n", member->PID, member->Position[0], member->Position[1], member->Position[2]);
    // }
    // sym_int.particles.writeBackMemberAll<Particle>();
    // for (Particle* member : Members) {
    //     fprintf(binout, "PID: %d. posx: %e, posy: %e, posz: %e\n", member->PID, member->Position[0], member->Position[1], member->Position[2]);
    // }

    if (bin_interrupt.status != AR::InterruptStatus::none) {

        FBTermination2(groupCM, bin_interrupt.time_now/EnzoTimeStep, particle);

        // bin_termination = true;

        return false;
    }

    sym_int.particles.shiftToCenterOfMassFrame();
    CurrentTime = next_time;
    // fprintf(binout, "AR done! time: %e\n", next_time*EnzoTimeStep*1e4);
    return true; 
}