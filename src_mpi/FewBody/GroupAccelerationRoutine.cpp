#include <iostream>
#include "../global.h"
#include "../def.h"
#include <cassert> // Eunwoo debug


void NewFBInitialization(int newOrder);
// void NewFBInitialization2(Particle* ptcl);
void MergeGroups();
void NewPrimordialBinaries(int newOrder);
void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void broadcastFromRoot(int &data);

/*
bool AddNewGroupsToList(std::vector<Particle*> &particle) {

	assert(GroupCandidateList.empty()); // Let's check GroupCandidateList is initially empty!

	for (Particle *ptcl : particle) {
        // if (ptcl->isCMptcl) continue; // Eunwoo: test
		// ptcl->isFBCandidate();

        // Eunwoo test
        if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 > tbin) // fiducial: 1e-5 but for rbin = 0.00025 pc, 1e-6 Myr seems good
            continue;
        // Eunwoo test
        
        ptcl->checkNewGroup();
	}

	if (GroupCandidateList.empty()) return true;

	MergeGroups(GroupCandidateList);	// Merge GroupCandidateList
										// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!

	for (Group *groupCandidate : GroupCandidateList) {

        // if (isNeighborInsideGroup(groupCandidate)) continue;    // Don't make a real group if it has no neighbor.
                                                                // Its TimeStepIrr will be so big and it can raise errors.

        NewFBInitialization(groupCandidate, particle);
				
	}
	GroupCandidateList.clear();
	GroupCandidateList.shrink_to_fit();

	return true;
}
*/
/*
bool AddNewGroupsToList2(std::vector<Particle*> &members, std::vector<Particle*> &particle) {

	assert(GroupCandidateList.empty()); // Let's check GroupCandidateList is initially empty!

	for (Particle *ptcl : members) {
        // if (ptcl->isCMptcl) continue; // Eunwoo: test
		// ptcl->isFBCandidate();
        ptcl->checkNewGroup2();
	}

	if (GroupCandidateList.empty()) return true;

	MergeGroups(GroupCandidateList);	// Merge GroupCandidateList
										// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!

	for (Group *groupCandidate : GroupCandidateList) {

        // if (isNeighborInsideGroup(groupCandidate)) continue;    // Don't make a real group if it has no neighbor.
                                                                // Its TimeStepIrr will be so big and it can raise errors.

        NewFBInitialization2(groupCandidate, particle);
				
	}
	GroupCandidateList.clear();
	GroupCandidateList.shrink_to_fit();

	return true;
}
*/

void SetPrimordialBinaries() {

	Particle* ptcl;
	Particle* NewCM;

	for (int i=0; i<NumberOfSingle; i++) {
		ptcl = &particles[i];
		if (ptcl->NewNumberOfNeighbor > 0) {
			NewCM = &particles[NumberOfParticle];
			NewCM->clear();
			for (int j = 0; j < ptcl->NewNumberOfNeighbor; j++) {
				NewCM->NewNeighbors[j] = ptcl->NewNeighbors[j];
			}
			NewCM->NewNeighbors[ptcl->NewNumberOfNeighbor] = ptcl->ParticleOrder;
			NewCM->NewNumberOfNeighbor = ptcl->NewNumberOfNeighbor + 1;

			NumberOfParticle++;
		}
	}

	if (NumberOfSingle == NumberOfParticle) return;

	MergeGroups();	// Merge group candidates
					// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!

	broadcastFromRoot(NumberOfParticle);

	for (int i=NumberOfSingle; i<NumberOfParticle; i++) {
		NewPrimordialBinaries(i);
	}
}

void SetBinaries(std::vector<int>& ParticleList) {

	Particle* ptcl;
	Particle* NewCM;

	int beforeNumberOfParticle = NumberOfParticle;

	for (int i=0; i<ParticleList.size(); i++) {
		ptcl = &particles[ParticleList[i]];
		if (ptcl->NewNumberOfNeighbor > 0) {
			NewCM = &particles[NumberOfParticle];
			NewCM->clear();
			for (int j = 0; j < ptcl->NewNumberOfNeighbor; j++) {
				NewCM->NewNeighbors[j] = ptcl->NewNeighbors[j];
			}
			NewCM->NewNeighbors[ptcl->NewNumberOfNeighbor] = ptcl->ParticleOrder;
			NewCM->NewNumberOfNeighbor = ptcl->NewNumberOfNeighbor + 1;

			NumberOfParticle++;
		}
	}

	if (beforeNumberOfParticle == NumberOfParticle) return;

	MergeGroups();	// Merge group candidates
					// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!

	for (int i=beforeNumberOfParticle; i<NumberOfParticle; i++) {
		NewFBInitialization(i);
	}
	/*
	for (int i=NumberOfSingle; i<NumberOfParticle; i++) {
		Particle* members = &particles[i];

		if (!members->GroupInfo) {
			int beforeOrder = particles[NumberOfParticle-1].ParticleOrder;
			int afterOrder = members->ParticleOrder;

			std::swap(particles[i], particles[NumberOfParticle-1]);
			particles[i].ParticleOrder = afterOrder;
			NumberOfParticle--;
			for (int j=0; j<NumberOfParticle; j++) {
				Particle* ptcl = &particles[j];
				
				std::replace(
					ptcl->Neighbors,
					ptcl->Neighbors + ptcl->NumberOfNeighbor,
					beforeOrder,
					afterOrder
				);
			}
		}
	}
	*/
}


void MergeGroups() {

    bool merged = true;

    while (merged) {
        merged = false;
		for (int i = NumberOfSingle; i < NumberOfParticle; i++) {

			Particle* currentCM = &particles[i];
			if (currentCM->NewNumberOfNeighbor == 0) continue; // Skip already deleted groups

			for (int j = i + 1; j < NumberOfParticle; j++) {

				Particle* otherCM = &particles[j];
				if (otherCM->NewNumberOfNeighbor == 0) continue; // Skip already deleted groups

                // Check if there's any common member between group1 and group2
                bool commonFound = false;
				for (int k=0; k < currentCM->NewNumberOfNeighbor; k++) {
					int member1 = currentCM->NewNeighbors[k];
					if (std::find(otherCM->NewNeighbors, otherCM->NewNeighbors + otherCM->NewNumberOfNeighbor, member1) != otherCM->NewNeighbors + otherCM->NewNumberOfNeighbor) {
                        commonFound = true;
                        break;
                    }
                }

                // If common members are found, merge group2 into group1
                if (commonFound) {
                    // Merge group2 into group1, avoiding duplicates
					for (int l=0; l < otherCM->NewNumberOfNeighbor; l++) {
						int member2 = otherCM->NewNeighbors[l];
						if (std::find(currentCM->NewNeighbors, currentCM->NewNeighbors + currentCM->NewNumberOfNeighbor, member2) == currentCM->NewNeighbors + currentCM->NewNumberOfNeighbor) {
							currentCM->NewNeighbors[currentCM->NewNumberOfNeighbor] = member2;
							currentCM->NewNumberOfNeighbor++;
                        }
                    }
					merged = true;

                    // Mark otherGroup for deletion after the loop
					particles[j] = particles[NumberOfParticle]; // Eunwoo: is this right?
					particles[NumberOfParticle].clear();
					NumberOfParticle--;
                }
            }
        }
    }
}

/*
bool isNeighborInsideGroup(Group* groupCandidate) {
    // Check if every neighbor of 'member' is also inside groupMembers
    for (Particle* member : groupCandidate->Members) {
        for (Particle* neighbor : member->ACList) {
            if (std::find(groupCandidate->Members.begin(), groupCandidate->Members.end(), neighbor) == groupCandidate->Members.end()) {
                // Found a neighbor that is not in groupMembers
                return false;
            }
        }
    }
    return true; // All neighbors are inside groupMembers
}
*/

// Make CM particles for primordial binaries
void NewPrimordialBinaries(int newOrder) {

	Particle* ptclCM;
	Group* ptclGroup;

	// Set ptclGroup members first; this will be very useful

	ptclGroup = new Group();

#ifdef SEVN
	ptclGroup->useSEVN = false;
	REAL dt_evolve_next = NUMERIC_FLOAT_MAX; // Myr
	for (Particle* members : group->Members) {
		if (members->star == nullptr || members->star->amiremnant())
			continue;
		else {
			ptclGroup->useSEVN = true;
			assert(members->EvolutionTime == 0.0);
			if (members->star->getp(Timestep::ID) < dt_evolve_next)
				dt_evolve_next = members->star->getp(Timestep::ID);
		}
	}
	if (ptclGroup->useSEVN) {
		ptclGroup->EvolutionTime = 0.0;
		ptclGroup->EvolutionTimeStep = dt_evolve_next;
		// fprintf(binout, "EvolutionTimeStep: %e Myr\n", ptclGroup->EvolutionTimeStep);
		// fflush(binout);
	}
#endif

	ptclCM = &particles[newOrder];
	ptclGroup->groupCM		= ptclCM;

	ptclCM->ParticleOrder = newOrder;
	ptclCM->PID = NewPID;
	NewPID++;
	ptclCM->isActive = true;

	// Let's link CM particle with the cm particles made in the binary tree (SDAR).

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(ptclCM->NewNumberOfNeighbor); // Binary tree is made and CM particle is made automatically.

	// ptclCM = &ptclGroup->sym_int.particles.cm;
	for (int dim=0; dim<Dim; dim++) {
		ptclCM->Position[dim] = ptclGroup->sym_int.particles.cm.Position[dim];
		ptclCM->Velocity[dim] = ptclGroup->sym_int.particles.cm.Velocity[dim];
	}

	// Set ptcl information like time, PID, etc.

	ptclCM->GroupInfo		= ptclGroup;

	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);


	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];

		fprintf(binout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(binout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(binout, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
		fprintf(binout, "PID: %d. Particle type: %d\n", members->PID, members->ParticleType);
        fprintf(binout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(binout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(binout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(binout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(binout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);

		// fprintf(binout, "PID: %d. Position - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0], members->Position[1], members->Position[2]);
		// fprintf(binout, "PID: %d. Velocity - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0], members->Velocity[1], members->Velocity[2]);
		// fprintf(binout, "PID: %d. Mass - %e, \n", members->PID, members->Mass);
		// fprintf(binout, "PID: %d. Particle type: %d\n", members->PID, members->ParticleType);
        // fprintf(binout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(binout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(binout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
    }

	ptclGroup->sym_int.initialIntegration(0); // This is primordial binary!
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);


	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	CalculateAcceleration01(ptclCM);
	CalculateAcceleration23(ptclCM);


// /* // Eunwoo test
	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->RadiusOfNeighbor));
		fprintf(binout, "Bound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(binout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
		fprintf(binout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(binout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(binout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
	else {
		ptclGroup->sym_int.info.r_break_crit = 2*bin_root.semi*(1-bin_root.ecc); // r_break_crit = 2*peri
		fprintf(binout, "Unbound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(binout, "period: %e Myr\n", bin_root.period*1e4);
		fprintf(binout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(binout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
// */ // Eunwoo test

	// Erase members in neighbors

	for (int i=0; i<newOrder; i++) {
		Particle* ptcl = &particles[i];
		auto newEnd = std::remove_if(
			ptcl->Neighbors, 
			ptcl->Neighbors + ptcl->NumberOfNeighbor, 
			[&particles](int j) {
				return !particles[j].isActive;
			}
		);

		if (newEnd != ptcl->Neighbors + ptcl->NumberOfNeighbor) {
			ptcl->NumberOfNeighbor = newEnd - ptcl->Neighbors;
			ptcl->Neighbors[ptcl->NumberOfNeighbor] = ptclCM->ParticleOrder;
			ptcl->NumberOfNeighbor++;
		}
	}

	fprintf(binout, "\nResult of CM particle value calculation from function NewPrimordialBinaries\n");

	fprintf(binout, "Position (pc) - x:%e, y:%e, z:%e, \n", ptclCM->Position[0]*position_unit, ptclCM->Position[1]*position_unit, ptclCM->Position[2]*position_unit);
	fprintf(binout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[1]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[2]*velocity_unit/yr*pc/1e5);
	fprintf(binout, "Mass (Msol) - %e, \n", ptclCM->Mass*mass_unit);
	fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	fprintf(binout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_reg[0][0], ptclCM->a_reg[1][0], ptclCM->a_reg[2][0]);
	// fprintf(binout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_reg[0][1], ptclCM->a_reg[1][1], ptclCM->a_reg[2][1]);
	// fprintf(binout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_reg[0][2], ptclCM->a_reg[1][2], ptclCM->a_reg[2][2]);
	// fprintf(binout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_reg[0][3], ptclCM->a_reg[1][3], ptclCM->a_reg[2][3]);
	fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_irr[0][0], ptclCM->a_irr[1][0], ptclCM->a_irr[2][0]);
	// fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_irr[0][1], ptclCM->a_irr[1][1], ptclCM->a_irr[2][1]);
	// fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_irr[0][2], ptclCM->a_irr[1][2], ptclCM->a_irr[2][2]);
	// fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_irr[0][3], ptclCM->a_irr[1][3], ptclCM->a_irr[2][3]);

	fprintf(binout, "------------------END-OF-NEW-PRIMORDIAL-BINARIES------------------\n\n");
	fflush(binout);
}