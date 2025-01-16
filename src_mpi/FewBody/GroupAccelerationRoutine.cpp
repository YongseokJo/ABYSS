#ifdef FEWBODY
#include <iostream>
#include <map>
#include <unordered_map>
#include "../global.h"
#include "../def.h"
#include <cassert> // Eunwoo debug


void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void deleteNeighbors(int newOrder);
void mergeGroupCandidates(int OriginalLastParticleIndex);

void formPrimordialBinaries(int OriginalLastParticleIndex) {

	Particle* ptcl;
	Particle* NewCM;

	for (int i=0; i<=OriginalLastParticleIndex; i++) {
		ptcl = &particles[i];
		if (ptcl->NewNumberOfNeighbor > 0) {
			NewCM = &particles[LastParticleIndex+1];
			NewCM->clear();
			NewCM->copyNewNeighbor(ptcl);		
			NewCM->NewNeighbors[ptcl->NewNumberOfNeighbor] = ptcl->ParticleIndex;
			NewCM->NewNumberOfNeighbor++;

			LastParticleIndex++;
		}
	}

	if (OriginalLastParticleIndex == LastParticleIndex) return;

	mergeGroupCandidates(OriginalLastParticleIndex);	// Merge group candidates
					// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!
	
	for (int i=OriginalLastParticleIndex+1; i<=LastParticleIndex; i++) {
		deleteNeighbors(i);
	}
	global_variable->LastParticleIndex = LastParticleIndex;
}

void formBinaries(std::vector<int>& ParticleList, std::vector<int>& newCMptcls,
 std::unordered_map<int,int>& existing, std::unordered_map<int,int>& terminated) {

	Particle* ptcl;
	Particle* NewCM;

	int OriginalLastParticleIndex = LastParticleIndex;

	for (int i=0; i<ParticleList.size(); i++) {
		ptcl = &particles[ParticleList[i]];
		if (ptcl->NewNumberOfNeighbor > 0) {
			fprintf(stdout, "GAR. Num: %d\n", ptcl->NewNumberOfNeighbor + 1);
			fprintf(stdout, "GAR. PID: %d\n", ptcl->PID);
			NewCM = &particles[LastParticleIndex+1];
			NewCM->clear();
			NewCM->copyNewNeighbor(ptcl);
			for (int j = 0; j < ptcl->NewNumberOfNeighbor; j++) {
				fprintf(stdout, "GAR. PID: %d\n", particles[ptcl->NewNeighbors[j]].PID);
			}
			fflush(stdout);
			NewCM->NewNeighbors[ptcl->NewNumberOfNeighbor] = ptcl->ParticleIndex;
			NewCM->NewNumberOfNeighbor++;

			LastParticleIndex++;
		}
	}

	if (OriginalLastParticleIndex == LastParticleIndex) return;

	// fprintf(stdout, "A. global_variable->NOP: %d, NOP: %d\n", global_variable->LastParticleIndex, LastParticleIndex);

	mergeGroupCandidates(OriginalLastParticleIndex);	// Merge group candidates
					// ex) A & B form a group and B & C form a group --> Merge so that A & B & C become one group!
	
	// fprintf(stdout, "B. global_variable->NOP: %d, NOP: %d\n", global_variable->LastParticleIndex, LastParticleIndex);

	assert(LastParticleIndex > OriginalLastParticleIndex); 

	while (terminated.size() != 0) {

		Particle* ptcl = &particles[LastParticleIndex];
		auto it = terminated.begin();
		particles[it->first].copyNewNeighbor(ptcl);

		existing.insert({it->first, it->second});
		newCMptcls.push_back(it->first);
		terminated.erase(it);
		LastParticleIndex--;
		if (OriginalLastParticleIndex == LastParticleIndex)
			break;
	}
	if (OriginalLastParticleIndex != LastParticleIndex) {
		for (int i = OriginalLastParticleIndex+1; i <= LastParticleIndex; i++) {
			existing.insert({i, existing.size() % NumberOfWorker + 1});
			newCMptcls.push_back(i);
		}
	}

	for (int i: newCMptcls) {
		deleteNeighbors(i);
	}
	global_variable->LastParticleIndex = LastParticleIndex;

	ParticleList.insert(ParticleList.end(), newCMptcls.begin(), newCMptcls.end());
}


void mergeGroupCandidates(int OriginalLastParticleIndex) {

    bool merged = true;

    while (merged) {
        merged = false;
		for (int i = OriginalLastParticleIndex+1; i <= LastParticleIndex; i++) {

			Particle* currentCM = &particles[i];
			if (currentCM->NewNumberOfNeighbor == 0) continue; // Skip already deleted groups

			for (int j = i + 1; j <= LastParticleIndex; j++) {

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
					if (j != LastParticleIndex) {
						particles[j].copyNewNeighbor(&particles[LastParticleIndex]);
					}
					particles[LastParticleIndex].clear();
					LastParticleIndex--;
                }
            }
        }
    }
}


void deleteNeighbors(int newOrder) {

	Particle* ptclCM;

	std::cout << "New particle index is " << newOrder << std::endl;
	ptclCM = &particles[newOrder];
	ptclCM->ParticleIndex = newOrder;
	ptclCM->PID = NewPID;
	NewPID++;
	ptclCM->isActive = true;
	ptclCM->isCMptcl = true;
	ptclCM->RadiusOfNeighbor = ACRadius*ACRadius;
	ptclCM->binary_state = 0;

	Particle* members;

	for (int i=0; i<ptclCM->NewNumberOfNeighbor; i++) {
		members = &particles[ptclCM->NewNeighbors[i]];
		members->isActive = false;
	}
	NumberOfParticle += 1 - ptclCM->NewNumberOfNeighbor;

	// Erase members and put CM particle in neighbors
	for (int i=0; i<=LastParticleIndex; i++) {

		Particle* ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;

		auto newEnd = std::remove_if(
			ptcl->Neighbors, 
			ptcl->Neighbors + ptcl->NumberOfNeighbor, 
			[&particles](int j) {
				return !particles[j].isActive;
			}
		);

		if (newEnd != ptcl->Neighbors + ptcl->NumberOfNeighbor) {
			ptcl->NumberOfNeighbor = newEnd - ptcl->Neighbors;
			ptcl->Neighbors[ptcl->NumberOfNeighbor] = ptclCM->ParticleIndex;
			ptcl->NumberOfNeighbor++;
		}
	}
}

void makePrimordialGroup(Particle* ptclCM) {

	Group* ptclGroup = new Group();

	ptclGroup->groupCM = ptclCM;

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

	// Let's link CM particle with the cm particles made in the binary tree (SDAR).

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(ptclCM->NewNumberOfNeighbor); // Binary tree is made and CM particle is made automatically.

	// ptclCM = &ptclGroup->sym_int.particles.cm;
	for (int dim=0; dim<Dim; dim++) {
		ptclCM->Position[dim] = ptclGroup->sym_int.particles.cm.Position[dim];
		ptclCM->Velocity[dim] = ptclGroup->sym_int.particles.cm.Velocity[dim];
		ptclCM->Mass = ptclGroup->sym_int.particles.cm.Mass;
	}


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

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	CalculateAcceleration01(ptclCM);
	CalculateAcceleration23(ptclCM);

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

#endif

