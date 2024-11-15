#include <iostream>
#include "global.h"
#include "defs.h"
#include <cassert> // Eunwoo debug

// Eunwoo edited

void NewFBInitialization(Group* group, std::vector<Particle*> &particle);
void NewFBInitialization2(Group* group, std::vector<Particle*> &particle);
void MergeGroups(std::vector<Group*> &groups);
bool isNeighborInsideGroup(Group* groupCandidate);
void NewPrimordialBinaries(Group*groupCandidate, std::vector<Particle*> &particle);
void InitializeFBParticle(Particle* FBParticle, std::vector<Particle*> &particle);

bool AddNewGroupsToList(std::vector<Particle*> &particle) {

	assert(GroupCandidateList.empty()); // Let's check GroupCandidateList is initially empty!

	for (Particle *ptcl : particle) {
        // if (ptcl->isCMptcl) continue; // Eunwoo: test
		// ptcl->isFBCandidate();

        // Eunwoo test
        if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 > 1e-5)
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

void FindPrimordialBinaries(std::vector<Particle*> &particle) {

    assert(GroupCandidateList.empty()); // Let's check GroupCandidateList is initially empty!

	for (Particle *ptcl : particle) 
        ptcl->checkNewGroup2();

    if (GroupCandidateList.empty()) 
        return;

	MergeGroups(GroupCandidateList);	// Merge GroupCandidateList
										// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!

	for (Group *groupCandidate : GroupCandidateList) 
        NewPrimordialBinaries(groupCandidate, particle);

	GroupCandidateList.clear();
	GroupCandidateList.shrink_to_fit();

	return;
}


void MergeGroups(std::vector<Group*> &groups) {

    bool merged = true;

    while (merged) {
        merged = false;
        for (size_t i = 0; i < groups.size(); ++i) {
            Group* currentGroup = groups[i];
			if (!currentGroup) continue; // Skip already deleted groups

            std::vector<Particle*> &group1 = currentGroup->Members;

            for (size_t j = i + 1; j < groups.size(); ++j) {
                Group* otherGroup = groups[j];
				if (!otherGroup) continue; // Skip already deleted groups

                std::vector<Particle*> &group2 = otherGroup->Members;

                // Check if there's any common member between group1 and group2
                bool commonFound = false;
                for (Particle* member1 : group1) {
                    if (std::find(group2.begin(), group2.end(), member1) != group2.end()) {
                        commonFound = true;
                        break;
                    }
                }

                // If common members are found, merge group2 into group1
                if (commonFound) {
                    // Merge group2 into group1, avoiding duplicates
                    for (Particle* member : group2) {
                        if (std::find(group1.begin(), group1.end(), member) == group1.end()) {
                            group1.push_back(member);
                        }
                    }

                    // Mark otherGroup for deletion after the loop
                    delete otherGroup;
                    groups[j] = nullptr; // Mark the position as null
                    merged = true;
                }
            }
        }

        // Clean up empty groups and delete them
        groups.erase(
            std::remove_if(groups.begin(), groups.end(), [](Group* g) {
                if (g == nullptr) return true; // If marked for deletion
                if (g->Members.empty()) {
                    delete g; // Deallocate memory for empty groups
                    return true;
                }
                return false;
            }),
            groups.end()
        );
    }
}

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

// Make CM particles for primordial binaries
void NewPrimordialBinaries(Group* group, std::vector<Particle*> &particle) {

	Particle *ptclCM;
	Group *ptclGroup;

	std::cout <<"\n\n\nStarting Routine NewPrimordialBinaries" << std::endl;

	// Set ptclGroup members first; this will be very useful

	ptclGroup = new Group();
	for (Particle* members : group->Members) 
        ptclGroup->Members.push_back(members);

	// Let's link CM particle with the cm particles made in the binary tree (SDAR).

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(); // Binary tree is made and CM particle is made automatically.

	ptclCM = &ptclGroup->sym_int.particles.cm;

	// Set ptcl information like time, PID, etc.

	ptclCM->PID             = -(std::abs(group->Members[0]->PID) + NNB);
	ptclCM->GroupInfo		= ptclGroup;
	ptclCM->isCMptcl        = true;

	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);


	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (Particle* members: ptclGroup->Members) {
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

	ptclGroup->groupCM		= ptclCM;

	ptclGroup->sym_int.initialIntegration(0); // This is primordial binary!
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	// Erase group particles from particle vector because now they should be replaced with CM particle

	for (Particle* members : group->Members) {
		members->isErase = true;
	}

	particle.erase(
			std::remove_if(particle.begin(), particle.end(),
				[](Particle* p) {
				bool to_remove = p->isErase;
				return to_remove;
				}),
			particle.end());

	// Update particle order and put CM particle into the last of particle vector
	for (int i=0; i<particle.size(); i++) {
		particle[i]->ParticleOrder = i;
	}
	ptclCM->ParticleOrder = particle.size();
	particle.push_back(ptclCM);

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	InitializeFBParticle(ptclCM, particle);


// /* // Eunwoo test
	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, ptclCM->RadiusOfAC);
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


	int size = 0;
	for (Particle* ptcl: particle) {
		size = ptcl->ACList.size();

		ptcl->ACList.erase(
				std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					return to_remove; }),
				ptcl->ACList.end());
		ptcl->NumberOfAC = ptcl->ACList.size();

		if (size != ptcl->NumberOfAC) {
			ptcl->ACList.push_back(ptclCM);
			ptcl->NumberOfAC++;
			// InitializeFBParticle(ptcl, particle); // Eunwoo added
			// ptcl->calculateTimeStepReg(); // Eunwoo added
			// ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr); // Eunwoo added
		}
	}

	for (Particle* members : group->Members) {
		members->isErase = false;
		if (members->isCMptcl) {
			// members->GroupInfo->isErase = true;
			// GroupList.erase(
			// 		std::remove_if(GroupList.begin(), GroupList.end(),
			// 			[](Group* p) {
			// 			bool to_remove = p->isErase;
			// 			//if (to_remove) delete p;
			// 			return to_remove;
			// 			}),
			// 		GroupList.end());
			delete members->GroupInfo;
			members->GroupInfo = nullptr;
			members  = nullptr;
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

	delete group;

	fprintf(binout, "------------------END-OF-NEW-PRIMORDIAL-BINARIES------------------\n\n");
	fflush(binout);
}