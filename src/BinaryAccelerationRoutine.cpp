#include <iostream>
#include "global.h"
#include "defs.h"

// Eunwoo edited

void NewFBInitialization(Group* group, std::vector<Particle*> &particle);
void FBModification(std::vector<Particle*> addMembers, std::vector<Particle*> &particle);
void MergeGroups(std::vector<Group*> &groups);
bool FindCMParticle(Group* group);

bool AddNewBinariesToList(std::vector<Particle*> &particle) {

	assert(GroupCandidateList.empty()); // Let's check GroupCandidateList is initially empty!

	for (Particle *ptcl : particle) {
		ptcl->isFBCandidate();
	}

	if (GroupCandidateList.empty())
		return true;

	MergeGroups(GroupCandidateList);	// Merge GroupCandidateList
										// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!

	for (Group *groupCandidate : GroupCandidateList) {

		bool cmInGroup = FindCMParticle(groupCandidate);

		if (cmInGroup == false) { // No CM particle inside a groupCandidate --> Simply just make new group!

			NewFBInitialization(groupCandidate, particle); // New group added.

		} else {	// 1. CM + particles case
					// 2. CM + CM + particles case --> Add newly detected members to the existing group!
					// In these cases, remove CM particles from the particle vector and ACList of particles.
					// Then remove groups related with CM particles from the GroupList
					// FBModification function --> Almost similar to the NewFBInitialization function!

			std::vector<Particle*> membersToAdd;	// members to add in the group
													// CM particles should be splitted up to its individual members, cmIndex too!!
													// group members of cmIndex should be included also because it is impossible to add particles to the existing group in SDAR

			for (Particle *member : groupCandidate->Members) {

				if (!member->isCMptcl) {
					membersToAdd.push_back(member); // If it is not the CM particle, just add to the membersToAdd!
													// They are not included in the already existing groups so not the problem at all
				} else {	// member is CM ptcl and we have to delete it from particle vector and ACList of particle in vector before we move on!!!
							// Because we are going to delete group which is connected to CM ptcl and CM ptcl will be deleted automatically.
					member->isErase = true;
					particle.erase(
							std::remove_if(particle.begin(), particle.end(),
								[](Particle* p) {
								bool to_remove = p->isErase;
								return to_remove;
								}),
							particle.end());

					for (int i = 0; i < particle.size(); i++) {
						particle[i]->ParticleOrder = i;

						auto& acList = particle[i]->ACList;

						acList.erase(
							std::remove_if(acList.begin(), acList.end(),
								[](Particle* p) {
									bool to_remove = p->isErase;
									return to_remove;       // Keep this member in ACList
									}),
							acList.end());
						particle[i]->NumberOfAC = particle[i]->ACList.size();
					}
					member->isErase = false;

					membersToAdd.insert(membersToAdd.end(), member->GroupInfo->Members.begin(), member->GroupInfo->Members.end());

					GroupList.erase(
							std::remove_if(GroupList.begin(), GroupList.end(), 
								[member](Group* g) {
								if (g->groupCM->PID == member->PID) {
									delete g;
									return true;
								}
								return false;
							}),
							GroupList.end()
					);
				}
			}

			// CM particle and its following group which is added to groupCandidate->Members[cmIndex] should be deleted!!!
			// I think I did it in the very above but I'm not sure it is correct or not...

			FBModification(membersToAdd, particle);

			membersToAdd.clear();
			membersToAdd.shrink_to_fit();
		}		
	}

	GroupCandidateList.clear();
	GroupCandidateList.shrink_to_fit();

	return true;
}

void BinaryAccelerationRoutine(REAL next_time) {

	if (next_time == 0) {
		return;
	}

	for (Group* ptclGroup: GroupList) {
		ptclGroup->ARIntegration(next_time);
	}
	// fprintf(stdout, "BinaryAccelerationRoutine ended, current time: %e\n", next_time); // Eunwoo added for debug
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

bool FindCMParticle(Group* group) {
    for (size_t i = 0; i < group->Members.size(); ++i) {
        if (group->Members[i]->isCMptcl) 
            return true; // Return the index of the first CM particle        
    }
    return false; // No CM particle found
}
