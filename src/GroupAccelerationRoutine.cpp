#include <iostream>
#include "global.h"
#include "defs.h"
#include <cassert> // Eunwoo debug

// Eunwoo edited

void NewFBInitialization(Group* group, std::vector<Particle*> &particle);
void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle);
void MergeGroups(std::vector<Group*> &groups);
bool FindCMParticle(Group* group);
bool isNeighborInsideGroup(Particle* member, const std::vector<Particle*>& groupMembers);

bool AddNewGroupsToList(std::vector<Particle*> &particle) {

	assert(GroupCandidateList.empty()); // Let's check GroupCandidateList is initially empty!

	for (Particle *ptcl : particle) {
		ptcl->isFBCandidate();
	}

	if (GroupCandidateList.empty())
		return true;

	MergeGroups(GroupCandidateList);	// Merge GroupCandidateList
										// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!

	for (Group *groupCandidate : GroupCandidateList) {

        // fprintf(binout, "groupCandidate Member: ");
        // for (Particle *p : groupCandidate->Members) {
        //     fprintf(binout, "PID: %d ", p->PID);
        // }
        // fprintf(binout, "\n");

		bool cmInGroup = FindCMParticle(groupCandidate);
        // fprint(binout, "cmInGroup>");

		if (!cmInGroup) { // No CM particle inside a groupCandidate --> Simply just make new group!

			NewFBInitialization(groupCandidate, particle); // New group added.
            continue;

		} else {	// 1. CM + particles case
					// 2. CM + CM + particles case --> Add newly detected members to the existing group!
					// In these cases, replace CM particles to their member particles inside the current groupCandidate vector.
					// Then, use FBTermination to those CM particles and use NewFBInitialization.


            std::vector<Particle*> particlesToErase;  // Collect particles to be erased later

            for (Particle *member : groupCandidate->Members) {

                if (!member->isCMptcl) {
                    continue;
                } else {
                    // Add all group members from member->GroupInfo->Members to groupCandidate
                    groupCandidate->Members.insert(groupCandidate->Members.end(), member->GroupInfo->Members.begin(), member->GroupInfo->Members.end());

                    // Mark this CM particle for erasure by adding it to a separate vector
                    particlesToErase.push_back(member);

                }
            }

            // Now erase particles and terminate them
            for (Particle *member : particlesToErase) {
                member->isErase = true;

                // Remove the marked particle from groupCandidate->Members
                groupCandidate->Members.erase(
                    std::remove_if(groupCandidate->Members.begin(), groupCandidate->Members.end(),
                        [](Particle* p) {
                            return p->isErase;
                        }),
                    groupCandidate->Members.end());

                member->isErase = false;  // Reset the erase flag (if necessary)

                // Terminate the member
                FBTermination(member, particle);
            }
            particlesToErase.clear();
            particlesToErase.shrink_to_fit();

			NewFBInitialization(groupCandidate, particle);
		}		
	}
	GroupCandidateList.clear();
	GroupCandidateList.shrink_to_fit();

	return true;
}

void GroupAccelerationRoutine(REAL next_time, std::vector<Particle*> &particle) {

	if (next_time == 0) {
		return;
	}

	for (Group* ptclGroup: GroupList) {
		ptclGroup->ARIntegration(next_time, particle);
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

bool isNeighborInsideGroup(Particle* member, const std::vector<Particle*>& groupMembers) {
    // Check if every neighbor of 'member' is also inside groupMembers
    for (Particle* neighbor : member->ACList) {
        if (std::find(groupMembers.begin(), groupMembers.end(), neighbor) == groupMembers.end()) {
            // Found a neighbor that is not in groupMembers
            return false;
        }
    }
    return true; // All neighbors are inside groupMembers
}