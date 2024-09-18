#include <iostream>
#include "global.h"
#include "defs.h"
#include <cassert> // Eunwoo debug

// Eunwoo edited

void NewFBInitialization(Group* group, std::vector<Particle*> &particle);
void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle);
void MergeGroups(std::vector<Group*> &groups);
bool FindCMParticle(Group* group);
bool isNeighborInsideGroup(Group* groupCandidate);

bool AddNewGroupsToList(std::vector<Particle*> &particle) {

	assert(GroupCandidateList.empty()); // Let's check GroupCandidateList is initially empty!

	for (Particle *ptcl : particle) {
        // if (ptcl->isCMptcl) continue; // Eunwoo: test
		ptcl->isFBCandidate();
	}

	if (GroupCandidateList.empty()) return true;

	MergeGroups(GroupCandidateList);	// Merge GroupCandidateList
										// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!

	for (Group *groupCandidate : GroupCandidateList) {

        if (isNeighborInsideGroup(groupCandidate)) continue;    // Don't make a real group if it has no neighbor.
                                                                // Its TimeStepIrr will be so big and it can raise errors.

        NewFBInitialization(groupCandidate, particle);
				
	}
	GroupCandidateList.clear();
	GroupCandidateList.shrink_to_fit();

	return true;
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