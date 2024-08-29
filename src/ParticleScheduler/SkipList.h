#ifndef SKIP_LIST_H
#define SKIP_LIST_H
#include "../global.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <vector>
#include <omp.h>

//class Node {
struct Node {
	ULL value;
	Node** forward;
	Pvector particle_list;

	Node(ULL value, int level, Particle* ptcl) {
		this->value = value;
		particle_list.push_back(ptcl);
		forward = new Node*[level + 1];
		memset(forward, 0, sizeof(Node*) * (level + 1));
	}

	Node(ULL value, int level) {
		this->value = value;
		forward = new Node*[level + 1];
		memset(forward, 0, sizeof(Node*) * (level + 1));
	}

	~Node() {
		delete[] forward;
		particle_list.clear();
		particle_list.shrink_to_fit();
	}
};

class SkipList {
	private:
		int maxLevel;
		int currentLevel;
		float probability;
		Node* header;

		int randomLevel() {
			int level = 0;
			while (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) < probability && level < maxLevel) {
				level++;
			}
			return level;
		}

	public:


		SkipList(int maxLevel, float probability) {
			this->maxLevel = maxLevel;
			this->probability = probability;
			currentLevel = 0;
			//header = nullptr;
			header = new Node(-1, maxLevel);
		}

		SkipList(int maxLevel, float probability, Particle *ptcl) {
			this->maxLevel = maxLevel;
			this->probability = probability;
			currentLevel = 0;
			header = new Node(-1, maxLevel, ptcl);
		}


		~SkipList() {
			delete header;
		}

		void insert(ULL value, Particle* ptcl) {
			//std::cout << "Inserting value: " << value << "\n";
			Node* current = header;
			Node* update[maxLevel + 1];

			for (int i = currentLevel; i >= 0; i--) {
				while (current->forward[i] != nullptr && current->forward[i]->value < value) {
					current = current->forward[i];
				}
				update[i] = current;
			}

			current = current->forward[0];

			if (current == nullptr || current->value != value) {
				int randomLevelGenerated = randomLevel();

				if (randomLevelGenerated > currentLevel) {
					for (int i = currentLevel + 1; i <= randomLevelGenerated; i++) {
						update[i] = header;
					}
					currentLevel = randomLevelGenerated;
				}

				Node* newNode = new Node(value, randomLevelGenerated, ptcl);

				for (int i = 0; i <= randomLevelGenerated; i++) {
					newNode->forward[i] = update[i]->forward[i];
					update[i]->forward[i] = newNode;
				}
				//std::cout << "Inserted value: " << value << "\n";
			}
		}


		void insert(ULL value) {
			Node* current = header;
			Node* update[maxLevel + 1];

			for (int i = currentLevel; i >= 0; i--) {
				while (current->forward[i] != nullptr && current->forward[i]->value < value) {
					current = current->forward[i];
				}
				update[i] = current;
			}

			current = current->forward[0];

			if (current == nullptr || current->value != value) {
				int randomLevelGenerated = randomLevel();

				if (randomLevelGenerated > currentLevel) {
					for (int i = currentLevel + 1; i <= randomLevelGenerated; i++) {
						update[i] = header;
					}
					currentLevel = randomLevelGenerated;
				}

				Node* newNode = new Node(value, randomLevelGenerated);

				for (int i = 0; i <= randomLevelGenerated; i++) {
					newNode->forward[i] = update[i]->forward[i];
					update[i]->forward[i] = newNode;
				}
				//std::cout << "Inserted value: " << value << "\n";
			}
		}





		// Search for a value in the skip list
		bool search(ULL value, Particle * ptcl) {
			//std::cout << "Parallel search started" <<  "\n";
			Node* current = header;
			bool found;

			// Start from the highest level and move forward
			for (int i = currentLevel; i >= 0; i--) {
				while (current->forward[i] != nullptr && current->forward[i]->value < value) {
					current = current->forward[i];
				}
			}

			current = current->forward[0];  // Move to level 0

			found =  current != nullptr && current->value == value;

			if (found)
				current->particle_list.push_back(ptcl);

			//std::cout << "Parallel search end" <<  std::endl;
			return found;  // Return true if found
		}



		bool parallel_search(ULL value, Particle* ptcl) {
			//std::cout << "Parallel search started" <<  "\n";
			bool found = false;
#pragma omp parallel
			{
#pragma for shared(found)
				for (int i = currentLevel; i >= 0; i--) {
					if (found) continue;

					Node* current = header;
					while (current->forward[i] != nullptr && current->forward[i]->value < value) {
						current = current->forward[i];
					}
					if (current->forward[i] != nullptr && current->forward[i]->value == value) {
#pragma omp critical
						{
							found = true;
							current->forward[i]->particle_list.push_back(ptcl);
						}
					}
				}
			}
			//std::cout << "Parallel search done" <<  "\n";
			return found;
		}



		Node *getFirstNode() {
			return header->forward[0];
		}



		void deleteFirstNode() {
			if (header->forward[0] == nullptr) {
				std::cout << "Skip list is empty\n";
				return;
			}

			Node* oldFirstNode = header->forward[0];

			for (int i = 0; i <= currentLevel; i++) {
				if (oldFirstNode == header->forward[i]) {
					header->forward[i] = header->forward[i]->forward[i];
				}
			}

			while (currentLevel > 0 && header->forward[currentLevel] == nullptr) {
				currentLevel--;
			}

			delete oldFirstNode;

			//std::cout << "Deleted the old header and updated the header to the next node\n";
		}



		size_t getNodeSize() const {
			return sizeof(Node) + sizeof(Node*) * (header->forward ? maxLevel + 1 : 0);
		}

		size_t getTotalSize() const {
			size_t totalSize = 0;
			Node* current = header->forward[0];

			while (current != nullptr) {
				totalSize += getNodeSize();
				current = current->forward[0];
			}

			return totalSize;
		}



		void display() {
			std::cout << "CurrentLevel " << currentLevel << std::endl;
			for (int i = 0; i <= currentLevel; i++) {
				Node* node = header->forward[i];
				std::cout << "Level " << i << ": ";
				while (node != nullptr) {
					std::cout << node->value << " ";
					std::cout << node->particle_list.size() << ":";
					for (Particle*ptcl:node->particle_list) {
						std::cout << ptcl->PID << " ";
					}
					std::cout << ", ";
					node = node->forward[i];
				}
				std::cout << std::endl;
			}
		}
};

#endif
