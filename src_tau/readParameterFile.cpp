#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cstring>
#include "global.h"

class Config {
	public:
		void read(const std::string& filename) {
			std::ifstream file(filename);
			if (!file.is_open()) {
				throw std::runtime_error("Unable to open configuration file");
			}

			std::string line;
			while (std::getline(file, line)) {
				if (line.empty() || line[0] == '#') {
					continue;
				}
				trim(line);
				size_t pos = line.find('=');
				if (pos != std::string::npos) {
					std::string key = line.substr(0, pos);
					std::string value = line.substr(pos + 1);
					trim(key);
					trim(value);
					config_data_[key] = value;
				}
			}

			file.close();
		}

		std::string getString(const std::string& key) const {
			if (config_data_.count(key)) {
				return config_data_.at(key);
			} else {
				throw std::runtime_error("Key not found: " + key);
			}
		}

		int getInt(const std::string& key) const {
			return std::stoi(getString(key));
		}

		double getDouble(const std::string& key) const {
			return std::stod(getString(key));
		}


    char* getChar(const std::string& key) const {
			std::string value = getString(key);
			char* result = new char[value.size() + 1];
			std::strcpy(result, value.c_str());
			return result;
		}


	private:
		void trim(std::string& str) const {
			size_t first = str.find_first_not_of(" \t\n\r");
			size_t last = str.find_last_not_of(" \t\n\r");

			if (first == std::string::npos || last == std::string::npos) {
				str.clear();
			} else {
				str = str.substr(first, (last - first + 1));
			}
		}

		std::unordered_map<std::string, std::string> config_data_;
};

int readParameterFile() {
	try {
		Config config;
		config.read("config.txt");

		fname				   = config.getChar("Filename");
		//int num_particles = config.getInt("num_particles");
		endTime        = config.getDouble("StopTime");
		outputTimeStep = config.getDouble("dtOutput");
		foutput			   = config.getChar("OutputDirectory");


		EnzoTimeStep   = endTime/1e10; // endTime should be Myr
		outputTimeStep = outputTimeStep/endTime; // endTime should be Myr

		if (MyRank == ROOT) {
			std::cout << "Starting ABYSS ..." <<  std::endl;
			std::cout << "File name: " << fname << std::endl;
			std::cout << "Output file name: " << foutput << std::endl;
			std::cout << "End Time: " << endTime << std::endl;
			std::cout << "EnzoTimeStep = "   << EnzoTimeStep   << std::endl;
			std::cout << "outputTimeStep = " << outputTimeStep << "Myr" <<std::endl;
		}


	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}

	return SUCCESS;
}

