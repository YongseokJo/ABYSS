#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include "global.h"

int getLineNumber();
void write_out(std::ofstream& outputFile, const Particle* ptcl);
void write_out_group(std::ofstream& outputFile, const Particle* ptclCM, const Particle* ptcl);
void write_neighbor(std::ofstream& outputFile, const Particle* ptcl); 
const int NUM_COLUMNS = 7; // Define the number of columns
const int width = 18;

int readData(std::vector<Particle*> &particle) {

	fprintf(stdout, "Opening %s ...\n", fname);
	std::ifstream inputFile(fname);

	if (!inputFile) {
		std::cerr << "Error: Could not open the file." << std::endl;
		return FAIL;
	}


	int NumParticle;
	NumParticle = getLineNumber();
	NNB = NumParticle;

	// Declaration
	//Particle *particle_temp;	
	//particle_temp = new Particle[NumParticle];
	REAL** data = new REAL*[NumParticle];

	for (int i = 0; i < NumParticle; ++i) {
		data[i] = new REAL[NUM_COLUMNS];
	}

	// Initialization
	for (int i = 0; i < NumParticle; ++i) {
		for (int j = 0; j < NUM_COLUMNS; ++j) {
			data[i][j] = 0;
		}
	}


	int row = 0;

	std::string line;
	while (std::getline(inputFile, line) && row < NumParticle) { // Read lines from the file
		std::istringstream iss(line); // Create a stringstream for each line

		REAL value;
		int col = 0;
		while (iss >> value && col < NUM_COLUMNS) { // Read values from the stringstream
			data[row][col] = value;
			++col;
		}
		//particle_temp[row].setParticleInfo(data[row], row);
		//particle.push_back(new Particle()particle_temp[row]);
		particle.push_back(new Particle(data[row],row));
		++row;
	}

	/*
	for (int col = 0; col < NUM_COLUMNS; ++col) {
		std::cout << "Column " << col + 1 << " values: ";
		for (int r = 0; r < row; ++r) {
			std::cout << data[col][r] << ' ';
		}
		std::cout << std::endl;
	}
	*/


	// Normalize particles
	std::cout << "Particle normalizing." << std::endl;
	for (Particle* element:particle) {
		element->normalizeParticle();
	}
	inputFile.close();

	for (int i=0; i<particle.size(); i++) {
		particle[i]->ParticleOrder = i;
	}



	// Deallocate memory
	for (int i = 0; i < NumParticle; ++i) {
		delete[] data[i];
	}
	delete[] data;


	return DONE;
}


int getLineNumber() {
	    std::ifstream inputFile(fname); // Open the file

			if (!inputFile) {
				std::cerr << "Error: Could not open the file." << std::endl;
				return 1;
			}

			int lineCount = 0;
			std::string line;
			while (std::getline(inputFile, line)) { // Read lines from the file
				lineCount++;
			}

			std::cout << "Number of lines in the file: " << lineCount << std::endl;

			inputFile.close(); // Close the file

			return lineCount;
}


int WriteData() {
	return DONE;
}



// Function to create a directory

bool createDirectory(const std::string& path) {
	// Create a folder with permissions 0777 (full access for user, group, others)
	int status = mkdir(path.c_str(), 0777);

	if (status == 0) {
		std::cout << "Folder created successfully." << std::endl;
	} else {
		std::cerr << "Error creating folder." << std::endl;
		// You can use perror to print the error message for more details
		perror("mkdir");
	}
	return true;
}



int writeParticle(std::vector<Particle*> &particle, REAL current_time, int outputNum) {

    std::cout << "Data is being written..." << std::endl;
    std::string directoryPath = "output";

    // Create the directory or check if it already exists
    if (!createDirectory(directoryPath)) {
        // Handle the error if necessary
        return 1;
    }


    // Now let's save the outputs in a new directory

    // Construct the filename with the timestamp
    std::string filename = directoryPath + "/" + foutput + "_" + std::to_string(outputNum) + ".txt";
    std::string nn_fname = directoryPath + "/neighbor/nn_" + std::to_string(outputNum) + ".txt";

    // Open a file for writing
    std::ofstream outputFile(filename);
    std::ofstream output_nn(nn_fname);


    // Check if the file is opened successfully
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        return 1;
    }

		outputFile << current_time*EnzoTimeStep*1e10/1e6 << " Myr, "; //
		//outputFile << global_time*EnzoTimeStep*1e10/1e6 << " Myr"; //
		outputFile << "\n";
		outputFile << outputTime << ", "; //
		outputFile << outputTimeStep << ", "; //
		outputFile << global_time << ""; //
		outputFile << "\n";
    outputFile << std::left 
			<< std::setw(width) << "PID"
			<< std::setw(width) << "Mass (Msun)"
			<< std::setw(width) << "X (pc)"
			<< std::setw(width) << "Y (pc)"
			<< std::setw(width) << "Z (pc)"
			<< std::setw(width) << "Vx (km/s)"
		 	<< std::setw(width) << "Vy (km/s)" 
			<< std::setw(width) << "Vz (km/s)" << "\n";


    // Write particle data to the file 
		for (Particle* ptcl:particle) {
			ptcl->predictParticleSecondOrderIrr(current_time);
			if (ptcl->isCMptcl)  { // Eunwoo edited
				for (Particle* members : ptcl->GroupInfo->Members) { // Eunwoo edited
					write_out_group(outputFile, ptcl, members); // Eunwoo edited
				} // Eunwoo edited
				// //write_neighbor(output_nn, ptcl->BinaryParticleI); // Eunwoo edited
				// write_out(outputFile, ptcl->BinaryParticleJ); // Eunwoo edited
				// //write_neighbor(output_nn, ptcl->BinaryParticleJ); // Eunwoo edited
			}
			else {
				write_out(outputFile, ptcl);
				//write_neighbor(output_nn, ptcl);
			}
			//write_out(outputFile, ptcl);
			//write_neighbor(output_nn, ptcl);
    }

    // Close the file
    outputFile.close();
    output_nn.close();

    std::cout << "Data written to output.txt successfully!" << std::endl;

    return 0;

}


void write_out(std::ofstream& outputFile, const Particle* ptcl) {
        outputFile  << std::left
										<< std::setw(width) << ptcl->PID
										<< std::setw(width) << ptcl->Mass*mass_unit
                    << std::setw(width) << ptcl->PredPosition[0]*position_unit
                    << std::setw(width) << ptcl->PredPosition[1]*position_unit
                    << std::setw(width) << ptcl->PredPosition[2]*position_unit
                    << std::setw(width) << ptcl->PredVelocity[0]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << ptcl->PredVelocity[1]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << ptcl->PredVelocity[2]*velocity_unit/yr*pc/1e5 << '\n';
}

// This function is for group members cause group members have pos, vel in original frame, not predicted values.
void write_out_group(std::ofstream& outputFile, const Particle* ptclCM, const Particle* ptcl) {
        outputFile  << std::left
										<< std::setw(width) << ptcl->PID
										<< std::setw(width) << ptcl->Mass*mass_unit
                    << std::setw(width) << (ptclCM->PredPosition[0] - ptclCM->Position[0] + ptcl->Position[0])*position_unit
                    << std::setw(width) << (ptclCM->PredPosition[1] - ptclCM->Position[1] + ptcl->Position[1])*position_unit
                    << std::setw(width) << (ptclCM->PredPosition[2] - ptclCM->Position[2] + ptcl->Position[2])*position_unit
                    << std::setw(width) << (ptclCM->PredVelocity[0] - ptclCM->Velocity[0] + ptcl->Velocity[0])*velocity_unit/yr*pc/1e5
                    << std::setw(width) << (ptclCM->PredVelocity[1] - ptclCM->Velocity[1] + ptcl->Velocity[1])*velocity_unit/yr*pc/1e5
                    << std::setw(width) << (ptclCM->PredVelocity[2] - ptclCM->Velocity[2] + ptcl->Velocity[2])*velocity_unit/yr*pc/1e5 << '\n';
}

void write_neighbor(std::ofstream& outputFile, const Particle* ptcl) {
	outputFile  << std::left\
			<< std::setw(width) << ptcl->PID << " = [ " ;
	for (Particle* nn:ptcl->ACList) {
			outputFile << nn->PID << "  ";
	}
	outputFile << "]\n";

}


#ifdef time_trace
void output_time_trace() {

}
#endif
