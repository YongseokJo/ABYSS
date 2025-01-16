#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include "global.h"

int getLineNumber();
//void write_out(std::ofstream& outputFile, const Particle* ptcl);
void write_out(std::ofstream& outputFile, const Particle* ptcl, const double *pos, const double *vel);
void write_out_group(std::ofstream& outputFile, const Particle* ptcl, const Particle* members, const double *pos, const double *vel);
void write_neighbor(std::ofstream& outputFile, const Particle* ptcl);
const int NUM_COLUMNS = 7; // Define the number of columns
const int width = 18;

int readData() {

	fprintf(stdout, "Opening %s ...\n", fname);
	std::ifstream inputFile(fname);

	if (!inputFile) {
		std::cerr << "Error: Could not open the file." << std::endl;
		return FAIL;
	}

	NumberOfParticle = getLineNumber();
	NewPID = NumberOfParticle;
	LastParticleIndex = NumberOfParticle - 1;

	// Declaration
	//Particle *particle_temp;	
	//particle_temp = new Particle[NumParticle];
	double** data = new double*[NumberOfParticle];

	for (int i = 0; i < NumberOfParticle; ++i) {
		data[i] = new double[NUM_COLUMNS];
	}

	// Initialization
	for (int i = 0; i < NumberOfParticle; ++i) {
		for (int j = 0; j < NUM_COLUMNS; ++j) {
			data[i][j] = 0;
		}
	}


	int row = 0;

	std::string line;
	while (std::getline(inputFile, line) && row < NumberOfParticle) { // Read lines from the file
		std::istringstream iss(line); // Create a stringstream for each line

		double value;
		int col = 0;
		while (iss >> value && col < NUM_COLUMNS) { // Read values from the stringstream
			data[row][col] = value;
			++col;
		}
		//particle_temp[row].setParticleInfo(data[row], row);
		//particle.push_back(new Particle()particle_temp[row]);
		particles_original[row].initialize(data[row],row);
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
	for (int i=0; i<NumberOfParticle; i++) {
		particles[i].normalizeParticle();
	}
	inputFile.close();

	/*
	for (int i=0; i<particle.size(); i++) {
		particle[i]->ParticleIndex = i;
	}
	*/



	// Deallocate memory
	for (int i = 0; i < NumberOfParticle; ++i) {
		delete[] data[i];
	}
	delete[] data;


	return SUCCESS;
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
	return SUCCESS;
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



int writeParticle(double current_time, int outputNum) {

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
    //std::string nn_fname = directoryPath + "/neighbor/nn_" + std::to_string(outputNum) + ".txt";

    // Open a file for writing
    std::ofstream outputFile(filename);
    //std::ofstream output_nn(nn_fname);


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
		outputFile << current_time << ""; //
		outputFile << "\n";
    outputFile << std::left 
			<< std::setw(width) << "PID"
			<< std::setw(width) << "Mass (Msun)"
			<< std::setw(width) << "X (pc)"
			<< std::setw(width) << "Y (pc)"
			<< std::setw(width) << "Z (pc)"
			<< std::setw(width) << "Vx (km/s)"
		 	<< std::setw(width) << "Vy (km/s)" 
#ifdef SEVN
			<< std::setw(width) << "Vz (km/s)"
			<< std::setw(width) << "Type" << "\n";
#else
			<< std::setw(width) << "Vz (km/s)" << "\n";
#endif 


    // Write particle data to the file
		Particle *ptcl;
		double pos[Dim], vel[Dim];
		for (int i=0; i<=LastParticleIndex; i++) {
			ptcl = &particles[i];
			if (ptcl->isCMptcl) continue;
#ifdef FEWBODY
			if (ptcl->isActive)
			{
				ptcl->predictParticleSecondOrder(current_time - ptcl->CurrentTimeIrr, pos, vel);
				write_out(outputFile, ptcl, pos, vel);
				// write_neighbor(output_nn, ptcl);
			}
			else
			{
				Particle* ptclCM = &particles[ptcl->CMPtclIndex];
				ptclCM->predictParticleSecondOrder(current_time - ptclCM->CurrentTimeIrr, pos, vel);
				write_out_group(outputFile, ptclCM, ptcl, pos, vel);
			}
#else
			ptcl->predictParticleSecondOrder(current_time - ptcl->CurrentTimeIrr, pos, vel);
			write_out(outputFile, ptcl, pos, vel);
// write_neighbor(output_nn, ptcl);
#endif
		}

		// Close the file
		outputFile.close();
		// output_nn.close();

		std::cout << "Data written to output.txt successfully!" << std::endl;

		return 0;

}


void write_out(std::ofstream& outputFile, const Particle* ptcl, const double *pos, const double *vel) {
        outputFile  << std::left << std::fixed << std::setprecision(8) // Eunwoo test
					<< std::setw(width) << ptcl->PID
					<< std::setw(width) << ptcl->Mass*mass_unit
                    << std::setw(width) << pos[0]*position_unit
                    << std::setw(width) << pos[1]*position_unit
                    << std::setw(width) << pos[2]*position_unit
                    << std::setw(width) << vel[0]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << vel[1]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << vel[2]*velocity_unit/yr*pc/1e5 << '\n';
}

// This function is for group members cause group members have pos, vel in original frame, not predicted values.
void write_out_group(std::ofstream& outputFile, const Particle* ptclCM, const Particle* ptcl, const double *pos, const double *vel) {
        // outputFile  << std::left
		outputFile  << std::left << std::fixed << std::setprecision(8) // Eunwoo test
					<< std::setw(width) << ptcl->PID
					<< std::setw(width) << ptcl->Mass*mass_unit
                    << std::setw(width) << (pos[0] - ptclCM->Position[0] + ptcl->Position[0])*position_unit
                    << std::setw(width) << (pos[1] - ptclCM->Position[1] + ptcl->Position[1])*position_unit
                    << std::setw(width) << (pos[2] - ptclCM->Position[2] + ptcl->Position[2])*position_unit
                    << std::setw(width) << (vel[0] - ptclCM->Velocity[0] + ptcl->Velocity[0])*velocity_unit/yr*pc/1e5
                    << std::setw(width) << (vel[1] - ptclCM->Velocity[1] + ptcl->Velocity[1])*velocity_unit/yr*pc/1e5;
#ifdef SEVN
		outputFile << std::setw(width) << (vel[2] - ptclCM->Velocity[2] + ptcl->Velocity[2])*velocity_unit/yr*pc/1e5;
		if (ptcl->star == nullptr)
			outputFile << std::setw(width) << "1" << '\n';
		else if (!ptcl->star->amiremnant())
			outputFile << std::setw(width) << int(ptcl->star->getp(Phase::ID)) << '\n';
		else
			outputFile << std::setw(width) << 8+int(ptcl->star->getp(RemnantType::ID)) << '\n';
#else
		outputFile << std::setw(width) << (vel[2] - ptclCM->Velocity[2] + ptcl->Velocity[2])*velocity_unit/yr*pc/1e5 << '\n';
#endif
}

/*
void write_neighbor(std::ofstream& outputFile, const Particle* ptcl) {
	outputFile  << std::left\
			<< std::setw(width) << ptcl->PID << " = [ " ;
	for (Particle* nn:ptcl->ACList) {
			outputFile << nn->PID << "  ";
	}
	outputFile << "]\n";

}
*/


#ifdef time_trace
void output_time_trace() {

}
#endif
