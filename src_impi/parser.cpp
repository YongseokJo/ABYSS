#include <iostream>
#include "global.h"


int Parser(int argc, char *argv[]) {
	for (int i = 1; i < argc; ++i) {

		std::string arg = argv[i];

		if (arg == "-h" || arg == "--help") {
			// Display help message
			std::cout << "Help message: Your program description here" << std::endl;

		} else if (arg == "-f" || arg == "--file") {
			// Check if the next argument exists
			if (i + 1 < argc) {
				fname = argv[i + 1];
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "Error: Missing argument for -f/--file option" << std::endl;
				return -1; // Error: Missing argument
			}

		} else if (arg == "-dtdump" || arg == "--timestep") {
			if (i + 1 < argc) {
				outputTimeStep = std::atof(argv[i + 1]);
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "Error: Missing argument for -t/--double option" << std::endl;
				return -2; // Error: Missing argument
			}

		} else if (arg == "-tend" || arg == "--end") {
			if (i + 1 < argc) {
				endTime = std::atof(argv[i + 1]);
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "Error: Missing argument for -tend/--double option" << std::endl;
				return -2; // Error: Missing argument
			}
		} else if (arg == "-d" || arg == "--dir") {
			if (i + 1 < argc) {
				foutput = argv[i + 1];
				i++; // Skip the next argument as it's already processed
			}
		} else {
			// Handle unrecognized arguments
			std::cerr << "Error: Unrecognized argument: " << arg << std::endl;
			return -11; // Error: Unrecognized argument
		}
	}

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
	else {
		std::cout << "End Time: " << endTime << std::endl;
		std::cout << "EnzoTimeStep = "   << EnzoTimeStep   << std::endl;
		std::cout << "outputTimeStep = " << outputTimeStep << "Myr" <<std::endl;
	}
	//broadcastFromRoot(endTime,)


	return _SUCCESS;
}
