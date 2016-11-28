#include <iostream>
#include "IntCond.h"
#include "RuntimeParameters.h"
#include "SolutionProcedure.h"
#include <vector>
#include "RuntimeJoke.h"
#include <ostream>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>

int main(int argc, char *argv[]) {

    // Test user specified input file
    if (argc != 4){
        std::cout << "Error: specify" << std::endl
                  << "Type of PDE" << std::endl
                  << "Inputfile" << std::endl
                  << "Outputfile" << std::endl;
    }

    // Use first string passed to determine what type of equation is being solved
    std::string whichPDE = argv[1];

    // Set PDE Type in SolutionProcedure
    SolutionProcedure* currentProcedure = new SolutionProcedure;
    // create output file
    std::string userOutputFile = argv[3];
    // pass from command line the filename to be used for procedure
    std::string userRuntimeArguments = argv[2];
    currentProcedure->start_procedure(userRuntimeArguments, userOutputFile, whichPDE);

    return 0;
}