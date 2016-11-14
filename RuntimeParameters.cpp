//
// Created by Malachi Phillips on 9/25/16.
//

// C++: where friends can access your private members!

#include "RuntimeParameters.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>

//============================================================================================================
/*
 * Runtime parser
 */
//============================================================================================================
void RuntimeParser::read_parameters_from_file(std::string &runtime_file) {
    /*
     * Runtime paramter parser::
     * User specifies in this order:
     * lo - LHS boundary on x
     * lf - RHS boundary on x
     * nl - Number of X steps
     * ho - Bottom boundary on y
     * hf - Top boundary on y
     * nh - Number of Y steps
     * tf - Simulation time
     * nt - Number of simulation time steps to be advanced
     * c - Constant wave speed, or diffusion coefficient
     * eps - Numerical stability for Laplace, maximum CFL for adaptive time codes
     * wfreq - Frequency with which to output to file
     * InitialCondition - name of initial condition, to be found from enum list
     */

    std::ifstream myInputFile;
    myInputFile.open(runtime_file);
    if (!myInputFile.is_open()){
        std::cout << "Error: cannot open file " << runtime_file << std::endl;
        std::cout << "Please ensure that the file name was properly entered.";
        // should throw runtime error here
    }
    std::string currentLine;
    std::vector<double> args;
    unsigned int i = 0;
    // iterate over lines of file
    while(std::getline(myInputFile,currentLine)){
        if (i < 11) {
            // string to double conversion
            double valueFromFile = std::stod(currentLine);
            // write current value to holder
            args.push_back(valueFromFile); //convert type afterward
        }
        if (i >= 11){
            initialCondition_ = currentLine;
        }
        i+=1;
    }


    //TODO: error trapping for inputs, ie. make sure that time does not reverse, etc.
}
