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
 * Factory function for creating runtime parameter handler
 */
//============================================================================================================
//RuntimeParameters* RuntimeParameters::create_runtime_param(std::string &whichPDE) {
//    if (whichPDE == "LinearWaveEquation"){
//        return new RuntimeParametersLinWave;
//    }
//
//    if (whichPDE == "DiffusionEquation"){
//        return new RuntimeParametersDiffusion;
//    }
//
//    else{
//        // Should never reach here -- this was already checked!
//        return NULL;
//    }
//}

//============================================================================================================
/*
 * Runtime parameter handler for linear wave equation
 */
//============================================================================================================

//============================================================================================================
void RuntimeParametersLinWave::read_parameters_from_file(std::string& runtime_file){


    /*
     * Read in runtime parameters.
     * Current method: user specifies the length of the wave in from of endpoints, xo and xf
     * User specifies begin time (let's say, to, although this is arbitrary)
     * User specifies number of time points
     * User specifies step size for time
     * From CFL (or safetyFactory * CFL), the step size is inferred
     */

    // Read and set parameters from file --  pass in order of
    /*
     * to
     * num t steps
     * dt
     * xo
     * xf
     * c
     * value at LHS wall
     *
     */
    std::ifstream myInputFile;
    myInputFile.open(runtime_file);
    if (!myInputFile.is_open()){
        std::cout << "Error: cannot open file " << runtime_file << std::endl;
        std::cout << "Please ensure that the file name was properly entered.";
        // should throw runtime error here
    }
    std::string currentLine;
    std::vector<double> inputParameters;
    // iterate over lines of file
    while(std::getline(myInputFile,currentLine)){
        // string to double conversion
        double valueFromFile = std::stod(currentLine);

        // write current value to holder
        inputParameters.push_back(valueFromFile);
    }

    // Check that size of input parameters is correct -- 8
    if (inputParameters.size() != 8){
        std::cout << "Error: number of inputs in file " << runtime_file << " is not 8" << std::endl;
        // should throw runtime error here
    }

    // In order of
    /*
     * to
     * num time steps
     * dt
     * xo
     * xf
     * c
     * value at LHS of wall
     * CFL
     *
     */

    set_to(inputParameters[0]);
    set_time_iterations((unsigned int)inputParameters[1]);
    set_dt(inputParameters[2]);
    set_xo(inputParameters[3]);
    set_xf(inputParameters[4]);
    set_c(inputParameters[5]);
    set_wall_boundary(inputParameters[6]);
    set_CFL(inputParameters[7]);

    // calculate important parameters from user specified
    //set_CFL();
    set_space_iterations();


    //TODO: error trapping for inputs, ie. make sure that time does not reverse, etc.
}

//============================================================================================================
void RuntimeParametersLinWave::set_CFL(){
    // For constant c ie. ut + cux = 0
    // Exact solution is obtained for CFL of 1
    // Ergo, force CFL to be 1 for this case, regardless of
    // user entry

    CFL_ = 1.0;

    // warn user of this decision
    std::cout << "CFL of: " << CFL_ << " is being used." << std::endl;
    std::cout << "For constant c, ut + cux = 0, CFL of 1 gives exact solution" << std::endl;

}

//============================================================================================================
void RuntimeParametersLinWave::set_space_iterations(){

    // will calculate needed dx to ensure CFL is unity
    dx_ = c_*dt_/CFL_;
    space_iterations_ = (int)((xf_-xo_)/dx_); // not sure if floor or ceiling - needs to be ceiling
    // Check this

}

//============================================================================================================
/*
 * Runtime parameter handler for diffusion PDE
 */
//============================================================================================================

//============================================================================================================
void RuntimeParametersDiffusion::read_parameters_from_file(std::string& runtime_file){


    /*
     * Read in runtime parameters.
     * For diffusion equation:
     * User picks to, tf, xo, xf, dx, v, value at LHS wall, safety factor
     * dt is picked from the stability analysis and safety factor
     *
     */

    // Read and set parameters from file --  pass in order of
    /*
     * to
     * tf
     * xo
     * xf
     * dx
     * v
     * value at LHS wall
     * value at RHS wall
     * safety factor
     */
    std::ifstream myInputFile;
    myInputFile.open(runtime_file);
    if (!myInputFile.is_open()){
        std::cout << "Error: cannot open file " << runtime_file << std::endl;
        std::cout << "Please ensure that the file name was properly entered.";
        // should throw runtime error here
    }
    std::string currentLine;
    std::vector<double> inputParameters;
    // iterate over lines of file
    while(std::getline(myInputFile,currentLine)){
        // string to double conversion
        double valueFromFile = std::stod(currentLine);

        // write current value to holder
        inputParameters.push_back(valueFromFile);
    }

    // Check that size of input parameters is correct -- 9
    if (inputParameters.size() != 9){
        std::cout << "Error: number of inputs in file " << runtime_file << " is not 9" << std::endl;
        // should throw runtime error here
    }

    // In order of
    /*
     * to
     * number iterations
     * xo
     * xf
     * dx
     * v
     * value at LHS wall
     * value at RHS wall
     * safety factor
     */

    set_to(inputParameters[0]);
    set_time_iterations(inputParameters[1]);
    set_xo(inputParameters[2]);
    set_xf(inputParameters[3]);
    set_dx(inputParameters[4]);
    set_nu(inputParameters[5]);
    add_wall_value(inputParameters[6]);
    add_wall_value(inputParameters[7]);
    set_safety_factor(inputParameters[8]);

    // calculate important parameters from user specified
    set_dt();
    set_space_iterations();

    set_alpha();

    //TODO: error trapping for inputs, ie. make sure that time does not reverse, etc.
}

//============================================================================================================
void RuntimeParametersDiffusion::set_space_iterations(){
    space_iterations_ = (int)((xf_-xo_)/dx_); // not sure if floor or ceiling - needs to be ceiling
    // Check this
}

//============================================================================================================
void RuntimeParametersDiffusion::set_dt(){
    dt_ = 0.5 * safetyFactor_ * dx_ * dx_ / nu_;
}

////============================================================================================================
//void RuntimeParametersDiffusion::set_time_iterations(){
//    space_iterations_ = (int)(ceil(((tf_-to_)/dt_)));
//}

//============================================================================================================
void RuntimeParametersDiffusion::set_alpha(){
    alpha_ = nu_ * dt_ / dx_ / dx_;
}

//============================================================================================================
void RuntimeParametersDiffusion::set_wall_boundary(std::vector<double>& wallVals){
    for(unsigned int i = 0 ; i < wallVals.size(); ++i){
        wallValues_[i] = wallVals[i];
    }
}


//============================================================================================================
/*
 * Runtime parameter handler for Burger PDE
 */
//============================================================================================================

//============================================================================================================
void RuntimeParametersBurger::read_parameters_from_file(std::string& runtime_file){


    /*
     * Read in runtime parameters.
     * For diffusion equation:
     * User picks to, tf, xo, xf, dx, v, value at LHS wall, safety factor
     * dt is picked from the stability analysis and safety factor
     *
     */

    // Read and set parameters from file --  pass in order of
    /*
     * to
     * tf
     * xo
     * xf
     * dx
     * v
     * value at LHS wall
     * value at RHS wall
     * safety factor
     */
    std::ifstream myInputFile;
    myInputFile.open(runtime_file);
    if (!myInputFile.is_open()){
        std::cout << "Error: cannot open file " << runtime_file << std::endl;
        std::cout << "Please ensure that the file name was properly entered.";
        // should throw runtime error here
    }
    std::string currentLine;
    std::vector<double> inputParameters;
    // iterate over lines of file
    while(std::getline(myInputFile,currentLine)){
        // string to double conversion
        double valueFromFile = std::stod(currentLine);

        // write current value to holder
        inputParameters.push_back(valueFromFile);
    }

    // Check that size of input parameters is correct -- 10
    if (inputParameters.size() != 10){
        std::cout << "Error: number of inputs in file " << runtime_file << " is not 10" << std::endl;
        // should throw runtime error here
    }

    // In order of
    /*
     * to
     * number iterations, time
     * xo
     * xf
     * dx
     * v
     * umax
     * value at LHS wall
     * value at RHS wall
     * safety factor
     */

    set_to(inputParameters[0]);
    set_time_iterations(inputParameters[1]);
    set_xo(inputParameters[2]);
    set_xf(inputParameters[3]);
    set_dx(inputParameters[4]);
    set_nu(inputParameters[5]);
    set_cmax(inputParameters[6]);
    add_wall_value(inputParameters[7]);
    add_wall_value(inputParameters[8]);
    set_safety_factor(inputParameters[9]);

    // calculate important parameters from user specified
    set_dt();
    set_space_iterations();

    set_alpha();

    //TODO: error trapping for inputs, ie. make sure that time does not reverse, etc.
}

//============================================================================================================
void RuntimeParametersBurger::set_space_iterations(){
    space_iterations_ = (int)((xf_-xo_)/dx_); // not sure if floor or ceiling - needs to be ceiling
    // Check this
}

//============================================================================================================
void RuntimeParametersBurger::set_dt(){
    dt_ = 0.5 * safetyFactor_ * dx_ * dx_ / nu_;
}

////============================================================================================================
//void RuntimeParametersDiffusion::set_time_iterations(){
//    space_iterations_ = (int)(ceil(((tf_-to_)/dt_)));
//}

//============================================================================================================
void RuntimeParametersBurger::set_alpha(){
    alpha_ = nu_ * dt_ / dx_ / dx_;
}

//============================================================================================================
void RuntimeParametersBurger::set_wall_boundary(std::vector<double>& wallVals){
    for(unsigned int i = 0 ; i < wallVals.size(); ++i){
        wallValues_[i] = wallVals[i];
    }
}