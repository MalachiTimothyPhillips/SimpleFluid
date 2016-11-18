//
// Created by Malachi Phillips on 9/30/16.
//

#include "SolutionProcedure.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iostream>
#include <iomanip>

//============================================================================================================
/*
 * Solution Procedure - shared by all other equations
 */
//============================================================================================================

//============================================================================================================
void SolutionProcedure::start_procedure(std::string& runtime_params, std::string& template_file_name, std::string& whichPDE ){
    // Universal parser

    RuntimeParser* runtime = new RuntimeParser;
    runtime->read_parameters_from_file(runtime_params);
    std::vector<double> myArgs;
    myArgs = runtime->return_args();

    std::string condition;
    condition = runtime->return_cond();

    solutionInitialCondition_ = solutionInitialCondition_->make_initial_condition(condition, myArgs);

    solutionInitialCondition_->fluidEquation_ =
            solutionInitialCondition_->fluidEquation_->make_fluid_equation(whichPDE, myArgs);

    solutionInitialCondition_->apply_initial_cond();
    solutionInitialCondition_->enforce_boundary();

    // execute procedure
    procedure(template_file_name);

    // end procedure
    end_procedure();

}

//============================================================================================================
void SolutionProcedure::procedure(std::string& template_file_name){
    // Write out initial step
    solutionInitialCondition_->fluidEquation_->write_to_file(template_file_name,0);
    // Apply step at every time step
    for(unsigned int t = 1; t < solutionInitialCondition_->fluidEquation_->get_nt(); ++t){
        // Enforce boundary prior to step, this *shouldn't* be needed, but just in case
        solutionInitialCondition_->enforce_boundary();
        solutionInitialCondition_->fluidEquation_->apply_step();
        // Write out solutions
        solutionInitialCondition_->fluidEquation_->write_to_file(template_file_name,t);
    }
}

//============================================================================================================
void SolutionProcedure::end_procedure(){
    // something important to finish writing the output file
    std::cout << "Done." << std::endl;
}