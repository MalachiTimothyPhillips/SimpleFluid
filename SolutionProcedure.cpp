//
// Created by Malachi Phillips on 9/30/16.
//

#include "SolutionProcedure.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

//============================================================================================================
/*
 * Factory function for solution procedures -- handles which equation to solve!
 */
//============================================================================================================

SolutionProcedure* SolutionProcedure::determine_solution_procedure(std::string &equationType){
    // From string, return correct solver equation
    if (equationType == "LinearWaveEquation"){
        return new SolutionProcedureLinWave;
    }
    if (equationType == "FixedLinearWaveEquation"){
        return new SolutionProcedureFixLinWave;
    }
    if (equationType == "InviscidBurger"){
        return new SolutionProcedureInviscid;
    }
    if (equationType == "Burger"){ // I'm quite hungry now
        return new SolutionProcedureBurger;
    }
    if (equationType == "DiffusionEquation") {
        return new SolutionProcedureDiffusion;
    }
    if (equationType == "LinearAdvectionEquation"){
        return new SolutionProcedureMultiDimAdvection;
    }
    else{
        std::cout << "Error: could not correctly find your PDE equation to solve!" << std::endl;
        // TODO: throw some error here that prevents user from continuing
    }

}

//============================================================================================================
/*
 * Solution procedure for linear wave equation
 */
//============================================================================================================

//============================================================================================================
void SolutionProcedureLinWave::start_procedure(std::string& runtime_params, std::string& template_file_name ){

    // Since type is already known, just hard code in the file reader
    runtime_param_ = new RuntimeParametersLinWave;

    runtime_param_->read_parameters_from_file(runtime_params);

    // Will take from input user to set boundary/initial condition type
    // Fow now, hard code it

    set_boundary();
    set_init_cond();

    // Set size of uSolutions_ vector
    uSolutions_.resize(runtime_param_->get_space_iterations());

    // Apply boundary condition, LHSWall
    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // For first step, need to apply the initial condition
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();
    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);

    // execute procedure
    procedure(template_file_name);

    // end procedure
    end_procedure();

}

//============================================================================================================
void SolutionProcedureLinWave::convert_idx_to_pos(unsigned int idx, double& pos){

    // Convert from index to position
    // idx * dx + xo = current position
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();

    pos = (double)(idx) * dx + xo;

}

//============================================================================================================
void SolutionProcedureLinWave::write_to_file()
{

}

//============================================================================================================
void SolutionProcedureLinWave::procedure(std::string& template_file_name){

    // write first point -- initial condition onto file
    double to = runtime_param_->get_to();
    std::ofstream outFile (template_file_name + "0");
    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
        double pos;
        convert_idx_to_pos(i,pos);
        outFile << pos << " ";
        outFile << uSolutions_[i] << std::endl;
    } // writes file in two column format: x and u(x)
    outFile.close();


    // Apply step at every time step
    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
        apply_step();

        // write results to output file, need t's, x's, and u's
        outFile.open(template_file_name + std::to_string(t));
        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
            double pos;
            convert_idx_to_pos(i,pos);
            outFile << pos << " ";
            outFile << uSolutions_[i] << std::endl;
            // close file when done
        }
        outFile.close();

    }

}

//============================================================================================================
void SolutionProcedureLinWave::end_procedure(){
    // something important to finish writing the output file
    std::cout << "Done." << std::endl;
}

//============================================================================================================
void SolutionProcedureLinWave::apply_step(){

    // Apply individual step to solutions

    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // loops over entries, save the one involved in the boundary condition
    // boundary hard coded to LHS of wall

    // need to save solution of u(i-1) at time n
    double uSolnm1 = uSolutions_[0];

    for (unsigned int i = 1; i < runtime_param_->get_space_iterations(); i++){
        // foward in time, backward in space
        uSolutions_[i] = uSolutions_[i] - runtime_param_->get_CFL()*(uSolutions_[i] - uSolnm1);

        // update u(i-1) at time n as current value
        uSolnm1 = uSolutions_[i];
    }

    // Need to write out solutions, but for now, leave that out
}

//============================================================================================================
void SolutionProcedureLinWave::set_boundary(){

    // hard code as LHSWall
    //std::string& myString;
    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
    solutionBoundaryCondition_ = new LHSWall;


    //TODO: should probably allow the user to go ahead and set what boundary condition to be used

}

//============================================================================================================
void SolutionProcedureLinWave::set_init_cond(){

    // hard code as wave type
    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
    //solutionInitialCondition_ = new StepWave;
    solutionInitialCondition_ = new SinWave;

    // TODO: should probably allow the user to go ahead and set the initial condition from the control file

}

//============================================================================================================
/*
 * Solution procedure for linear wave equation, where c may be positive or negative
 */
//============================================================================================================

//============================================================================================================
void SolutionProcedureFixLinWave::start_procedure(std::string& runtime_params, std::string& template_file_name ){

    // Since type is already known, just hard code in the file reader
    runtime_param_ = new RuntimeParametersLinWave;

    runtime_param_->read_parameters_from_file(runtime_params);

    // Will take from input user to set boundary/initial condition type
    // Fow now, hard code it

    set_boundary();
    set_init_cond();

    // Set size of uSolutions_ vector
    uSolutions_.resize(runtime_param_->get_space_iterations());

    // Apply boundary condition, LHSWall
    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // For first step, need to apply the initial condition
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();
    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);

    // execute procedure
    procedure(template_file_name);

    // end procedure
    end_procedure();

}

//============================================================================================================
void SolutionProcedureFixLinWave::convert_idx_to_pos(unsigned int idx, double& pos){

    // Convert from index to position
    // idx * dx + xo = current position
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();

    pos = (double)(idx) * dx + xo;

}

//============================================================================================================
void SolutionProcedureFixLinWave::write_to_file()
{

}

//============================================================================================================
void SolutionProcedureFixLinWave::procedure(std::string& template_file_name){

    // write first point -- initial condition onto file
    double to = runtime_param_->get_to();
    std::ofstream outFile (template_file_name + "0");
    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
        double pos;
        convert_idx_to_pos(i,pos);
        outFile << pos << " ";
        outFile << uSolutions_[i] << std::endl;
    } // writes file in two column format: x and u(x)
    outFile.close();


    // Apply step at every time step
    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
        apply_step();

        // write results to output file, need t's, x's, and u's
        outFile.open(template_file_name + std::to_string(t));
        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
            double pos;
            convert_idx_to_pos(i,pos);
            outFile << pos << " ";
            outFile << uSolutions_[i] << std::endl;
            // close file when done
        }
        outFile.close();

    }

}

//============================================================================================================
void SolutionProcedureFixLinWave::end_procedure(){
    // something important to finish writing the output file
    std::cout << "Done." << std::endl;
}

//============================================================================================================
void SolutionProcedureFixLinWave::apply_step(){

    // Apply individual step to solutions

    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // loops over entries, save the one involved in the boundary condition
    // boundary hard coded to LHS of wall

    double uSolnm1 = uSolutions_[0]; // Store solution of u_(i-1) at time n

    for (unsigned int i = 1; i < runtime_param_->get_space_iterations(); i++){
        // foward in time, backward in space

        // make sure to use u_(i-1) state at time n, not n+1, so have to save previous value

        uSolutions_[i] = uSolutions_[i] - runtime_param_->get_CFL()/2. * (uSolutions_[i+1] - uSolnm1)
                + std::abs(runtime_param_->get_CFL()/2.) * (uSolutions_[i+1] - 2 * uSolutions_[i] + uSolnm1);

        // write new position, will be i-1 at next pass
        uSolnm1 = uSolutions_[i];
    }

    // Need to write out solutions, but for now, leave that out
}

//============================================================================================================
void SolutionProcedureFixLinWave::set_boundary(){

    // hard code as LHSWall
    //std::string& myString;
    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
    solutionBoundaryCondition_ = new LHSWall;


    //TODO: should probably allow the user to go ahead and set what boundary condition to be used

}

//============================================================================================================
void SolutionProcedureFixLinWave::set_init_cond(){

    // hard code as wave type
    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
    //solutionInitialCondition_ = new StepWave;
    solutionInitialCondition_ = new SinWave;

    // TODO: should probably allow the user to go ahead and set the initial condition from the control file

}

//============================================================================================================
/*
 * Solution procedure for the diffusion equation
 */
//============================================================================================================

//============================================================================================================
void SolutionProcedureDiffusion::start_procedure(std::string& runtime_params, std::string& template_file_name ){

    // Type is known, go ahead and force it
    runtime_param_ = new RuntimeParametersDiffusion;
    // Take from file to set parameters
    runtime_param_->read_parameters_from_file(runtime_params);

    // Will take from input user to set boundary/initial condition type
    // Fow now, hard code it

    set_boundary();
    set_init_cond();

    // Set size of uSolutions_ vector
    uSolutions_.resize(runtime_param_->get_space_iterations());

    // Apply boundary condition, LHSWall
    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // For first step, need to apply the initial condition
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();
    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);

    // execute procedure
    procedure(template_file_name);

    // end procedure
    end_procedure();

}

//============================================================================================================
void SolutionProcedureDiffusion::convert_idx_to_pos(unsigned int idx, double& pos){

    // Convert from index to position
    // idx * dx + xo = current position
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();

    pos = (double)(idx) * dx + xo;

}

//============================================================================================================
void SolutionProcedureDiffusion::write_to_file()
{

}

//============================================================================================================
void SolutionProcedureDiffusion::procedure(std::string& template_file_name){

    // write first point -- initial condition onto file
    double to = runtime_param_->get_to();
    std::ofstream outFile (template_file_name + "0");
    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
        double pos;
        convert_idx_to_pos(i,pos);
        outFile << pos << " ";
        outFile << uSolutions_[i] << std::endl;
    } // writes file in two column format: x and u(x)
    outFile.close();

    // Apply step at every time step
    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
        apply_step();

        // write results to output file, need t's, x's, and u's
        outFile.open(template_file_name + std::to_string(t));
        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
            double pos;
            convert_idx_to_pos(i,pos);
            outFile << pos << " ";
            outFile << uSolutions_[i] << std::endl;
            // close file when done
        }
        outFile.close();

    }

}

//============================================================================================================
void SolutionProcedureDiffusion::end_procedure(){
    // something important to finish writing the output file
    std::cout << "Done." << std::endl;
}

//============================================================================================================
void SolutionProcedureDiffusion::apply_step(){

    // Apply individual step to solutions

    std::vector<double> myWallValues = runtime_param_->get_wall_vals();
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // loops over entries, save the one involved in the boundary condition
    // boundary hard coded to LHS of wall

    // need to store u_(i-1) at time n
    double uSolnm1 = uSolutions_[0];

    for (unsigned int i = 1; i < runtime_param_->get_space_iterations()-1; i++){ // ensure it won't go too far
        // foward in time, central in space
        uSolutions_[i] = (1.0 - 2.0 * runtime_param_->get_alpha()) * uSolutions_[i] +
                         runtime_param_->get_alpha() * (uSolutions_[i+1] + uSolnm1);
        // values at end points are handled by the specified boundary conditions

        // update previous state
        uSolnm1 = uSolutions_[i];
    }

    // Need to write out solutions, but for now, leave that out
}

//============================================================================================================
void SolutionProcedureDiffusion::set_boundary(){

    // hard code as LHSWall
    //std::string& myString;
    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
    solutionBoundaryCondition_ = new TwoWall;


    //TODO: should probably allow the user to go ahead and set what boundary condition to be used

}

//============================================================================================================
void SolutionProcedureDiffusion::set_init_cond(){

    // hard code as wave type
    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
    //solutionInitialCondition_ = new StepWave;
    solutionInitialCondition_ = new SinWave;

    // TODO: should probably allow the user to go ahead and set the initial condition from the control file

}

//============================================================================================================
/*
 * Solution procedure for inviscid burger equation
 */
//============================================================================================================

//============================================================================================================
void SolutionProcedureInviscid::start_procedure(std::string& runtime_params, std::string& template_file_name ){

    // Since type is already known, just hard code in the file reader
    runtime_param_ = new RuntimeParametersLinWave;

    runtime_param_->read_parameters_from_file(runtime_params);

    // Will take from input user to set boundary/initial condition type
    // Fow now, hard code it

    set_boundary();
    set_init_cond();

    // Set size of uSolutions_ vector
    uSolutions_.resize(runtime_param_->get_space_iterations());

    // Apply boundary condition, LHSWall
    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // For first step, need to apply the initial condition
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();
    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);

    // execute procedure
    procedure(template_file_name);

    // end procedure
    end_procedure();

}

//============================================================================================================
void SolutionProcedureInviscid::convert_idx_to_pos(unsigned int idx, double& pos){

    // Convert from index to position
    // idx * dx + xo = current position
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();

    pos = (double)(idx) * dx + xo;

}

//============================================================================================================
void SolutionProcedureInviscid::write_to_file()
{

}

//============================================================================================================
void SolutionProcedureInviscid::procedure(std::string& template_file_name){

    // write first point -- initial condition onto file
    double to = runtime_param_->get_to();
    std::ofstream outFile (template_file_name + "0");
    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
        double pos;
        convert_idx_to_pos(i,pos);
        outFile << pos << " ";
        outFile << uSolutions_[i] << std::endl;
    } // writes file in two column format: x and u(x)
    outFile.close();


    // Apply step at every time step
    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
        apply_step();

        // write results to output file, need t's, x's, and u's
        outFile.open(template_file_name + std::to_string(t));
        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
            double pos;
            convert_idx_to_pos(i,pos);
            outFile << pos << " ";
            outFile << uSolutions_[i] << std::endl;
            // close file when done
        }
        outFile.close();

    }

}

//============================================================================================================
void SolutionProcedureInviscid::end_procedure(){
    // something important to finish writing the output file
    std::cout << "Done." << std::endl;
}

//============================================================================================================
void SolutionProcedureInviscid::apply_step(){

    // Apply individual step to solutions

    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // loops over entries, save the one involved in the boundary condition
    // boundary hard coded to LHS of wall

    // Record u(i-1) at time step n
    double uSolnm1 = uSolutions_[0];

    for (unsigned int i = 1; i < runtime_param_->get_space_iterations(); i++){
        // foward in time, backward in space
        double CFL = uSolutions_[i] * runtime_param_->get_dt() / runtime_param_->get_dx();
        uSolutions_[i] = uSolutions_[i] - CFL *(uSolutions_[i] - uSolnm1);

        uSolnm1 = uSolutions_[i];
    }

    // Need to write out solutions, but for now, leave that out
}

//============================================================================================================
void SolutionProcedureInviscid::set_boundary(){

    // hard code as LHSWall
    //std::string& myString;
    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
    solutionBoundaryCondition_ = new LHSWall;


    //TODO: should probably allow the user to go ahead and set what boundary condition to be used

}

//============================================================================================================
void SolutionProcedureInviscid::set_init_cond(){

    // hard code as wave type
    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
    //solutionInitialCondition_ = new StepWave;
    solutionInitialCondition_ = new PositiveWave;

    // TODO: should probably allow the user to go ahead and set the initial condition from the control file

}



//============================================================================================================
/*
 * Solution procedure for full burger equation
 */
//============================================================================================================

//============================================================================================================
void SolutionProcedureBurger::start_procedure(std::string& runtime_params, std::string& template_file_name ){

    // Since type is already known, just hard code in the file reader
    runtime_param_ = new RuntimeParametersBurger;

    runtime_param_->read_parameters_from_file(runtime_params);

    // Will take from input user to set boundary/initial condition type
    // Fow now, hard code it

    set_boundary();
    set_init_cond();

    // Set size of uSolutions_ vector
    uSolutions_.resize(runtime_param_->get_space_iterations());

    // Apply boundary condition, LHSWall
    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // For first step, need to apply the initial condition
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();
    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);

    // execute procedure
    procedure(template_file_name);

    // end procedure
    end_procedure();

}

//============================================================================================================
void SolutionProcedureBurger::convert_idx_to_pos(unsigned int idx, double& pos){

    // Convert from index to position
    // idx * dx + xo = current position
    double dx = runtime_param_->get_dx();
    double xo = runtime_param_->get_xo();

    pos = (double)(idx) * dx + xo;

}

//============================================================================================================
void SolutionProcedureBurger::write_to_file()
{

}

//============================================================================================================
void SolutionProcedureBurger::procedure(std::string& template_file_name){

    // write first point -- initial condition onto file
    double to = runtime_param_->get_to();
    std::ofstream outFile (template_file_name + "0");
    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
        double pos;
        convert_idx_to_pos(i,pos);
        outFile << pos << " ";
        outFile << uSolutions_[i] << std::endl;
    } // writes file in two column format: x and u(x)
    outFile.close();


    // Apply step at every time step
    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
        apply_step();

        // write results to output file, need t's, x's, and u's
        outFile.open(template_file_name + std::to_string(t));
        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
            double pos;
            convert_idx_to_pos(i,pos);
            outFile << pos << " ";
            outFile << uSolutions_[i] << std::endl;
            // close file when done
        }
        outFile.close();

    }

}

//============================================================================================================
void SolutionProcedureBurger::end_procedure(){
    // something important to finish writing the output file
    std::cout << "Done." << std::endl;
}

//============================================================================================================
void SolutionProcedureBurger::apply_step(){

    // Apply individual step to solutions

    std::vector<double> myWallValues;
    myWallValues.push_back(runtime_param_->get_wall_value());
    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);

    // loops over entries, save the one involved in the boundary condition
    // boundary hard coded to LHS of wall
    uSolutions_[0] = 2.0;
    uSolutions_[uSolutions_.size()-1] = 2.544; // Need to remove hard code

    // store u(i-1) at time step n
    double uSolnm1 = uSolutions_[0];

    for (unsigned int i = 1; i < runtime_param_->get_space_iterations()-1; i++){
        // foward in time, backward in space
        double CFL = uSolutions_[i] * runtime_param_->get_dt() / runtime_param_->get_dx();
        double alpha = runtime_param_->get_alpha();

        uSolutions_[i] = uSolutions_[i] - CFL *(uSolutions_[i] - uSolnm1) +
                alpha * (uSolutions_[i+1] - 2 * uSolutions_[i] + uSolnm1);

        uSolnm1 = uSolutions_[i];
    }

    // Need to write out solutions, but for now, leave that out
}

//============================================================================================================
void SolutionProcedureBurger::set_boundary(){

    // hard code as LHSWall
    //std::string& myString;
    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
    solutionBoundaryCondition_ = new TwoWall;


    //TODO: should probably allow the user to go ahead and set what boundary condition to be used

}

//============================================================================================================
void SolutionProcedureBurger::set_init_cond(){

    // hard code as wave type
    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
    //solutionInitialCondition_ = new StepWave;
    solutionInitialCondition_ = new PositiveWave;

    // TODO: should probably allow the user to go ahead and set the initial condition from the control file

}

//============================================================================================================
/*
 * Solution procedure for linear advection equation
 */
//============================================================================================================

//============================================================================================================
void SolutionProcedureMultiDimAdvection::start_procedure(std::string& runtime_params, std::string& template_file_name ){

    // Since type is already known, just hard code in the file reader
    runtime_args_ = new RuntimeParamMultiDim;

    runtime_args_->read_parameters_from_file(runtime_params);

    // Will take from input user to set boundary/initial condition type
    // Fow now, hard code it

    set_boundary();
    set_init_cond();

    // set size of vector

    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));

    // For first step, need to apply the initial condition
    double dx = runtime_args_->get_dx();
    double xo = runtime_args_->get_xo();
    multiSolutionInitialCondition_->apply_initial_cond(multiUSolutions_, *(runtime_args_));

    // execute procedure
    procedure(template_file_name);

    // end procedure
    end_procedure();

}

//============================================================================================================
void SolutionProcedureMultiDimAdvection::convert_idx_to_pos_x(unsigned int idx, double& pos){

    // Convert from index to position
    // idx * dx + xo = current position
    double dx = runtime_args_->get_dx();
    double xo = runtime_args_->get_xo();

    pos = (double)(idx) * dx + xo;

}
//============================================================================================================
void SolutionProcedureMultiDimAdvection::convert_idx_to_pos_y(unsigned int idx, double& pos){

    // Convert from index to position
    // idx * dx + xo = current position
    double dy = runtime_args_->get_dy();
    double yo = runtime_args_->get_yo();

    pos = (double)(idx) * dy + yo;

}
//============================================================================================================
void SolutionProcedureMultiDimAdvection::write_to_file()
{

}

//============================================================================================================
void SolutionProcedureMultiDimAdvection::procedure(std::string& template_file_name){

    // write first point -- initial condition onto file
    double to = runtime_args_->get_to();
    std::ofstream outFile (template_file_name + "0");

    // Needed nested for loops
    for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
            //write to file (x,y) pair, with u(x,y) at current time step
            double posx;
            convert_idx_to_pos_x(i, posx);
            double posy;
            convert_idx_to_pos_y(j, posy);
            outFile << posx << "     ";
            outFile << posy << "     ";
            outFile << multiUSolutions_[i][j] << std::endl;
        }
    }

    outFile.close();


    // Apply step at every time step
    for(unsigned int t = 1; t < runtime_args_->get_time_iterations(); ++t){
        apply_step();

        // write results to output file, need x's, y's, and t's
        outFile.open(template_file_name + std::to_string(t));

        // I hate the idea of a triply nested for loop, but this just prints out data
        for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
            for (unsigned int  j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
                //write to file (x,y) pair, with u(x,y) at current time step
                double posx;
                convert_idx_to_pos_x(i,posx);
                double posy;
                convert_idx_to_pos_y(j,posy);
                outFile << posx << "     ";
                outFile << posy << "     ";
                outFile << multiUSolutions_[i][j] << std::endl;
            }
        }

        outFile.close();

    }

}

//============================================================================================================
void SolutionProcedureMultiDimAdvection::end_procedure(){
    // something important to finish writing the output file
    std::cout << "Done." << std::endl;
}

//============================================================================================================
void SolutionProcedureMultiDimAdvection::apply_step(){

    // Apply individual step to solutions

    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));

    // loops over entries, save the one involved in the boundary condition
    // boundary hard coded to LHS of wall


    double c = runtime_args_->get_c();
    double dx = runtime_args_->get_dx();
    double dy = runtime_args_->get_dy();
    double dt = runtime_args_->get_dt();
    // start by looping over x-variable (space)

    // save previous column in holder
    std::vector<double> previous_col (runtime_args_->get_y_iterations(), 0.0);
    for (unsigned int i = 0 ; i < runtime_args_->get_y_iterations(); ++i){
        previous_col[i] = multiUSolutions_[0][i];
    }

    for (unsigned int  i = 1 ; i < runtime_args_->get_x_iterations()-1; ++i){
        double uSolnym1 = multiUSolutions_[i][0];
        for (unsigned int j = 1 ; j < runtime_args_->get_y_iterations()-1; ++j){
            multiUSolutions_[i][j] = multiUSolutions_[i][j] -
                                     c*dt/dx * (multiUSolutions_[i][j] - previous_col[j]) -
                                     c*dt/dy * (multiUSolutions_[i][j] - uSolnym1);
            uSolnym1 = multiUSolutions_[i][j];
        }

        //rewrite the previous_col
        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
            previous_col[j] = multiUSolutions_[i][j];
        }

    }
     // Need to save the (n) state so as to not use (n+1) in solutions

    // Need to write out solutions, but for now, leave that out
}

//============================================================================================================
void SolutionProcedureMultiDimAdvection::set_boundary(){

    multiSolutionBoundaryCondition_ = new BoxBoundaryCondition;


    //TODO: should probably allow the user to go ahead and set what boundary condition to be used

}

//============================================================================================================
void SolutionProcedureMultiDimAdvection::set_init_cond(){

    // hard code as wave type
    multiSolutionInitialCondition_ = new Curvilinear;

    // TODO: should probably allow the user to go ahead and set the initial condition from the control file

}
