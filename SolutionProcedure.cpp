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

    solutionInitialCondition_->set_args(runtime->return_args());
    std::string condition;
    condition = runtime->return_cond();

    solutionInitialCondition_ = solutionInitialCondition_->make_initial_condition(condition);

    solutionInitialCondition_->fluidEquation_ = solutionInitialCondition_->make_fluid_equation(whichPDE);

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


////============================================================================================================
///*
// * Solution procedure for linear wave equation, where c may be positive or negative
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureFixLinWave::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Since type is already known, just hard code in the file reader
//    runtime_param_ = new RuntimeParametersLinWave;
//
//    runtime_param_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // Set size of uSolutions_ vector
//    uSolutions_.resize(runtime_param_->get_space_iterations());
//
//    // Apply boundary condition, LHSWall
//    std::vector<double> myWallValues;
//    myWallValues.push_back(runtime_param_->get_wall_value());
//    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);
//
//    // For first step, need to apply the initial condition
//    double dx = runtime_param_->get_dx();
//    double xo = runtime_param_->get_xo();
//    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureFixLinWave::convert_idx_to_pos(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_param_->get_dx();
//    double xo = runtime_param_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
//
////============================================================================================================
//void SolutionProcedureFixLinWave::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureFixLinWave::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_param_->get_to();
//    std::ofstream outFile (template_file_name + "0");
//    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
//        double pos;
//        convert_idx_to_pos(i,pos);
//        outFile << pos << " ";
//        outFile << uSolutions_[i] << std::endl;
//    } // writes file in two column format: x and u(x)
//    outFile.close();
//
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need t's, x's, and u's
//        outFile.open(template_file_name + std::to_string(t));
//        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
//            double pos;
//            convert_idx_to_pos(i,pos);
//            outFile << pos << " ";
//            outFile << uSolutions_[i] << std::endl;
//            // close file when done
//        }
//        outFile.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureFixLinWave::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureFixLinWave::apply_step(){
//
//    // Apply individual step to solutions
//
//    std::vector<double> myWallValues;
//    myWallValues.push_back(runtime_param_->get_wall_value());
//    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//
//    double uSolnm1 = uSolutions_[0]; // Store solution of u_(i-1) at time n
//
//    for (unsigned int i = 1; i < runtime_param_->get_space_iterations(); i++){
//        // foward in time, backward in space
//
//        // make sure to use u_(i-1) state at time n, not n+1, so have to save previous value
//
//        uSolutions_[i] = uSolutions_[i] - runtime_param_->get_CFL()/2. * (uSolutions_[i+1] - uSolnm1)
//                + std::abs(runtime_param_->get_CFL()/2.) * (uSolutions_[i+1] - 2 * uSolutions_[i] + uSolnm1);
//
//        // write new position, will be i-1 at next pass
//        uSolnm1 = uSolutions_[i];
//    }
//
//    // Need to write out solutions, but for now, leave that out
//}
//
////============================================================================================================
//void SolutionProcedureFixLinWave::set_boundary(){
//
//    // hard code as LHSWall
//    //std::string& myString;
//    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
//    solutionBoundaryCondition_ = new LHSWall;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureFixLinWave::set_init_cond(){
//
//    // hard code as wave type
//    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
//    //solutionInitialCondition_ = new StepWave;
//    solutionInitialCondition_ = new SinWave;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}
//
////============================================================================================================
///*
// * Solution procedure for the diffusion equation
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureDiffusion::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Type is known, go ahead and force it
//    runtime_param_ = new RuntimeParametersDiffusion;
//    // Take from file to set parameters
//    runtime_param_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // Set size of uSolutions_ vector
//    uSolutions_.resize(runtime_param_->get_space_iterations());
//
//    // Apply boundary condition, LHSWall
//    std::vector<double> myWallValues;
//    myWallValues.push_back(runtime_param_->get_wall_value());
//    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);
//
//    // For first step, need to apply the initial condition
//    double dx = runtime_param_->get_dx();
//    double xo = runtime_param_->get_xo();
//    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureDiffusion::convert_idx_to_pos(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_param_->get_dx();
//    double xo = runtime_param_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
//
////============================================================================================================
//void SolutionProcedureDiffusion::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureDiffusion::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_param_->get_to();
//    std::ofstream outFile (template_file_name + "0");
//    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
//        double pos;
//        convert_idx_to_pos(i,pos);
//        outFile << pos << " ";
//        outFile << uSolutions_[i] << std::endl;
//    } // writes file in two column format: x and u(x)
//    outFile.close();
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need t's, x's, and u's
//        outFile.open(template_file_name + std::to_string(t));
//        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
//            double pos;
//            convert_idx_to_pos(i,pos);
//            outFile << pos << " ";
//            outFile << uSolutions_[i] << std::endl;
//            // close file when done
//        }
//        outFile.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureDiffusion::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureDiffusion::apply_step(){
//
//    // Apply individual step to solutions
//
//    std::vector<double> myWallValues = runtime_param_->get_wall_vals();
//    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//
//    // need to store u_(i-1) at time n
//    double uSolnm1 = uSolutions_[0];
//
//    for (unsigned int i = 1; i < runtime_param_->get_space_iterations()-1; i++){ // ensure it won't go too far
//        // foward in time, central in space
//        uSolutions_[i] = (1.0 - 2.0 * runtime_param_->get_alpha()) * uSolutions_[i] +
//                         runtime_param_->get_alpha() * (uSolutions_[i+1] + uSolnm1);
//        // values at end points are handled by the specified boundary conditions
//
//        // update previous state
//        uSolnm1 = uSolutions_[i];
//    }
//
//    // Need to write out solutions, but for now, leave that out
//}
//
////============================================================================================================
//void SolutionProcedureDiffusion::set_boundary(){
//
//    // hard code as LHSWall
//    //std::string& myString;
//    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
//    solutionBoundaryCondition_ = new TwoWall;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureDiffusion::set_init_cond(){
//
//    // hard code as wave type
//    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
//    //solutionInitialCondition_ = new StepWave;
//    solutionInitialCondition_ = new SinWave;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}
//
////============================================================================================================
///*
// * Solution procedure for inviscid burger equation
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureInviscid::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Since type is already known, just hard code in the file reader
//    runtime_param_ = new RuntimeParametersLinWave;
//
//    runtime_param_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // Set size of uSolutions_ vector
//    uSolutions_.resize(runtime_param_->get_space_iterations());
//
//    // Apply boundary condition, LHSWall
//    std::vector<double> myWallValues;
//    myWallValues.push_back(runtime_param_->get_wall_value());
//    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);
//
//    // For first step, need to apply the initial condition
//    double dx = runtime_param_->get_dx();
//    double xo = runtime_param_->get_xo();
//    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureInviscid::convert_idx_to_pos(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_param_->get_dx();
//    double xo = runtime_param_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
//
////============================================================================================================
//void SolutionProcedureInviscid::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureInviscid::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_param_->get_to();
//    std::ofstream outFile (template_file_name + "0");
//    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
//        double pos;
//        convert_idx_to_pos(i,pos);
//        outFile << pos << " ";
//        outFile << uSolutions_[i] << std::endl;
//    } // writes file in two column format: x and u(x)
//    outFile.close();
//
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need t's, x's, and u's
//        outFile.open(template_file_name + std::to_string(t));
//        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
//            double pos;
//            convert_idx_to_pos(i,pos);
//            outFile << pos << " ";
//            outFile << uSolutions_[i] << std::endl;
//            // close file when done
//        }
//        outFile.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureInviscid::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureInviscid::apply_step(){
//
//    // Apply individual step to solutions
//
//    std::vector<double> myWallValues;
//    myWallValues.push_back(runtime_param_->get_wall_value());
//    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//
//    // Record u(i-1) at time step n
//    double uSolnm1 = uSolutions_[0];
//
//    for (unsigned int i = 1; i < runtime_param_->get_space_iterations(); i++){
//        // foward in time, backward in space
//        double CFL = uSolutions_[i] * runtime_param_->get_dt() / runtime_param_->get_dx();
//        uSolutions_[i] = uSolutions_[i] - CFL *(uSolutions_[i] - uSolnm1);
//
//        uSolnm1 = uSolutions_[i];
//    }
//
//    // Need to write out solutions, but for now, leave that out
//}
//
////============================================================================================================
//void SolutionProcedureInviscid::set_boundary(){
//
//    // hard code as LHSWall
//    //std::string& myString;
//    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
//    solutionBoundaryCondition_ = new LHSWall;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureInviscid::set_init_cond(){
//
//    // hard code as wave type
//    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
//    //solutionInitialCondition_ = new StepWave;
//    solutionInitialCondition_ = new PositiveWave;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}
//
//
//
////============================================================================================================
///*
// * Solution procedure for full burger equation
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureBurger::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Since type is already known, just hard code in the file reader
//    runtime_param_ = new RuntimeParametersBurger;
//
//    runtime_param_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // Set size of uSolutions_ vector
//    uSolutions_.resize(runtime_param_->get_space_iterations());
//
//    // Apply boundary condition, LHSWall
//    std::vector<double> myWallValues;
//    myWallValues.push_back(runtime_param_->get_wall_value());
//    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);
//
//    // For first step, need to apply the initial condition
//    double dx = runtime_param_->get_dx();
//    double xo = runtime_param_->get_xo();
//    solutionInitialCondition_->apply_initial_cond(uSolutions_, dx, xo);
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureBurger::convert_idx_to_pos(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_param_->get_dx();
//    double xo = runtime_param_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
//
////============================================================================================================
//void SolutionProcedureBurger::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureBurger::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_param_->get_to();
//    std::ofstream outFile (template_file_name + "0");
//    for(unsigned int i = 0 ; i < uSolutions_.size(); ++i){
//        double pos;
//        convert_idx_to_pos(i,pos);
//        outFile << pos << " ";
//        outFile << uSolutions_[i] << std::endl;
//    } // writes file in two column format: x and u(x)
//    outFile.close();
//
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_param_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need t's, x's, and u's
//        outFile.open(template_file_name + std::to_string(t));
//        for (unsigned int i = 0 ; i < uSolutions_.size(); ++i){
//            double pos;
//            convert_idx_to_pos(i,pos);
//            outFile << pos << " ";
//            outFile << uSolutions_[i] << std::endl;
//            // close file when done
//        }
//        outFile.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureBurger::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureBurger::apply_step(){
//
//    // Apply individual step to solutions
//
//    std::vector<double> myWallValues;
//    myWallValues.push_back(runtime_param_->get_wall_value());
//    solutionBoundaryCondition_->enforce_boundary_conditions(uSolutions_, myWallValues);
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//    uSolutions_[0] = 2.0;
//    uSolutions_[uSolutions_.size()-1] = 2.544; // Need to remove hard code
//
//    // store u(i-1) at time step n
//    double uSolnm1 = uSolutions_[0];
//
//    for (unsigned int i = 1; i < runtime_param_->get_space_iterations()-1; i++){
//        // foward in time, backward in space
//        double CFL = uSolutions_[i] * runtime_param_->get_dt() / runtime_param_->get_dx();
//        double alpha = runtime_param_->get_alpha();
//
//        uSolutions_[i] = uSolutions_[i] - CFL *(uSolutions_[i] - uSolnm1) +
//                alpha * (uSolutions_[i+1] - 2 * uSolutions_[i] + uSolnm1);
//
//        uSolnm1 = uSolutions_[i];
//    }
//
//    // Need to write out solutions, but for now, leave that out
//}
//
////============================================================================================================
//void SolutionProcedureBurger::set_boundary(){
//
//    // hard code as LHSWall
//    //std::string& myString;
//    //solutionBoundaryCondition_ = BoundaryCondition::make_boundary_condition(myString);
//    solutionBoundaryCondition_ = new TwoWall;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureBurger::set_init_cond(){
//
//    // hard code as wave type
//    //solutionInitialCondition_ = BoundaryCondition::make_init_cond();
//    //solutionInitialCondition_ = new StepWave;
//    solutionInitialCondition_ = new PositiveWave;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}
//
////============================================================================================================
///*
// * Solution procedure for linear advection equation
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Since type is already known, just hard code in the file reader
//    runtime_args_ = new RuntimeParamMultiDim;
//
//    runtime_args_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // set size of vector
//    unsigned int XSIZE = runtime_args_->get_x_iterations();
//    unsigned int YSIZE = runtime_args_->get_y_iterations();
//    multiUSolutions_.resize(XSIZE, std::vector<double>(YSIZE));
//
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//
//    // For first step, need to apply the initial condition
//    double dx = runtime_args_->get_dx();
//    double xo = runtime_args_->get_xo();
//    multiSolutionInitialCondition_->apply_initial_cond(multiUSolutions_, *(runtime_args_));
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::convert_idx_to_pos_x(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_args_->get_dx();
//    double xo = runtime_args_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::convert_idx_to_pos_y(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dy = runtime_args_->get_dy();
//    double yo = runtime_args_->get_yo();
//
//    pos = (double)(idx) * dy + yo;
//
//}
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_args_->get_to();
//    std::ofstream outFile (template_file_name + "0");
//    std::cout.fill(' ');
//
//    // Needed nested for loops
//    for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            //write to file (x,y) pair, with u(x,y) at current time step
//            double posx;
//            convert_idx_to_pos_x(i, posx);
//            double posy;
//            convert_idx_to_pos_y(j, posy);
//            outFile << std::setw(4) << posx << "    ";
//            outFile << posy << "     ";
//            outFile << multiUSolutions_[i][j] << std::endl;
//        }
//    }
//
//    outFile.close();
//
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_args_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need x's, y's, and t's
//        outFile.open(template_file_name + std::to_string(t));
//
//        // I hate the idea of a triply nested for loop, but this just prints out data
//        for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//            for (unsigned int  j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//                //write to file (x,y) pair, with u(x,y) at current time step
//                double posx;
//                convert_idx_to_pos_x(i,posx);
//                double posy;
//                convert_idx_to_pos_y(j,posy);
//                outFile << std::setw(4) << posx << "    ";
//                outFile << posy << "    ";
//                outFile << multiUSolutions_[i][j] << std::endl;
//            }
//        }
//
//        outFile.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::apply_step(){
//
//    // Apply individual step to solutions
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//
//
//    double c = runtime_args_->get_c();
//    double dx = runtime_args_->get_dx();
//    double dy = runtime_args_->get_dy();
//    double dt = runtime_args_->get_dt();
//    // start by looping over x-variable (space)
//
//    // save previous column in holder
//    std::vector<double> previous_col (runtime_args_->get_y_iterations(), 0.0);
//    for (unsigned int i = 0 ; i < runtime_args_->get_y_iterations(); ++i){
//        previous_col[i] = multiUSolutions_[0][i];
//    }
//
//    for (unsigned int  i = 1 ; i < runtime_args_->get_x_iterations()-2; ++i){
//        double uSolnym1 = multiUSolutions_[i][0];
//        for (unsigned int j = 1 ; j < runtime_args_->get_y_iterations()-2; ++j){
//            double uij = multiUSolutions_[i][j];
//            double uip1j = multiUSolutions_[i+1][j];
//            double uim1j = previous_col[j];
//            double uijp1 = multiUSolutions_[i][j+1];
//            double uijm1 = uSolnym1;
//            double fx = 0.5*dt/dx;
//            double fy = 0.5*dt/dy;
//
//            // Change to upshift/downshift equation form
//            multiUSolutions_[i][j] = uij - c * fx * (uip1j - uim1j) + fabs(c) * fx * (uip1j - 2. * uij + uim1j)
//                    - c * fy * (uijp1 - uijm1) + fabs(c) * fy * (uijp1 - 2. * uij + uijm1);
//            uSolnym1 = multiUSolutions_[i][j];
//        }
//
//        //rewrite the previous_col
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            previous_col[j] = multiUSolutions_[i][j];
//        }
//
//    }
//     // Need to save the (n) state so as to not use (n+1) in solutions
//
//    // Need to write out solutions, but for now, leave that out
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::set_boundary(){
//
//    multiSolutionBoundaryCondition_ = new BoxBoundaryCondition;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimAdvection::set_init_cond(){
//
//    // hard code as wave type
//    multiSolutionInitialCondition_ = new Curvilinear;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}
//
////============================================================================================================
///*
// * Solution procedure for multi dimension non linear coupled advective equation
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Since type is already known, just hard code in the file reader
//    runtime_args_ = new RuntimeParamMultiDim;
//
//    runtime_args_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // set size of vector
//    unsigned int XSIZE = runtime_args_->get_x_iterations();
//    unsigned int YSIZE = runtime_args_->get_y_iterations();
//    multiUSolutions_.resize(XSIZE, std::vector<double>(YSIZE));
//    multiVSolutions_.resize(XSIZE, std::vector<double>(YSIZE));
//
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiVSolutions_, *(runtime_args_));
//
//    // For first step, need to apply the initial condition
//    multiSolutionInitialCondition_->apply_initial_cond(multiUSolutions_, *(runtime_args_));
//    multiSolutionInitialCondition_->apply_initial_cond(multiVSolutions_, *(runtime_args_));
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::convert_idx_to_pos_x(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_args_->get_dx();
//    double xo = runtime_args_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::convert_idx_to_pos_y(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dy = runtime_args_->get_dy();
//    double yo = runtime_args_->get_yo();
//
//    pos = (double)(idx) * dy + yo;
//
//}
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_args_->get_to();
//    std::ofstream outFileU (template_file_name + "U0");
//    std::ofstream outFileV (template_file_name + "V0");
//
//    // Needed nested for loops
//    for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            //write to file (x,y) pair, with u(x,y) at current time step
//            double posx;
//            convert_idx_to_pos_x(i, posx);
//            double posy;
//            convert_idx_to_pos_y(j, posy);
//            outFileU << std::setw(4) << posx << "     ";
//            outFileU << posy << "     ";
//            outFileU << multiUSolutions_[i][j] << std::endl;
//            outFileV << std::setw(4) << posx << "     ";
//            outFileV << posy << "     ";
//            outFileV << multiVSolutions_[i][j] << std::endl;
//
//        }
//    }
//
//    outFileU.close();
//    outFileV.close();
//
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_args_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need x's, y's, and t's
//        outFileU.open(template_file_name + "U" + std::to_string(t));// + "U-" + std::to_string(currentTime_));
//        outFileV.open(template_file_name + "V" + std::to_string(t));// + "V-" + std::to_string(currentTime_));
//
//        // I hate the idea of a triply nested for loop, but this just prints out data
//        for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//            for (unsigned int  j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//                //write to file (x,y) pair, with u(x,y) at current time step
//                double posx;
//                convert_idx_to_pos_x(i,posx);
//                double posy;
//                convert_idx_to_pos_y(j,posy);
//                outFileU << std::setw(4) << posx << "     ";
//                outFileU << posy << "     ";
//                outFileU << multiUSolutions_[i][j] << std::endl;
//                outFileV << std::setw(4) << posx << "     ";
//                outFileV << posy << "     ";
//                outFileV << multiUSolutions_[i][j] << std::endl;
//            }
//        }
//
//        outFileU.close();
//        outFileV.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::apply_step(){
//
//    // Apply individual step to solutions
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//
//
//    double c = runtime_args_->get_c(); // No longer wave speed, but maximum allowable CFL number
//    double dx = runtime_args_->get_dx();
//    double dy = runtime_args_->get_dy();
//    double dt;
//
//    /*
//     * calculate current dt needed, new current time after update
//     */
//
//    // Find the maximum U and V, this is actually costly
//    std::vector<double> maxUs;
//    std::vector<double> maxVs;
//    for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//        maxUs.push_back(*std::max_element(std::begin(multiUSolutions_[i]), std::end(multiUSolutions_[i])));
//        maxVs.push_back(*std::max_element(std::begin(multiVSolutions_[i]), std::end(multiVSolutions_[i])));
//    }
//    double maxU = *std::max_element(std::begin(maxUs), std::end(maxUs));
//    double maxV = *std::max_element(std::begin(maxVs), std::end(maxVs));
//
//    double maxC = std::max(maxU, maxV); // Largest value in wave speed
//
//    /*
//     * Have Udt/dx + Udt/dy = SF -> dt = SF/(U/dx + U/dy)
//     */
//
//    dt = c/(maxC/dx + maxC/dy); // Recall c: is in this case the max CFL (safety), likely 0.9 for safety
//
//    // start by looping over x-variable (space)
//
//    // save previous column in holder
//    std::vector<double> previous_colU (runtime_args_->get_y_iterations(), 0.0);
//    std::vector<double> previous_colV (runtime_args_->get_y_iterations(), 0.0);
//    for (unsigned int i = 0 ; i < runtime_args_->get_y_iterations(); ++i){
//        previous_colU[i] = multiUSolutions_[0][i];
//        previous_colV[i] = multiVSolutions_[0][i];
//    }
//
//    for (unsigned int  i = 1 ; i < runtime_args_->get_x_iterations()-1; ++i){
//        double uSolnym1 = multiUSolutions_[i][0];
//        double vSolynm1 = multiVSolutions_[i][0];
//        for (unsigned int j = 1 ; j < runtime_args_->get_y_iterations()-1; ++j){
//            // Solve U, V into temporary holders first
//            double temp_u_at_ij;
//            double temp_v_at_ij;
//
//            // For easier reference
//            double uij = multiUSolutions_[i][j];
//            double uip1j = multiUSolutions_[i+1][j];
//            double uim1j = previous_colU[j];
//            double uijp1 = multiUSolutions_[i][j+1];
//            double uijm1 = uSolnym1;
//            double vij = multiVSolutions_[i][j];
//            double vip1j = multiVSolutions_[i+1][j];
//            double vim1j = previous_colV[j];
//            double vijp1 = multiVSolutions_[i][j+1];
//            double vijm1 = vSolynm1;
//            double fx = 0.5*dt/dx;
//            double fy = 0.5*dt/dy;
//
//            temp_u_at_ij = uij - uij * fx * (uip1j - uim1j) + fabs(uij) * fx * (uip1j - 2.*uij + uim1j)
//                    - vij * fy * (uijp1 - uijm1) + fabs(vij) * fy * (uijp1 - 2.*uij + uijm1);
//
//            temp_v_at_ij = vij - uij * fx * (vip1j - vim1j) + fabs(uij) * fx * (vip1j - 2.*vij + vim1j)
//                    - vij * fy * (vijp1 - vijm1) + fabs(vij) * fy * (vijp1 - 2.*vij + vijm1);
//
//            // Write U, V into new form
//            multiUSolutions_[i][j] = temp_u_at_ij;
//            multiVSolutions_[i][j] = temp_v_at_ij;
//
//            uSolnym1 = multiUSolutions_[i][j];
//            vSolynm1 = multiVSolutions_[i][j];
//
//        }
//
//        //rewrite the previous_col
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            previous_colU[j] = multiUSolutions_[i][j];
//            previous_colV[j] = multiVSolutions_[i][j];
//        }
//
//    }
//
//    // Add to extra time
//    currentTime_ += dt;
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::set_boundary(){
//
//    multiSolutionBoundaryCondition_ = new BoxBoundaryCondition;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimNonLinAdvEqn::set_init_cond(){
//
//    // hard code as wave type
//    multiSolutionInitialCondition_ = new Curvilinear;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}
//
////============================================================================================================
///*
// * Solution procedure for multidimensional diffusion equation
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Since type is already known, just hard code in the file reader
//    runtime_args_ = new RuntimeParamMultiDim;
//
//    runtime_args_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // set size of vector
//    unsigned int XSIZE = runtime_args_->get_x_iterations();
//    unsigned int YSIZE = runtime_args_->get_y_iterations();
//    multiUSolutions_.resize(XSIZE, std::vector<double>(YSIZE));
//
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//
//    // For first step, need to apply the initial condition
//    double dx = runtime_args_->get_dx();
//    double xo = runtime_args_->get_xo();
//    multiSolutionInitialCondition_->apply_initial_cond(multiUSolutions_, *(runtime_args_));
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::convert_idx_to_pos_x(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_args_->get_dx();
//    double xo = runtime_args_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::convert_idx_to_pos_y(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dy = runtime_args_->get_dy();
//    double yo = runtime_args_->get_yo();
//
//    pos = (double)(idx) * dy + yo;
//
//}
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_args_->get_to();
//    std::ofstream outFile (template_file_name + "0");
//    std::cout.fill(' ');
//
//    // Needed nested for loops
//    for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            //write to file (x,y) pair, with u(x,y) at current time step
//            double posx;
//            convert_idx_to_pos_x(i, posx);
//            double posy;
//            convert_idx_to_pos_y(j, posy);
//            outFile << std::setw(4) << posx << "     ";
//            outFile << posy << "     ";
//            outFile << multiUSolutions_[i][j] << std::endl;
//        }
//    }
//
//    outFile.close();
//
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_args_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need x's, y's, and t's
//        outFile.open(template_file_name + std::to_string(t));
//
//        // I hate the idea of a triply nested for loop, but this just prints out data
//        for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//            for (unsigned int  j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//                //write to file (x,y) pair, with u(x,y) at current time step
//                double posx;
//                convert_idx_to_pos_x(i,posx);
//                double posy;
//                convert_idx_to_pos_y(j,posy);
//                outFile << std::setw(4) << posx << "     ";
//                outFile << posy << "     ";
//                outFile << multiUSolutions_[i][j] << std::endl;
//            }
//        }
//
//        outFile.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::apply_step(){
//
//    // Apply individual step to solutions
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//
//
//    double c = runtime_args_->get_c(); // In this case, the nu coefficient in the equation
//    double dx = runtime_args_->get_dx();
//    double dy = runtime_args_->get_dy();
//    double dt = runtime_args_->get_dt();
//
//    // Assume the user has specified a reasonable set of conditions
//    // Else, throw a warning to the user that this has occured.
//
//    // start by looping over x-variable (space)
//
//    // save previous column in holder
//    std::vector<double> previous_col (runtime_args_->get_y_iterations(), 0.0);
//    for (unsigned int i = 0 ; i < runtime_args_->get_y_iterations(); ++i){
//        previous_col[i] = multiUSolutions_[0][i];
//    }
//
//    for (unsigned int  i = 1 ; i < runtime_args_->get_x_iterations()-1; ++i){
//        double uSolnym1 = multiUSolutions_[i][0];
//        for (unsigned int j = 1 ; j < runtime_args_->get_y_iterations()-1; ++j){
//            double uij = multiUSolutions_[i][j];
//            double uip1j = multiUSolutions_[i+1][j];
//            double uim1j = previous_col[j];
//            double uijp1 = multiUSolutions_[i][j+1];
//            double uijm1 = uSolnym1;
//            double fx = c*dt/(dx*dx);
//            double fy = c*dt/(dy*dy);
//
//            multiUSolutions_[i][j] = uij + fx * (uip1j - 2. * uij + uim1j)
//                    + fy * (uijp1 - 2. * uij + uijm1);
//
//
//            uSolnym1 = multiUSolutions_[i][j];
//        }
//
//        //rewrite the previous_col
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            previous_col[j] = multiUSolutions_[i][j];
//        }
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::set_boundary(){
//
//    multiSolutionBoundaryCondition_ = new BoxBoundaryCondition;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::set_init_cond(){
//
//    // hard code as wave type
//    multiSolutionInitialCondition_ = new Curvilinear;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}
//
////============================================================================================================
///*
// * Solution procedure for multidimensional burger equation
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureMultiDimBurger::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Since type is already known, just hard code in the file reader
//    runtime_args_ = new RuntimeParamMultiDim;
//
//    runtime_args_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // set size of vector
//    unsigned int XSIZE = runtime_args_->get_x_iterations();
//    unsigned int YSIZE = runtime_args_->get_y_iterations();
//    multiUSolutions_.resize(XSIZE, std::vector<double>(YSIZE));
//    multiVSolutions_.resize(XSIZE, std::vector<double>(YSIZE));
//
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiVSolutions_, *(runtime_args_));
//
//    // For first step, need to apply the initial condition
//    multiSolutionInitialCondition_->apply_initial_cond(multiUSolutions_, *(runtime_args_));
//    multiSolutionInitialCondition_->apply_initial_cond(multiVSolutions_, *(runtime_args_));
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimBurger::convert_idx_to_pos_x(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_args_->get_dx();
//    double xo = runtime_args_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
////============================================================================================================
//void SolutionProcedureMultiDimBurger::convert_idx_to_pos_y(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dy = runtime_args_->get_dy();
//    double yo = runtime_args_->get_yo();
//
//    pos = (double)(idx) * dy + yo;
//
//}
////============================================================================================================
//void SolutionProcedureMultiDimBurger::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimBurger::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_args_->get_to();
//    std::ofstream outFileU (template_file_name + "U-0");
//    std::ofstream outFileV (template_file_name + "V-0");
//
//    // Needed nested for loops
//    for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            //write to file (x,y) pair, with u(x,y) at current time step
//            double posx;
//            convert_idx_to_pos_x(i, posx);
//            double posy;
//            convert_idx_to_pos_y(j, posy);
//            outFileU << std::setw(4) << posx << "     ";
//            outFileU << posy << "     ";
//            outFileU << multiUSolutions_[i][j] << std::endl;
//            outFileV << std::setw(4) << posx << "     ";
//            outFileV << posy << "     ";
//            outFileV << multiVSolutions_[i][j] << std::endl;
//
//        }
//    }
//
//    outFileU.close();
//    outFileV.close();
//
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_args_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need x's, y's, and t's
//        outFileU.open(template_file_name + std::to_string(currentTime_));
//        outFileV.open(template_file_name + std::to_string(currentTime_));
//
//        // I hate the idea of a triply nested for loop, but this just prints out data
//        for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//            for (unsigned int  j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//                //write to file (x,y) pair, with u(x,y) at current time step
//                double posx;
//                convert_idx_to_pos_x(i,posx);
//                double posy;
//                convert_idx_to_pos_y(j,posy);
//                outFileU << std::setw(4) << posx << "     ";
//                outFileU << posy << "     ";
//                outFileU << multiUSolutions_[i][j] << std::endl;
//                outFileV << std::setw(4) << posx << "     ";
//                outFileV << posy << "     ";
//                outFileV << multiUSolutions_[i][j] << std::endl;
//            }
//        }
//
//        outFileU.close();
//        outFileV.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimBurger::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimBurger::apply_step(){
//
//    // Apply individual step to solutions
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//
//
//    double c = runtime_args_->get_c(); // No longer wave speed, but maximum allowable CFL number
//    double dx = runtime_args_->get_dx();
//    double dy = runtime_args_->get_dy();
//    double nu = runtime_args_->get_tf(); // No longer tf, but nu on the diffusion portion
//    double dt;
//
//
//    /*
//     * calculate current dt needed, new current time after update
//     */
//
//    // Find the maximum U and V, this is actually costly
//    std::vector<double> maxUs;
//    std::vector<double> maxVs;
//    for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//        maxUs.push_back(*std::max_element(std::begin(multiUSolutions_[i]), std::end(multiUSolutions_[i])));
//        maxVs.push_back(*std::max_element(std::begin(multiVSolutions_[i]), std::end(multiVSolutions_[i])));
//    }
//    double maxU = *std::max_element(std::begin(maxUs), std::end(maxUs));
//    double maxV = *std::max_element(std::begin(maxVs), std::end(maxVs));
//
//    double maxC = std::max(maxU, maxV); // Largest value in wave speed
//
//    /*
//     * Have Udt/dx + Udt/dy = SF -> dt = SF/(U/dx + U/dy)
//     */
//    double dt_adv;
//    dt_adv = c/(maxC/dx + maxC/dy); // Recall c: is in this case the max CFL (safety), likely 0.9 for safety
//
//    // dt also needs to be calculated another means to ensure diffusive part is O.K.
//    double dt_diff;
//    double den = nu * (1./dx/dx + 1./dy/dy);
//    dt_diff = 0.5*c/den;
//
//    // Take minimum dt from list
//    dt = std::min(dt_diff,dt_adv);
//
//    // start by looping over x-variable (space)
//
//    // save previous column in holder
//    std::vector<double> previous_colU (runtime_args_->get_y_iterations(), 0.0);
//    std::vector<double> previous_colV (runtime_args_->get_y_iterations(), 0.0);
//    for (unsigned int i = 0 ; i < runtime_args_->get_y_iterations(); ++i){
//        previous_colU[i] = multiUSolutions_[0][i];
//        previous_colV[i] = multiVSolutions_[0][i];
//    }
//
//    for (unsigned int  i = 1 ; i < runtime_args_->get_x_iterations()-1; ++i){
//        double uSolnym1 = multiUSolutions_[i][0];
//        double vSolynm1 = multiVSolutions_[i][0];
//        for (unsigned int j = 1 ; j < runtime_args_->get_y_iterations()-1; ++j){
//            // Solve U, V into temporary holders first
//            double temp_u_at_ij;
//            double temp_v_at_ij;
//
//            // For easier reference
//            double uij = multiUSolutions_[i][j];
//            double uip1j = multiUSolutions_[i+1][j];
//            double uim1j = previous_colU[j];
//            double uijp1 = multiUSolutions_[i][j+1];
//            double uijm1 = uSolnym1;
//            double vij = multiVSolutions_[i][j];
//            double vip1j = multiVSolutions_[i+1][j];
//            double vim1j = previous_colV[j];
//            double vijp1 = multiVSolutions_[i][j+1];
//            double vijm1 = vSolynm1;
//            double fx = 0.5*dt/dx;
//            double fy = 0.5*dt/dy;
//            double fxdiff = nu*dt/(dx*dx);
//            double fydiff = nu*dt/(dy*dy);
//
//
//            temp_u_at_ij = uij - uij * fx * (uip1j - uim1j) + fabs(uij) * fx * (uip1j - 2.*uij + uim1j)
//                           - vij * fy * (uijp1 - uijm1) + fabs(vij) * fy * (uijp1 - 2.*uij + uijm1); // From burger
//
//            temp_u_at_ij += fxdiff * (uip1j - 2. * uij + uim1j) + fydiff * (uijp1 - 2. * uij + uijm1); // From diffusion
//
//            temp_v_at_ij = vij - uij * fx * (vip1j - vim1j) + fabs(uij) * fx * (vip1j - 2.*vij + vim1j)
//                           - vij * fy * (vijp1 - vijm1) + fabs(vij) * fy * (vijp1 - 2.*vij + vijm1); // From burger
//
//            temp_v_at_ij += fxdiff * (vip1j - 2. * vij + vim1j) + fydiff * (vijp1 - 2. * vij + vijm1); // From diffusion
//
//            // Merely add on separate portions from before.
//
//            // Write U, V into new form
//            multiUSolutions_[i][j] = temp_u_at_ij;
//            multiVSolutions_[i][j] = temp_v_at_ij;
//
//            uSolnym1 = multiUSolutions_[i][j];
//            vSolynm1 = multiVSolutions_[i][j];
//
//        }
//
//        //rewrite the previous_col
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            previous_colU[j] = multiUSolutions_[i][j];
//            previous_colV[j] = multiVSolutions_[i][j];
//        }
//
//    }
//
//    // Add to extra time
//    currentTime_ += dt;
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimBurger::set_boundary(){
//
//    multiSolutionBoundaryCondition_ = new BoxBoundaryCondition;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimBurger::set_init_cond(){
//
//    // hard code as wave type
//    multiSolutionInitialCondition_ = new Curvilinear;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}
//
////============================================================================================================
///*
// * Solution procedure for steady state laplace equation
// */
////============================================================================================================
//
////============================================================================================================
//void SolutionProcedureLaplace::start_procedure(std::string& runtime_params, std::string& template_file_name ){
//
//    // Since type is already known, just hard code in the file reader
//    runtime_args_ = new RuntimeParamMultiDim;
//
//    runtime_args_->read_parameters_from_file(runtime_params);
//
//    // Will take from input user to set boundary/initial condition type
//    // Fow now, hard code it
//
//    set_boundary();
//    set_init_cond();
//
//    // set size of vector
//    unsigned int XSIZE = runtime_args_->get_x_iterations();
//    unsigned int YSIZE = runtime_args_->get_y_iterations();
//    multiUSolutions_.resize(XSIZE, std::vector<double>(YSIZE));
//
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//
//    // For first step, need to apply the initial condition
//    double dx = runtime_args_->get_dx();
//    double xo = runtime_args_->get_xo();
//    multiSolutionInitialCondition_->apply_initial_cond(multiUSolutions_, *(runtime_args_));
//
//    // execute procedure
//    procedure(template_file_name);
//
//    // end procedure
//    end_procedure();
//
//}
//
////============================================================================================================
//void SolutionProcedureLaplace::convert_idx_to_pos_x(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dx = runtime_args_->get_dx();
//    double xo = runtime_args_->get_xo();
//
//    pos = (double)(idx) * dx + xo;
//
//}
////============================================================================================================
//void SolutionProcedureLaplace::convert_idx_to_pos_y(unsigned int idx, double& pos){
//
//    // Convert from index to position
//    // idx * dx + xo = current position
//    double dy = runtime_args_->get_dy();
//    double yo = runtime_args_->get_yo();
//
//    pos = (double)(idx) * dy + yo;
//
//}
////============================================================================================================
//void SolutionProcedureLaplace::write_to_file()
//{
//
//}
//
////============================================================================================================
//void SolutionProcedureMultiDimDiffusion::procedure(std::string& template_file_name){
//
//    // write first point -- initial condition onto file
//    double to = runtime_args_->get_to();
//    std::ofstream outFile (template_file_name + "0");
//    std::cout.fill(' ');
//
//    // Needed nested for loops
//    for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            //write to file (x,y) pair, with u(x,y) at current time step
//            double posx;
//            convert_idx_to_pos_x(i, posx);
//            double posy;
//            convert_idx_to_pos_y(j, posy);
//            outFile << std::setw(4) << posx << "     ";
//            outFile << posy << "     ";
//            outFile << multiUSolutions_[i][j] << std::endl;
//        }
//    }
//
//    outFile.close();
//
//
//    // Apply step at every time step
//    for(unsigned int t = 1; t < runtime_args_->get_time_iterations(); ++t){
//        apply_step();
//
//        // write results to output file, need x's, y's, and t's
//        outFile.open(template_file_name + std::to_string(t));
//
//        // I hate the idea of a triply nested for loop, but this just prints out data
//        for (unsigned int i = 0 ; i < runtime_args_->get_x_iterations(); ++i){
//            for (unsigned int  j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//                //write to file (x,y) pair, with u(x,y) at current time step
//                double posx;
//                convert_idx_to_pos_x(i,posx);
//                double posy;
//                convert_idx_to_pos_y(j,posy);
//                outFile << std::setw(4) << posx << "     ";
//                outFile << posy << "     ";
//                outFile << multiUSolutions_[i][j] << std::endl;
//            }
//        }
//
//        outFile.close();
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureLaplace::end_procedure(){
//    // something important to finish writing the output file
//    std::cout << "Done." << std::endl;
//}
//
////============================================================================================================
//void SolutionProcedureLaplace::apply_step(){
//
//    // Apply individual step to solutions
//
//    multiSolutionBoundaryCondition_->enforce_boundary_conditions(multiUSolutions_, *(runtime_args_));
//
//    // loops over entries, save the one involved in the boundary condition
//    // boundary hard coded to LHS of wall
//
//
//    double c = runtime_args_->get_c(); // In this case, the nu coefficient in the equation
//    double dx = runtime_args_->get_dx();
//    double dy = runtime_args_->get_dy();
//    double dt = runtime_args_->get_dt();
//
//    // Assume the user has specified a reasonable set of conditions
//    // Else, throw a warning to the user that this has occured.
//
//    // start by looping over x-variable (space)
//
//    // save previous column in holder
//    std::vector<double> previous_col (runtime_args_->get_y_iterations(), 0.0);
//    for (unsigned int i = 0 ; i < runtime_args_->get_y_iterations(); ++i){
//        previous_col[i] = multiUSolutions_[0][i];
//    }
//
//    for (unsigned int  i = 1 ; i < runtime_args_->get_x_iterations()-1; ++i){
//        double uSolnym1 = multiUSolutions_[i][0];
//        for (unsigned int j = 1 ; j < runtime_args_->get_y_iterations()-1; ++j){
//            double uij = multiUSolutions_[i][j];
//            double uip1j = multiUSolutions_[i+1][j];
//            double uim1j = previous_col[j];
//            double uijp1 = multiUSolutions_[i][j+1];
//            double uijm1 = uSolnym1;
//            double fx = c*dt/(dx*dx);
//            double fy = c*dt/(dy*dy);
//
//            multiUSolutions_[i][j] = uij + fx * (uip1j - 2. * uij + uim1j)
//                                     + fy * (uijp1 - 2. * uij + uijm1);
//
//
//            uSolnym1 = multiUSolutions_[i][j];
//        }
//
//        //rewrite the previous_col
//        for (unsigned int j = 0 ; j < runtime_args_->get_y_iterations(); ++j){
//            previous_col[j] = multiUSolutions_[i][j];
//        }
//
//    }
//
//}
//
////============================================================================================================
//void SolutionProcedureLaplace::set_boundary(){
//
//    multiSolutionBoundaryCondition_ = new BoxBoundaryCondition;
//
//
//    //TODO: should probably allow the user to go ahead and set what boundary condition to be used
//
//}
//
////============================================================================================================
//void SolutionProcedureLaplace::set_init_cond(){
//
//    // hard code as wave type
//    multiSolutionInitialCondition_ = new Curvilinear;
//
//    // TODO: should probably allow the user to go ahead and set the initial condition from the control file
//
//}