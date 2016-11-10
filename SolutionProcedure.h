//
// Created by Malachi Phillips on 9/30/16.
//

#ifndef CFD_HW_SOLUTIONPROCEDURE_H
#define CFD_HW_SOLUTIONPROCEDURE_H

#include <vector>
#include "RuntimeParameters.h"
#include "BoundaryCond.h"
#include "IntCond.h"
#include <boost/multi_array.hpp>
#include <cassert>


/*
 * Base class procedure
 */

class SolutionProcedure{
public:

    static SolutionProcedure* determine_solution_procedure(std::string& equationType);
    virtual void apply_step() = 0;
    virtual void end_procedure() = 0;
    virtual void start_procedure(std::string& runtime_params, std::string& template_file_name) = 0;
    virtual void procedure(std::string& template_file_name) = 0;
protected:
    BoundaryCondition* solutionBoundaryCondition_;
    InitCond* solutionInitialCondition_;

    std::vector<double> uSolutions_;

    virtual void set_boundary() = 0;
    virtual void set_init_cond() = 0;

};

/*
 * Derived class for 2D cases
 */
class MultiDimProcedure : public SolutionProcedure {
protected:
    std::vector<std::vector<double>> multiUSolutions_;
    MultiDimensionBoundaryCondition* multiSolutionBoundaryCondition_;
    MutliDimInitCond* multiSolutionInitialCondition_;
};


/*
 * Handle solution procedure for solving ut + cux = 0
 */

class SolutionProcedureLinWave : public SolutionProcedure{
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos(unsigned int idx, double& pos);
    void write_to_file();

protected:

    RuntimeParametersLinWave* runtime_param_;

    void set_boundary();
    void set_init_cond();

private:
    //~virtual SolutionProcedure(){};
};

/*
 * Handle solution procedure for solving ut + cux = 0
 * Handles case when c may be positive or negative
 */
class SolutionProcedureFixLinWave : public SolutionProcedure{
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos(unsigned int idx, double& pos);
    void write_to_file();

protected:

    RuntimeParametersLinWave* runtime_param_; // same params as above

    void set_boundary();
    void set_init_cond();

private:
    //~virtual SolutionProcedure(){};
};


/*
 * Handle solution procedure for solving ut = vuxx
 */

class SolutionProcedureDiffusion : public SolutionProcedureLinWave{
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos(unsigned int idx, double& pos);
    void write_to_file();
protected:

    RuntimeParametersDiffusion* runtime_param_;

    void set_boundary();
    void set_init_cond();
};

/*
 * Handle solution procedure for solving ut + uux = 0
 */

class SolutionProcedureInviscid : public SolutionProcedureLinWave  {
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos(unsigned int idx, double& pos);
    void write_to_file();
protected:

    RuntimeParametersLinWave* runtime_param_;

    void set_boundary();
    void set_init_cond();
};

/*
 * Handle solution procedure for solving ut + uux = vuxx
 */

class SolutionProcedureBurger : public SolutionProcedureLinWave  {
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos(unsigned int idx, double& pos);
    void write_to_file();
protected:

    RuntimeParametersBurger* runtime_param_;

    void set_boundary();
    void set_init_cond();
};

/*
 * Solution procedure for 2D advection equation
 */
class SolutionProcedureMultiDimAdvection : public MultiDimProcedure  {
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos_x(unsigned int idx, double& pos);
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
    void write_to_file();
protected:

    // Multi dimensional parameters
    RuntimeParamMultiDim* runtime_args_;

    void set_boundary();
    void set_init_cond();

};


/*
 * Solution procedure for 2D advection equation
 */
class SolutionProcedureMultiDimNonLinAdvEqn : public MultiDimProcedure  {
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos_x(unsigned int idx, double& pos);
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
    void write_to_file();
    double get_current_time() {return currentTime_; } ;
protected:

    // Multi dimensional parameters
    RuntimeParamMultiDim* runtime_args_;

    // For my other sub problem
    std::vector<std::vector<double>> multiVSolutions_;

    double currentTime_ = 0;

    void set_boundary();
    void set_init_cond();

};

/*
 * Solution procedure for 2D diffusion equation
 */
class SolutionProcedureMultiDimDiffusion : public MultiDimProcedure  {
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos_x(unsigned int idx, double& pos);
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
    void write_to_file();
protected:

    // Multi dimensional parameters
    RuntimeParamMultiDim* runtime_args_;

    void set_boundary();
    void set_init_cond();

};

/*
 * Solution procedure for 2D viscous burger equation
 */
class SolutionProcedureMultiDimBurger : public MultiDimProcedure  {
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos_x(unsigned int idx, double& pos);
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
    void write_to_file();
    double get_current_time() {return currentTime_; } ;
protected:

    // Multi dimensional parameters
    RuntimeParamMultiDim* runtime_args_;

    // For my other sub problem
    std::vector<std::vector<double>> multiVSolutions_;

    double currentTime_ = 0;

    void set_boundary();
    void set_init_cond();

};

/*
 * Solution procedure for 2D Laplace equation
 */
class SolutionProcedureLaplace : public MultiDimProcedure  {
public:
    void apply_step();
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name);
    void procedure(std::string& template_file_name);
    void convert_idx_to_pos_x(unsigned int idx, double& pos);
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
    void write_to_file();
protected:

    // Multi dimensional parameters
    RuntimeParamMultiDim* runtime_args_;

    void set_boundary();
    void set_init_cond();

};




#endif //CFD_HW_SOLUTIONPROCEDURE_H
