//
// Created by Malachi Phillips on 9/30/16.
//

#ifndef CFD_HW_SOLUTIONPROCEDURE_H
#define CFD_HW_SOLUTIONPROCEDURE_H

#include <vector>
#include "RuntimeParameters.h"
#include "BoundaryCond.h"
#include "IntCond.h"

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
    std::vector<double> uSolutions_; // solutions vector, always exists
    BoundaryCondition* solutionBoundaryCondition_;
    InitCond* solutionInitialCondition_;

    virtual void set_boundary() = 0;
    virtual void set_init_cond() = 0;

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

#endif //CFD_HW_SOLUTIONPROCEDURE_H
