//
// Created by Malachi Phillips on 9/25/16.
//

#ifndef CFD_HW_RUNTIMEPARAMETERS_H
#define CFD_HW_RUNTIMEPARAMETERS_H

// Controls solution parameters, from user
/*
 * User must input grid spacing (dx), length, and total time.
 * CFL, dt are inferred from these settings to ensure a reasonable PDE solution
 */

#include <string>
#include <vector>

//============================================================================================================
/*
 * RuntimeParameter: parse file from input deck
 */
//============================================================================================================
class RuntimeParser{
public:
    // Read in parameters from input file
    void read_parameters_from_file(std::string& runtime_args);
    std::vector<double>& return_args(){return args_;};
    std::string& return_cond(){return initialCondition_;};
protected:
    std::vector<double> args_;
    std::string initialCondition_;
private:
};

#endif //CFD_HW_RUNTIMEPARAMETERS_H