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
#include "FluidEquation.h"


/*
 * Procedure
 */
class SolutionProcedure{
public:
    void end_procedure();
    void start_procedure(std::string& runtime_params, std::string& template_file_name, std::string& whichPDE);
    void procedure(std::string& template_file_name);
protected:
    InitCond* solutionInitialCondition_;
};
#endif //CFD_HW_SOLUTIONPROCEDURE_H
