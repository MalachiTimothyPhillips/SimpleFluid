//
// Created by Malachi Phillips on 9/25/16.
//

/*
 * Initial condition u(x,o) = f(x) for integration
 */

#ifndef CFD_HW_INTCOND_H
#define CFD_HW_INTCOND_H

#include <vector>
#include "RuntimeParameters.h"
#include "FluidEquation.h"


#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iostream>
#include <iomanip>

/*
 * Forward declaration
 */

class FluidEquation;

typedef std::vector<std::vector<double>> matrix;

// Base class for initial condition
class InitCond{
public:
    InitCond *make_initial_condition(std::string& init_cond);
    FluidEquation* make_fluid_equation(std::string& equationType);
    virtual void apply_initial_cond() = 0;
    void convert_idx_to_pos(unsigned int idx, double& pos);
    virtual void enforce_boundary() = 0;

    std::vector<double>& return_args(){return args_;};

    void set_args();

    FluidEquation* fluidEquation_;

protected:
    std::vector<double> args_;// = std::vector<double>(11); //pre-determined size
private:
};

// Sin Wave
class SinWave : public InitCond{
public:
    void apply_initial_cond();
    void enforce_boundary();
protected:
private:

};


//
///*
// * Sin wave
// */
//
//class SinWave: public InitCond{
//public:
//    void apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo);
//    void convert_idx_to_position(int idx, double dx, double xo, double& pos);
//    void sin_func(double pos, double& value);
//protected:
//private:
//
//};
//
///*
// * Square wave profile
// */
//
//class StepWave : public InitCond{
//public:
//    void apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo);
//    void convert_idx_to_position(int idx, double dx, double xo, double& pos);
//    void wave_func(double pos, double& value);
//
//protected:
//
//
//private:
//    //virtual ~StepWave(){};
//};
//
///*
// * Positive Wave
// */
//
//class PositiveWave : public InitCond{
//public:
//    void apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo);
//    void convert_idx_to_position(int idx, double dx, double xo, double& pos);
//    void pos_func(double pos, double& value);
//
//protected:
//
//
//private:
//    //virtual ~StepWave(){};
//};
//
//class Curvilinear : public MutliDimInitCond{
//public:
//    void apply_initial_cond(std::vector<std::vector<double>>& uSol, RuntimeParamMultiDim& runtime);
//    void convert_idx_to_position(int idx, double dx, double xo, double& pos);
//    void pos_func(double x, double y, double& u);
//protected:
//private:
//};
#endif //CFD_HW_INTCOND_H
