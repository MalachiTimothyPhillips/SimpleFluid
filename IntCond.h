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

    InitCond(std::vector<double>& args);

    InitCond *make_initial_condition(std::string& init_cond, std::vector<double>& args);
    virtual void apply_initial_cond() = 0;
    void convert_idx_to_pos(unsigned int idx, double& pos);
    virtual void enforce_boundary() = 0;

    FluidEquation* fluidEquation_;

protected:

    // Initial condition (constructed via std::vector<double>& args) should always have these base parameters
    double lo_;
    unsigned int nl_;
    double dx_;
    double lf_;

private:
};

// Sin Wave
class SinWave : public InitCond{
public:
    SinWave(std::vector<double>& args);
    void apply_initial_cond();
    void enforce_boundary();
protected:
private:

};

// Sod Shock Tube
class SodShockTube : public InitCond{
public:
    SodShockTube(std::vector<double>& args);
    void apply_initial_cond();
    void enforce_boundary();
protected:

    double gamma_ = 1.4; // ideal gas

    // Because magic numbers are bad
    double P_l_ = 1.0;
    double P_r_ = 0.1;

    double rho_l_ = 1.0;
    double rho_r_ = 0.125;

    double u_l_ = 0.0;
    double u_r_ = 0.0;

private:
};

// 2D stepwave at 0
class Curvilinear : public InitCond{
public:
    Curvilinear(std::vector<double>& args);
    void apply_initial_cond();
    void enforce_boundary();
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
    void pos_func(double x, double y, double& u);
protected:

    // Additional initial conditions from 2D
    double ho_;
    unsigned int nh_;
    double dy_;
    double hf_;

private:
};

// 2D stepwave at 0, but applied for u and v
class ExtendedCurvilinear : public InitCond{
public:
    ExtendedCurvilinear(std::vector<double>& args);
    void apply_initial_cond();
    void enforce_boundary();
    void pos_func(double x, double y, double& u);
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
protected:
    // Additional initial conditions from 2D
    double ho_;
    unsigned int nh_;
    double dy_;
    double hf_;
private:

};

// 2D stepwave at 0, but applied for u and v
class LaplacianBoundary : public InitCond{
public:
    LaplacianBoundary(std::vector<double>& args);
    void apply_initial_cond();
    void enforce_boundary();
    //void pos_func(double x, double y, double& u);
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
protected:
    // Additional initial conditions from 2D
    double ho_;
    unsigned int nh_;
    double dy_;
    double hf_;
private:
};

// Do absolutely nothing
class DoNothing : public InitCond{
public:
    DoNothing(std::vector<double>& args);
    void apply_initial_cond();
    void enforce_boundary();
    //void pos_func(double x, double y, double& u);
    void convert_idx_to_pos_y(unsigned int idx, double& pos);
protected:
    // Additional initial conditions from 2D
    double ho_;
    unsigned int nh_;
    double dy_;
    double hf_;
private:
};

#endif //CFD_HW_INTCOND_H
