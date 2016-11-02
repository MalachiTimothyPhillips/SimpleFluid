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

// Base class for initial condition
class InitCond{
public:
    // factory to make initial conditions
    static InitCond *make_initial_condition();
    virtual void apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo) = 0;
protected:
private:
};

// Base class for multidimensional initial condition
class MutliDimInitCond{
public:
    // factory to make initial conditions
    static InitCond *make_initial_condition();
    virtual void apply_initial_cond(std::vector<std::vector<double>>& uSol, RuntimeParamMultiDim& runtime) = 0;
protected:
private:
};


/*
 * Sin wave
 */

class SinWave: public InitCond{
public:
    void apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo);
    void convert_idx_to_position(int idx, double dx, double xo, double& pos);
    void sin_func(double pos, double& value);
protected:
private:

};

/*
 * Square wave profile
 */

class StepWave : public InitCond{
public:
    void apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo);
    void convert_idx_to_position(int idx, double dx, double xo, double& pos);
    void wave_func(double pos, double& value);

protected:


private:
    //virtual ~StepWave(){};
};

/*
 * Positive Wave
 */

class PositiveWave : public InitCond{
public:
    void apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo);
    void convert_idx_to_position(int idx, double dx, double xo, double& pos);
    void pos_func(double pos, double& value);

protected:


private:
    //virtual ~StepWave(){};
};

class Curvilinear : public MutliDimInitCond{
public:
    void apply_initial_cond(std::vector<std::vector<double>>& uSol, RuntimeParamMultiDim& runtime);
    void convert_idx_to_position(int idx, double dx, double xo, double& pos);
    void pos_func(double x, double y, double& u);
protected:
private:
};
#endif //CFD_HW_INTCOND_H
