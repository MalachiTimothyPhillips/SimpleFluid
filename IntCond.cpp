//
// Created by Malachi Phillips on 9/25/16.
//

#include "IntCond.h"
#include <math.h>
////Factory
//InitCond* InitCond::make_initial_condition()
//{
//    return new StepWave; //Hard code in single initial condition, for now
//}

// Essentially, what the user has to specify is the function
void StepWave::wave_func(double pos, double& value){
    /*
     * u(x,o) = f(x)
     * Do:: from x = 1 to x = 3, u(x,o) = 3. Else, u(x,o) = 1
     *
     */

    if (pos <= 3.0 && pos >= 1.0){
        value = 3.0;
    }
    else{
        value = 1.0;
    }

    //TODO:: avoid conflicts with Boundary Condition

}

void StepWave::convert_idx_to_position(int idx, double dx, double xo, double& pos)
{
    // Convert from index to position
    // idx * dx + xo = current position
    pos = (double)(idx) * dx + xo;
}

void StepWave::apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo){
    // u_solutions is already sized, so find number of points to apply from there
//    double dx = runtime.get_dx();
//    double xo = runtime.get_xo();

    // Runtime parameters already set in SolutionProcedure
    // Ensure no issues with boundary conditions, leave off first point
    for (unsigned int i = 1 ; i < u_solutions.size(); ++i){
        double currentPos;
        convert_idx_to_position(i, dx, xo, currentPos);
        // use current position to get function value
        double waveVal;
        wave_func(currentPos, waveVal);
        u_solutions[i] = waveVal;
    }

}

void SinWave::apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo){
    // u_solutions is already sized, so find number of points to apply from there
//    double dx = runtime.get_dx();
//    double xo = runtime.get_xo();

    // Runtime parameters already set in SolutionProcedure
    // Ensure no issues with boundary conditions, leave off first point
    for (unsigned int i = 1 ; i < u_solutions.size(); ++i){
        double currentPos;
        convert_idx_to_position(i, dx, xo, currentPos);
        // use current position to get function value
        double waveVal;
        sin_func(currentPos, waveVal);
        u_solutions[i] = waveVal;
    }

    // TODO: basic form of this is shared. Perhaps template this class?
}

void SinWave::convert_idx_to_position(int idx, double dx, double xo, double& pos)
{
    // Convert from index to position
    // idx * dx + xo = current position
    pos = (double)(idx) * dx + xo;

    // TODO: shared between both, should probably tuck away in the base class for now
}

void SinWave::sin_func(double pos, double &value) {

    /*
     * u(x,0) = f(x)
     * f(x) = 10sin(0.1x) //radian mode
     *
     */

    value = sin(0.1 * pos) * 10.0;

}

void PositiveWave::convert_idx_to_position(int idx, double dx, double xo, double& pos)
{
    // Convert from index to position
    // idx * dx + xo = current position
    pos = (double)(idx) * dx + xo;

    // TODO: shared between both, should probably tuck away in the base class for now
}

void PositiveWave::pos_func(double pos, double& wave_func){
    wave_func = 2. - sin(0.1 * pos);
}

void PositiveWave::apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo){
    // u_solutions is already sized, so find number of points to apply from there
//    double dx = runtime.get_dx();
//    double xo = runtime.get_xo();

    // Runtime parameters already set in SolutionProcedure
    // Ensure no issues with boundary conditions, leave off first point
    for (unsigned int i = 1 ; i < u_solutions.size(); ++i){
        double currentPos;
        convert_idx_to_position(i, dx, xo, currentPos);
        // use current position to get function value
        double waveVal;
        pos_func(currentPos, waveVal);
        u_solutions[i] = waveVal;
    }

    // TODO: basic form of this is shared. Perhaps template this class?
}

void Curvilinear::convert_idx_to_position(int idx, double dx, double xo, double& pos)
{
    // Convert from index to position
    // idx * dx + xo = current position
    pos = (double)(idx) * dx + xo;

    // TODO: shared between both, should probably tuck away in the base class for now
}

void Curvilinear::pos_func(double x, double y, double& u){
    u = x * x * cos(y) + sin(y);
}

void Curvilinear::apply_initial_cond(std::vector<std::vector<double>>& uSoln, RuntimeParamMultiDim& runtime){
    double dx = runtime.get_dx();
    double dy = runtime.get_dy();
    double xo = runtime.get_xo();
    double yo = runtime.get_yo();
    unsigned int sizeX = runtime.get_x_iterations();
    unsigned int sizeY = runtime.get_y_iterations();
    for (unsigned int i = 0 ; i < sizeX; ++i){
        for (unsigned int j = 0 ; j < sizeY; ++j){
            double x_curr;
            convert_idx_to_position(i, dx, xo, x_curr);
            double y_curr;
            convert_idx_to_position(j, dy, yo, y_curr);

            double waveVal;
            pos_func(x_curr, y_curr, waveVal);
            uSoln[i][j] = waveVal;

        }
    }

    // TODO: basic form of this is shared. Perhaps template this class?
}