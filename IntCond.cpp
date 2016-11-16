//
// Created by Malachi Phillips on 9/25/16.
//

#include "IntCond.h"
#include "FluidEquation.h"
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iostream>
#include <iomanip>

//Factory
InitCond* InitCond::make_initial_condition(std::string& init_cond)
{
    // Handle initial condition
    if(init_cond == "SinWave"){
        return new SinWave;
    }

}

//Factory
FluidEquation* InitCond::make_fluid_equation(std::string &equationType) {

    std::cout << "Reach here?" << std::endl;
    if(equationType == "UpwindLinWave"){
        std::cout << "Inside comparison" << std::endl;
        return new UpwindLinWave(args_);
    }

    std::cout << "Handled it correctly" << std::endl;
}

void InitCond::set_args() {
    args_.resize(11);
    std::cout << "Success" << std::endl;
//    for (unsigned int i = 0; i < args.size(); ++i) {
//        std::cout << args.size() << std::endl;
//        std::cout << args_.size() << std::endl;
//        std::cout << args[i] << std::endl;
//        args_.push_back((args[i]));// = args[i];
//        //args_[i] = args[i];
//        std::cout << args[i] << std::endl;
//        std::cout << "Did it work?" << std::endl;
//
//    }
}

void InitCond::convert_idx_to_pos(unsigned int idx, double &pos) {
    pos = (double)idx * fluidEquation_->get_dx() + fluidEquation_->get_lo();
}

void SinWave::apply_initial_cond() {
    /*
     * Sin wave, 10sin(0.1x)
     */
    for (unsigned int i = 0; i < fluidEquation_->uSolutions_.size(); ++i){
        double pos;
        convert_idx_to_pos(i,pos);
        fluidEquation_->uSolutions_[i] = 10.0 * sin(0.1*pos);
    }
}

void SinWave::enforce_boundary() {
    /*
     * At end points, evaluate initial condition function
     */
    fluidEquation_->uSolutions_[0] = 10.0 * sin(0.1 * fluidEquation_->get_lo());
    fluidEquation_->uSolutions_[fluidEquation_->get_nl()] = 10.0 * sin(0.1 * fluidEquation_->get_lf());
}

//// Essentially, what the user has to specify is the function
//void StepWave::wave_func(double pos, double& value){
//    /*
//     * u(x,o) = f(x)
//     * Do:: from x = 1 to x = 3, u(x,o) = 3. Else, u(x,o) = 1
//     *
//     */
//
//    if (pos <= 3.0 && pos >= 1.0){
//        value = 3.0;
//    }
//    else{
//        value = 1.0;
//    }
//
//    //TODO:: avoid conflicts with Boundary Condition
//
//}
//
//void StepWave::convert_idx_to_position(int idx, double dx, double xo, double& pos)
//{
//    // Convert from index to position
//    // idx * dx + xo = current position
//    pos = (double)(idx) * dx + xo;
//}
//
//void StepWave::apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo){
//    // u_solutions is already sized, so find number of points to apply from there
////    double dx = runtime.get_dx();
////    double xo = runtime.get_xo();
//
//    // Runtime parameters already set in SolutionProcedure
//    // Ensure no issues with boundary conditions, leave off first point
//    for (unsigned int i = 1 ; i < u_solutions.size(); ++i){
//        double currentPos;
//        convert_idx_to_position(i, dx, xo, currentPos);
//        // use current position to get function value
//        double waveVal;
//        wave_func(currentPos, waveVal);
//        u_solutions[i] = waveVal;
//    }
//
//}
//
//void SinWave::apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo){
//    // u_solutions is already sized, so find number of points to apply from there
////    double dx = runtime.get_dx();
////    double xo = runtime.get_xo();
//
//    // Runtime parameters already set in SolutionProcedure
//    // Ensure no issues with boundary conditions, leave off first point
//    for (unsigned int i = 1 ; i < u_solutions.size(); ++i){
//        double currentPos;
//        convert_idx_to_position(i, dx, xo, currentPos);
//        // use current position to get function value
//        double waveVal;
//        sin_func(currentPos, waveVal);
//        u_solutions[i] = waveVal;
//    }
//
//    // TODO: basic form of this is shared. Perhaps template this class?
//}
//
//void SinWave::convert_idx_to_position(int idx, double dx, double xo, double& pos)
//{
//    // Convert from index to position
//    // idx * dx + xo = current position
//    pos = (double)(idx) * dx + xo;
//
//    // TODO: shared between both, should probably tuck away in the base class for now
//}
//
//void SinWave::sin_func(double pos, double &value) {
//
//    /*
//     * u(x,0) = f(x)
//     * f(x) = 10sin(0.1x) //radian mode
//     *
//     */
//
//    value = sin(0.1 * pos) * 10.0;
//
//}
//
//void PositiveWave::convert_idx_to_position(int idx, double dx, double xo, double& pos)
//{
//    // Convert from index to position
//    // idx * dx + xo = current position
//    pos = (double)(idx) * dx + xo;
//
//    // TODO: shared between both, should probably tuck away in the base class for now
//}
//
//void PositiveWave::pos_func(double pos, double& wave_func){
//    wave_func = 2. - sin(0.1 * pos);
//}
//
//void PositiveWave::apply_initial_cond(std::vector<double>& u_solutions, double dx, double xo){
//    // u_solutions is already sized, so find number of points to apply from there
////    double dx = runtime.get_dx();
////    double xo = runtime.get_xo();
//
//    // Runtime parameters already set in SolutionProcedure
//    // Ensure no issues with boundary conditions, leave off first point
//    for (unsigned int i = 1 ; i < u_solutions.size(); ++i){
//        double currentPos;
//        convert_idx_to_position(i, dx, xo, currentPos);
//        // use current position to get function value
//        double waveVal;
//        pos_func(currentPos, waveVal);
//        u_solutions[i] = waveVal;
//    }
//
//    // TODO: basic form of this is shared. Perhaps template this class?
//}
//
//void Curvilinear::convert_idx_to_position(int idx, double dx, double xo, double& pos)
//{
//    // Convert from index to position
//    // idx * dx + xo = current position
//    pos = (double)(idx) * dx + xo;
//
//    // TODO: shared between both, should probably tuck away in the base class for now
//}
//
//void Curvilinear::pos_func(double x, double y, double& u) {
////    double x_shift = x-5.;
////    double y_shift = y-5.;
////
////    double dist = sqrt(x_shift * x_shift + y_shift * y_shift);
////
////    u = fabs(sin(dist)/(dist));
////
////    // if radius is greater than 2.85234173635078, set to 0.1
////    if (dist >=  2.85234173635078){
////        u = 0.1;
////    }
//
//    //u = sin(x)*sin(y);
//
//    // Square wave
//    // draw square wave size (0,0),(1,0),(0,1)(1,1)
//    if ( x < 1.0 && y < 1.0){
//        u = 5.0;
//    }
//    else{
//        u = 1.0;
//    }
//
//}
//
//void Curvilinear::apply_initial_cond(std::vector<std::vector<double>>& uSoln, RuntimeParamMultiDim& runtime){
//    double dx = runtime.get_dx();
//    double dy = runtime.get_dy();
//    double xo = runtime.get_xo();
//    double yo = runtime.get_yo();
//    unsigned int sizeX = runtime.get_x_iterations();
//    unsigned int sizeY = runtime.get_y_iterations();
//    for (unsigned int i = 0 ; i < sizeX; ++i){
//        for (unsigned int j = 0 ; j < sizeY; ++j){
//            double x_curr;
//            convert_idx_to_position(i, dx, xo, x_curr);
//            double y_curr;
//            convert_idx_to_position(j, dy, yo, y_curr);
//
//            double waveVal;
//            pos_func(x_curr, y_curr, waveVal);
//            uSoln[i][j] = waveVal;
//
//        }
//    }
//
//    // TODO: basic form of this is shared. Perhaps template this class?
//}