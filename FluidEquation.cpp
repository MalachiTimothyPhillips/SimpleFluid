//
// Created by Malachi Phillips on 11/10/16.
//

#include "FluidEquation.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iostream>
#include <iomanip>
//============================================================================================================
/*
 * Fluid Equation base class
 */
//============================================================================================================
FluidEquation::FluidEquation(std::vector<double> &args) {
    std::cout << "Does it ever reach here?" << std::endl;
    std::cout << DOF_IDS::lo << std::endl;
    lo_ = args[DOF_IDS::lo];
    std::cout << "Work?" << std::endl;

    lf_ = args[DOF_IDS::lf];
    std::cout << "Check?" << std::endl;

    nl_ = (unsigned int) args[DOF_IDS::nl];
    nt_ = (unsigned int) args[DOF_IDS::nt];

    // Calculate dx
    dx_ = (lf_ - lo_) / (double) nl_;

} // constructor
//-----------------------------------------------------------------------------------------------------------
void FluidEquation::write_to_file(std::string& template_file_name, unsigned int currentStep){
 //Write current time step to the file
    std::ofstream outFile(template_file_name + std::to_string(currentStep));
    for (unsigned int i = 0; i < uSolutions_.size(); ++i) {
        double pos;
        convert_idx_to_pos(i, pos);
        outFile << pos << " ";
        outFile << uSolutions_[i] << std::endl;
    } // writes file in two column format: x and u(x)
    outFile.close();
}

//-----------------------------------------------------------------------------------------------------------
FluidEquation* FluidEquation::make_fluid_equation(std::string& equationType, std::vector<double>& args){
    std::cout << "Reach here?" << std::endl;
    if(equationType == "UpwindLinWave"){
        std::cout << "Inside comparison" << std::endl;
        return new UpwindLinWave(args);
    }

    std::cout << "Handled it correctly" << std::endl;
}

//-----------------------------------------------------------------------------------------------------------
void FluidEquation::convert_idx_to_pos(unsigned int idx, double &pos) {
    pos = (double) idx * dx_ + lo_;
}

////============================================================================================================
///*
// * One dimensional fluid equation base class
// */
////============================================================================================================
//
//OneDimFluidEquation::OneDimFluidEquation(std::vector<double>& args) : FluidEquation(args){
//        std::cout << "How about here?" << std::endl;
//        // Resize uSolutions_ to appropriate size
//        uSolutions_.assign(nl_, 0.0);
//} //  constructor
//
////-----------------------------------------------------------------------------------------------------------
//void OneDimFluidEquation::write_to_file(std::string &template_file_name, unsigned int currentStep) {
//    // Write current time step to the file
//    std::ofstream outFile(template_file_name + std::to_string(currentStep));
//    for (unsigned int i = 0; i < uSolutions_.size(); ++i) {
//        double pos;
//        convert_idx_to_pos(i, pos);
//        outFile << pos << " ";
//        outFile << uSolutions_[i] << std::endl;
//    } // writes file in two column format: x and u(x)
//    outFile.close();
//}
//
////-----------------------------------------------------------------------------------------------------------
//void OneDimFluidEquation::convert_idx_to_pos(unsigned int idx, double &pos) {
//    pos = (double) idx * dx_ + lo_;
//}
//
////============================================================================================================
///*
// * Linear Wave Equation Base Class
// */
////============================================================================================================
//
//LinearWaveEquation::LinearWaveEquation(std::vector<double>&args) : OneDimFluidEquation(args){
//    std::cout << "Does it ever reach this bit?" << std::endl;
//    c_ = args[DOF_IDS::c];
//    tf_ = args[DOF_IDS::tf];
//
//    // Calculate the dt, CFL
//    dt_ = tf_ / (double) nt_;
//    CFL_ = c_*dt_/dx_;
//}

//============================================================================================================
/*
 * Linear Wave Equation, solved using FTBS/FTFS upwind/downwind shifter
 */
//============================================================================================================

UpwindLinWave::UpwindLinWave(std::vector<double>& args) : FluidEquation(args){
    std::cout << "How about here?" << std::endl;
    // Resize uSolutions_ to appropriate size
    uSolutions_.assign(nl_, 0.0);
    c_ = args[DOF_IDS::c];
    tf_ = args[DOF_IDS::tf];

    // Calculate the dt, CFL
    dt_ = tf_ / (double) nt_;
    CFL_ = c_*dt_/dx_;
}

//------------------------------------------------------------------------------------------------------------
void UpwindLinWave::apply_step() {
    // Apply single step
    // Do not include boundary points in solution
    double uSolnm1 = uSolutions_[0]; // Store solution of u_(i-1) at time n

    for (unsigned int i = 1; i < nl_; i++) {
        // foward in time, backward in space

        // make sure to use u_(i-1) state at time n, not n+1, so have to save previous value

        uSolutions_[i] = uSolutions_[i] - CFL_ / 2. * (uSolutions_[i + 1] - uSolnm1)
                         + std::abs(CFL_ / 2.) * (uSolutions_[i + 1] - 2 * uSolutions_[i] + uSolnm1);

        // write new position, will be i-1 at next pass
        uSolnm1 = uSolutions_[i];
    }
}

//------------------------------------------------------------------------------------------------------------
void UpwindLinWave::write_to_file(std::string &template_file_name, unsigned int currentStep) {
    FluidEquation::write_to_file(template_file_name, currentStep);
}

////============================================================================================================
///*
// * One dimensional diffusion equation - base case
// */
////============================================================================================================
//
//// Constructor hidden in away in header file
//
////============================================================================================================
///*
// * Diffusion Equation, solved in FTCS
// */
////============================================================================================================
//
//void DiffusionEquationFTCS::apply_step() {
//    // need to store u_(i-1) at time n
//    double uSolnm1 = uSolutions_[0];
//
//    for (unsigned int i = 1; i < nl_ - 1; i++) { // ensure it won't go too far
//        // foward in time, central in space
//        uSolutions_[i] = (1.0 - 2.0 * alpha_) * uSolutions_[i] +
//                         alpha_ * (uSolutions_[i + 1] + uSolnm1);
//        // values at end points are handled by the specified boundary conditions
//
//        // update previous state
//        uSolnm1 = uSolutions_[i];
//    }
//}
//
////------------------------------------------------------------------------------------------------------------
//void DiffusionEquationFTCS::write_to_file(std::string &template_file_name, unsigned int currentStep) {
//    OneDimFluidEquation::write_to_file(template_file_name, currentStep);
//}
//
////============================================================================================================
///*
// * Burger Equation Base Class
// */
////============================================================================================================
//
//// Constructor hidden away in base class
//
////============================================================================================================
///*
// * Burger Equation, Solved FTBS/FTCS in upwind/downwindshift with adaptive time step
// */
////============================================================================================================
//void BurgerEquationFTCS::apply_step() {
//
//    // Record u(i-1) at time step n
//    double uSolnm1 = uSolutions_[0];
//
//    // Calculate the time step needed to ensure stability, using the numeric limit eps
//    // In order to do this, need the maximum speed at the current step
//    double Cmax = *(std::max_element(std::begin(uSolutions_), std::end(uSolutions_)));
//
//    dt_ = dx_ * eps_ / Cmax;
//
//
//    for (unsigned int i = 1; i < nl_; i++) {
//        // foward in time, backward in space
//        double CFL = uSolutions_[i] * dt_ / dx_;
//        uSolutions_[i] = uSolutions_[i] - CFL * (uSolutions_[i] - uSolnm1);
//
//        uSolnm1 = uSolutions_[i];
//    }
//
//    // Add dt_ to find current time
//    T_ += dt_;
//}
//
////------------------------------------------------------------------------------------------------------------
//void BurgerEquationFTCS::write_to_file(std::string &template_file_name, unsigned int currentStep) {
//    OneDimFluidEquation::write_to_file(template_file_name, currentStep);
//}
//
//////============================================================================================================
/////*
//// * Two dimensional fluid equation base class
//// */
//////============================================================================================================
////    TwoDimFluidEquation::TwoDimFluidEquation(std::vector<double> &args) : FluidEquation(args) {
////        ho_ = args[DOF_IDS::ho];
////        hf_ = args[DOF_IDS::hf];
////        nh_ = (unsigned int) args[DOF_IDS::nh];
////
////        // Calculate dy
////        dy_ = (hf_ - ho_) / (double) nh_;
////
////    } // constructor
