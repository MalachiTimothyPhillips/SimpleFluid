//
// Created by Malachi Phillips on 11/10/16.
//

#include "FluidEquation.h"

//============================================================================================================
/*
 * Fluid Equation base class
 */
//============================================================================================================
FluidEquation::FluidEquation(std::vector<double>& args){
    lo_ = args[DOF_IDS::lo];
    lf_ = args[DOF_IDS::lf];
    nl_ = args[DOF_IDS::nl];
    nt_ = args[DOF_IDS::nt];

    // Calculate dx
    dx_ = (lf_-lo_)/(double)nl_;

} // constructor

//============================================================================================================
/*
 * One dimensional fluid equation base class
 */
//============================================================================================================

// Constructor hidden away in header file

//============================================================================================================
/*
 * Linear Wave Equation Base Class
 */
//============================================================================================================

// Constructor hidden away in header file

//============================================================================================================
/*
 * Linear Wave Equation, solved using FTBS/FTFS upwind/downwind shifter
 */
//============================================================================================================

// Constructor ommitted - same as LinearWaveEquation

//------------------------------------------------------------------------------------------------------------
void UpwindLinWave::apply_step(){
    // Apply single step
    // Do not include boundary points in solution
    double uSolnm1 = uSolutions_[0]; // Store solution of u_(i-1) at time n

    for (unsigned int i = 1; i < nl_; i++){
        // foward in time, backward in space

        // make sure to use u_(i-1) state at time n, not n+1, so have to save previous value

        uSolutions_[i] = uSolutions_[i] - CFL_/2. * (uSolutions_[i+1] - uSolnm1)
                         + std::abs(CFL_/2.) * (uSolutions_[i+1] - 2 * uSolutions_[i] + uSolnm1);

        // write new position, will be i-1 at next pass
        uSolnm1 = uSolutions_[i];
    }
}

//============================================================================================================
/*
 * One dimensional diffusion equation - base case
 */
//============================================================================================================

// Constructor hidden in away in header file

//============================================================================================================
/*
 * Diffusion Equation, solved in FTCS
 */
//============================================================================================================

void DiffusionEquationFTCS::apply_step(){
    // need to store u_(i-1) at time n
    double uSolnm1 = uSolutions_[0];

    for (unsigned int i = 1; i < nl_-1; i++){ // ensure it won't go too far
        // foward in time, central in space
        uSolutions_[i] = (1.0 - 2.0 * alpha_) * uSolutions_[i] +
                         alpha_ * (uSolutions_[i+1] + uSolnm1);
        // values at end points are handled by the specified boundary conditions

        // update previous state
        uSolnm1 = uSolutions_[i];
    }
}

//============================================================================================================
/*
 * Burger Equation Base Class
 */
//============================================================================================================

// Constructor hidden away in base class

//============================================================================================================
/*
 * Burger Equation, Solved FTBS/FTCS in upwind/downwindshift
 */
//============================================================================================================
void BurgerEquationFTCS::apply_step(){

}

//============================================================================================================
/*
 * Two dimensional fluid equation base class
 */
//============================================================================================================
TwoDimFluidEquation::TwoDimFluidEquation(std::vector<double>& args) : FluidEquation(args){
    ho_ = args[DOF_IDS::ho];
    hf_ = args[DOF_IDS::hf];
    nh_ = args[DOF_IDS::nh];

    // Calculate dy
    dy_ = (hf_-ho_)/(double) nh_;

} // constructor

