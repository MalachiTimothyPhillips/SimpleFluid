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

//-----------------------------------------------------------------------------------------------------------
FluidEquation* FluidEquation::make_fluid_equation(std::string& equationType, std::vector<double>& args){
    if(equationType == "UpwindLinWave"){
        return new UpwindLinWave(args);
    } // FTBS/FTFS first order time, space upwind/downwind shifter method
    if(equationType == "DiffusionEquationFTCS"){
        return new DiffusionEquationFTCS(args);
    } // FTCS first order time, second order space
    if(equationType == "BurgerEquationFTCS"){
        return new BurgerEquationFTCS(args);
    } // FTBS/FTFS first order time, space upwind/downwind shifter with adapative timestep
    if(equationType == "TwoDimLinAdvectionEquation"){
        return new TwoDimLinAdvectionEquation(args);
    }
    if(equationType == "MultiDimDiffusion"){
        return new MultiDimDiffusion(args);
    }
    if(equationType == "MultiDimBurger"){
        return new MultiDimBurger(args);
    }

}


//============================================================================================================
/*
 * Fluid Equation base class
 */
//============================================================================================================
FluidEquation::FluidEquation(std::vector<double> &args) {
    lo_ = args[DOF_IDS::lo];

    lf_ = args[DOF_IDS::lf];
    nl_ = (unsigned int) args[DOF_IDS::nl];
    nt_ = (unsigned int) args[DOF_IDS::nt];

    // Calculate dx
    dx_ = (lf_ - lo_) / (double) nl_;

} // constructor
////-----------------------------------------------------------------------------------------------------------
//void FluidEquation::write_to_file(std::string& template_file_name, unsigned int currentStep){
// //Write current time step to the file
//    std::ofstream outFile(template_file_name + std::to_string(currentStep));
//    for (unsigned int i = 0; i < uSolutions_.size(); ++i) {
//        double pos;
//        convert_idx_to_pos(i, pos);
//        outFile << pos << " ";
//        outFile << uSolutions_[i] << std::endl;
//    } // writes file in two column format: x and u(x)
//    outFile.close();
//}

//-----------------------------------------------------------------------------------------------------------
void FluidEquation::convert_idx_to_pos(unsigned int idx, double &pos) {
    pos = (double) idx * dx_ + lo_;
}

//============================================================================================================
/*
 * Linear Wave Equation, solved using FTBS/FTFS upwind/downwind shifter
 */
//============================================================================================================

UpwindLinWave::UpwindLinWave(std::vector<double>& args) : FluidEquation(args){
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

    for (unsigned int i = 1; i < nl_-1; i++) { // Avoid overwriting the boundary value
        // foward in time, backward/forward in space

        // make sure to use u_(i-1) state at time n, not n+1, so have to save previous value

        uSolutions_[i] = uSolutions_[i] - CFL_ / 2. * (uSolutions_[i + 1] - uSolnm1)
                         + std::abs(CFL_ / 2.) * (uSolutions_[i + 1] - 2 * uSolutions_[i] + uSolnm1);

        // write new position, will be i-1 at next pass
        uSolnm1 = uSolutions_[i];
    }
}

//------------------------------------------------------------------------------------------------------------
void UpwindLinWave::write_to_file(std::string &template_file_name, unsigned int currentStep) {
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

//============================================================================================================
/*
 * Diffusion Equation, solved in FTCS
 */
//============================================================================================================
DiffusionEquationFTCS::DiffusionEquationFTCS(std::vector<double>& args) : FluidEquation(args){
        // Resize uSolutions_ to appropriate size
        uSolutions_.assign(nl_, 0.0);
        nu_ = args[DOF_IDS::c];
        tf_ = args[DOF_IDS::tf];
        dt_ = tf_/(double) nt_;
        // Compute alpha
        alpha_ = nu_*dt_/dx_/dx_; // For stability, alpha <= 0.5
} // constructor

//------------------------------------------------------------------------------------------------------------
void DiffusionEquationFTCS::apply_step() {
    // need to store u_(i-1) at time n
    double uSolnm1 = uSolutions_[0];

    for (unsigned int i = 1; i < nl_ - 1; i++) { // ensure it won't go too far
        // foward in time, central in space
        uSolutions_[i] = (1.0 - 2.0 * alpha_) * uSolutions_[i] +
                         alpha_ * (uSolutions_[i + 1] + uSolnm1);
        // values at end points are handled by the specified boundary conditions

        // update previous state
        uSolnm1 = uSolutions_[i];
    }
}

//------------------------------------------------------------------------------------------------------------
void DiffusionEquationFTCS::write_to_file(std::string &template_file_name, unsigned int currentStep) {
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

//============================================================================================================
/*
 * Burger Equation, Solved FTBS/FTCS in upwind/downwindshift with adaptive time step
 */
//============================================================================================================

BurgerEquationFTCS::BurgerEquationFTCS(std::vector<double>& args) : FluidEquation(args){
    uSolutions_.assign(nl_, 0.0);
    eps_ = args[DOF_IDS::eps];
}
void BurgerEquationFTCS::apply_step() {

    // Record u(i-1) at time step n
    double uSolnm1 = uSolutions_[0];

    // Calculate the time step needed to ensure stability, using the numeric limit eps
    // In order to do this, need the maximum speed at the current step
    // This operation is not cheap, and is O(nx)
    double Cmax = *(std::max_element(std::begin(uSolutions_), std::end(uSolutions_)));

    dt_ = dx_ * eps_ / Cmax;

    for (unsigned int i = 1; i < nl_-1; i++) { // May be forward, do not mess with boundary value
        // foward in time, backward/forward in space
        double CFL = uSolutions_[i] * dt_/dx_;
        uSolutions_[i] = uSolutions_[i] - CFL / 2. * (uSolutions_[i + 1] - uSolnm1)
                         + std::abs(CFL / 2.) * (uSolutions_[i + 1] - 2 * uSolutions_[i] + uSolnm1);

        // write new position, will be i-1 at next pass
        uSolnm1 = uSolutions_[i];
    }

    // Add dt_ to find current time
    T_ += dt_;
}

//------------------------------------------------------------------------------------------------------------
void BurgerEquationFTCS::write_to_file(std::string &template_file_name, unsigned int currentStep) {
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

//============================================================================================================
/*
 * ViscousBurger Equation, Solved FTBS/FTCS in upwind/downwindshift with adaptive time step
 */
//============================================================================================================

ViscousBurgerEquationFTCS::ViscousBurgerEquationFTCS(std::vector<double>& args) : FluidEquation(args){
    uSolutions_.resize(nl_, 0.0);
    eps_ = args[DOF_IDS::eps];
    nu_ = args[DOF_IDS::c];
} // constructor

//------------------------------------------------------------------------------------------------------------
void ViscousBurgerEquationFTCS::write_to_file(std::string &template_file_name, unsigned int currentStep) {
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

//------------------------------------------------------------------------------------------------------------
void ViscousBurgerEquationFTCS::apply_step() {

    // Record u(i-1) at time step n
    double uSolnm1 = uSolutions_[0];

    // Calculate the time step needed to ensure stability, using the numeric limit eps
    // In order to do this, need the maximum speed at the current step
    // This operation is not cheap, and is O(nx)
    double Cmax = *(std::max_element(std::begin(uSolutions_), std::end(uSolutions_)));
    double dt_adv = dx_ * eps_/Cmax; // dt as calculated from the advective portion
    double dt_diff = 0.5 * eps_/nu_ * dx_ * dx_; // dt as calculated from the diffusive portion
    dt_ = std::min(dt_adv,dt_diff); // The correct dt is the minimum - depends on parameters if
    // the diffusive of the advective portion dominates the stability criterion

    for (unsigned int i = 1; i < nl_-1; i++) { // Central, so do not overwrite boundary value
        // foward in time, backward/forward in space
        double CFL = uSolutions_[i] * dt_/dx_; // CFL number
        uSolutions_[i] = uSolutions_[i] - CFL / 2. * (uSolutions_[i + 1] - uSolnm1)
                         + std::abs(CFL / 2.) * (uSolutions_[i + 1] - 2 * uSolutions_[i] + uSolnm1); // Advective
                                                                                                     // Portion
        double alpha = nu_ * dt_ /dx_/dx_; // Diffusive von Neumann number
        uSolutions_[i] += (1.0 - 2.0 * alpha) * uSolutions_[i] +
                         alpha * (uSolutions_[i + 1] + uSolnm1); // Diffusive portion
        // write new position, will be i-1 at next pass
        uSolnm1 = uSolutions_[i];
    }

    // Add dt_ to find current time
    T_ += dt_;
}

//============================================================================================================
/*
 * Two dimensional fluid equation base class
 */
//============================================================================================================
TwoDimFluidEquation::TwoDimFluidEquation(std::vector<double> &args) : FluidEquation(args) {
    ho_ = args[DOF_IDS::ho];
    hf_ = args[DOF_IDS::hf];
    nh_ = (unsigned int) args[DOF_IDS::nh];

    // Calculate dy
    dy_ = (hf_ - ho_) / (double) nh_;

    // Resize solutions vector
    uSolutionsMatrix_.resize(nl_, std::vector<double>(nh_));
} // constructor

//------------------------------------------------------------------------------------------------------------
void TwoDimFluidEquation::write_to_file(std::string &template_file_name, unsigned int currentStep) {
    //Write current time step to the file
    std::ofstream outFile(template_file_name + std::to_string(currentStep));
    for (unsigned int i = 0; i < nl_; ++i) {
        double pos_x;
        convert_idx_to_pos(i, pos_x);
        for (unsigned int j = 0; j < nh_; ++j){
            double pos_y;
            convert_idx_to_pos_y(j,pos_y);
            outFile << pos_x << "     ";
            outFile << pos_y << "     ";
            outFile << uSolutionsMatrix_[i][j] << std::endl;
        }
    } // writes file in two column format: x and u(x)
    outFile.close();
}

//-----------------------------------------------------------------------------------------------------------
void TwoDimFluidEquation::convert_idx_to_pos_y(unsigned int idx, double &pos) {
    pos = (double) idx * dy_ + ho_;
}

//============================================================================================================
/*
 * Linear Advection Equation, FTBS/FTFS
 */
//============================================================================================================

TwoDimLinAdvectionEquation::TwoDimLinAdvectionEquation(std::vector<double>& args) : TwoDimFluidEquation(args){
    c_ = args[DOF_IDS::c];
    tf_ = args[DOF_IDS::tf];

    // Calculate the dt, CFL
    dt_ = tf_ / (double) nt_;
}

//------------------------------------------------------------------------------------------------------------
void TwoDimLinAdvectionEquation::write_to_file(std::string &template_file_name, unsigned int currentStep) {
    TwoDimFluidEquation::write_to_file(template_file_name, currentStep);
}

//------------------------------------------------------------------------------------------------------------
void TwoDimLinAdvectionEquation::apply_step() {

        // save previous column in holder
    std::vector<double> previous_col (nh_, 0.0);
    for (unsigned int i = 0 ; i < nh_; ++i){
        previous_col[i] = uSolutionsMatrix_[0][i];
    }

    for (unsigned int  i = 1 ; i < nl_-1; ++i){
        double uSolnym1 = uSolutionsMatrix_[i][0];
        for (unsigned int j = 1 ; j < nh_-1; ++j){
            double uij = uSolutionsMatrix_[i][j];
            double uip1j = uSolutionsMatrix_[i+1][j];
            double uim1j = previous_col[j];
            double uijp1 = uSolutionsMatrix_[i][j+1];
            double uijm1 = uSolnym1;
            double fx = 0.5*dt_/dx_;
            double fy = 0.5*dt_/dy_;

            // Change to upshift/downshift equation form
            uSolutionsMatrix_[i][j] = uij - c_ * fx * (uip1j - uim1j) + fabs(c_) * fx * (uip1j - 2. * uij + uim1j)
                    - c_ * fy * (uijp1 - uijm1) + fabs(c_) * fy * (uijp1 - 2. * uij + uijm1);
            uSolnym1 = uSolutionsMatrix_[i][j];
        }

        //rewrite the previous_col
        for (unsigned int j = 0 ; j < nh_; ++j){
            previous_col[j] = uSolutionsMatrix_[i][j];
        }

    }
     // Need to save the (n) state so as to not use (n+1) in solutions

    // Need to write out solutions, but for now, leave that out
}

//============================================================================================================
/*
 * Non Linear Advection Equation, FTBS/FTFS
 */
//============================================================================================================
MultiDimNonLinAdvEqn::MultiDimNonLinAdvEqn(std::vector<double>& args) : TwoDimFluidEquation(args){
    eps_ = args[DOF_IDS::eps];

    // set size of vSolutionsMatrix_
    vSolutionsMatrix_.resize(nl_, std::vector<double>(nh_));
}

//------------------------------------------------------------------------------------------------------------
void MultiDimNonLinAdvEqn::apply_step(){

    /*
     * calculate current dt needed, new current time after update
     */

    // Find the maximum U and V, this is actually costly
    std::vector<double> maxUs;
    std::vector<double> maxVs;
    for (unsigned int i = 0 ; i < nl_; ++i){
        maxUs.push_back(*std::max_element(std::begin(uSolutionsMatrix_[i]), std::end(uSolutionsMatrix_[i])));
        maxVs.push_back(*std::max_element(std::begin(vSolutionsMatrix_[i]), std::end(vSolutionsMatrix_[i])));
    }
    double maxU = *std::max_element(std::begin(maxUs), std::end(maxUs));
    double maxV = *std::max_element(std::begin(maxVs), std::end(maxVs));

    double maxC = std::max(maxU, maxV); // Largest value in wave speed

    /*
     * Have Udt/dx + Udt/dy = SF -> dt = SF/(U/dx + U/dy)
     */

    dt_ = eps_/(maxC/dx_ + maxC/dy_); // Eps is safety factor, typically 0.9

    // start by looping over x-variable (space)

    // save previous column in holder
    std::vector<double> previous_colU (nh_, 0.0);
    std::vector<double> previous_colV (nh_, 0.0);
    for (unsigned int i = 0 ; i < nh_; ++i){
        previous_colU[i] = uSolutionsMatrix_[0][i];
        previous_colV[i] = vSolutionsMatrix_[0][i];
    }

    for (unsigned int  i = 1 ; i < nl_-1; ++i){
        double uSolnym1 = uSolutionsMatrix_[i][0];
        double vSolynm1 = vSolutionsMatrix_[i][0];
        for (unsigned int j = 1 ; j < nh_-1; ++j){
            // Solve U, V into temporary holders first
            double temp_u_at_ij;
            double temp_v_at_ij;

            // For easier reference
            double uij = uSolutionsMatrix_[i][j];
            double uip1j = uSolutionsMatrix_[i+1][j];
            double uim1j = previous_colU[j];
            double uijp1 = uSolutionsMatrix_[i][j+1];
            double uijm1 = uSolnym1;
            double vij = vSolutionsMatrix_[i][j];
            double vip1j = vSolutionsMatrix_[i+1][j];
            double vim1j = previous_colV[j];
            double vijp1 = vSolutionsMatrix_[i][j+1];
            double vijm1 = vSolynm1;
            double fx = 0.5*dt_/dx_;
            double fy = 0.5*dt_/dy_;

            temp_u_at_ij = uij - uij * fx * (uip1j - uim1j) + fabs(uij) * fx * (uip1j - 2.*uij + uim1j)
                    - vij * fy * (uijp1 - uijm1) + fabs(vij) * fy * (uijp1 - 2.*uij + uijm1);

            temp_v_at_ij = vij - uij * fx * (vip1j - vim1j) + fabs(uij) * fx * (vip1j - 2.*vij + vim1j)
                    - vij * fy * (vijp1 - vijm1) + fabs(vij) * fy * (vijp1 - 2.*vij + vijm1);

            // Write U, V into new form
            uSolutionsMatrix_[i][j] = temp_u_at_ij;
            vSolutionsMatrix_[i][j] = temp_v_at_ij;

            uSolnym1 = uSolutionsMatrix_[i][j];
            vSolynm1 = vSolutionsMatrix_[i][j];

        }

        //rewrite the previous_col
        for (unsigned int j = 0 ; j < nh_; ++j){
            previous_colU[j] = uSolutionsMatrix_[i][j];
            previous_colV[j] = vSolutionsMatrix_[i][j];
        }

    }

    // Add to extra time
    T_ += dt_;
}

//------------------------------------------------------------------------------------------------------------
void MultiDimNonLinAdvEqn::write_to_file(std::string &template_file_name, unsigned int currentStep) {
    //Write current time step to the file
    std::ofstream outFile(template_file_name + "U" + std::to_string(currentStep));
    for (unsigned int i = 0; i < nl_; ++i) {
        double pos_x;
        convert_idx_to_pos(i, pos_x);
        for (unsigned int j = 0; j < nh_; ++j){
            double pos_y;
            convert_idx_to_pos_y(j,pos_y);
            outFile << pos_x << "     ";
            outFile << pos_y << "     ";
            outFile << uSolutionsMatrix_[i][j] << std::endl;
        }
    } // writes file in two column format: x and u(x)
    outFile.close();

    //Write current time step to the file
    std::ofstream outVFile(template_file_name + "V" + std::to_string(currentStep));
    for (unsigned int i = 0; i < nl_; ++i) {
        double pos_x;
        convert_idx_to_pos(i, pos_x);
        for (unsigned int j = 0; j < nh_; ++j){
            double pos_y;
            convert_idx_to_pos_y(j,pos_y);
            outVFile << pos_x << "     ";
            outVFile << pos_y << "     ";
            outVFile << vSolutionsMatrix_[i][j] << std::endl;
        }
    } // writes file in two column format: x and u(x)
    outVFile.close();
}

//============================================================================================================
/*
 * 2D Diffusion Equation
 */
//============================================================================================================
MultiDimDiffusion::MultiDimDiffusion(std::vector<double>& args) : TwoDimFluidEquation(args){
    nu_ = args[DOF_IDS::c];
    tf_ = args[DOF_IDS::tf];
    dt_ = tf_/(double) nt_;
} // constructor

//------------------------------------------------------------------------------------------------------------
void MultiDimDiffusion::apply_step(){
    // save previous column in holder
    std::vector<double> previous_col (nh_, 0.0);
    for (unsigned int i = 0 ; i < nh_; ++i){
        previous_col[i] = uSolutionsMatrix_[0][i];
    }

    for (unsigned int  i = 1 ; i < nl_-1; ++i){
        double uSolnym1 = uSolutionsMatrix_[i][0];
        for (unsigned int j = 1 ; j < nh_-1; ++j){
            double uij = uSolutionsMatrix_[i][j];
            double uip1j = uSolutionsMatrix_[i+1][j];
            double uim1j = previous_col[j];
            double uijp1 = uSolutionsMatrix_[i][j+1];
            double uijm1 = uSolnym1;
            double fx = nu_*dt_/(dx_*dx_);
            double fy = nu_*dt_/(dy_*dy_);

            uSolutionsMatrix_[i][j] = uij + fx * (uip1j - 2. * uij + uim1j)
                    + fy * (uijp1 - 2. * uij + uijm1);


            uSolnym1 = uSolutionsMatrix_[i][j];
        }

        //rewrite the previous_col
        for (unsigned int j = 0 ; j < nh_; ++j){
            previous_col[j] = uSolutionsMatrix_[i][j];
        }
    }
}

//------------------------------------------------------------------------------------------------------------
void MultiDimDiffusion::write_to_file(std::string &template_file_name, unsigned int currentStep) {
    TwoDimFluidEquation::write_to_file(template_file_name, currentStep);
}

//============================================================================================================
/*
 * 2D Burger Equation (Non Linear Advective Equation + Diffusive Equation)
 */
//============================================================================================================
MultiDimBurger::MultiDimBurger(std::vector<double>& args) : TwoDimFluidEquation(args){
    eps_ = args[DOF_IDS::eps];
    nu_ = args[DOF_IDS::c];
} // constructor

//------------------------------------------------------------------------------------------------------------
void MultiDimBurger::apply_step(){

    /*
     * calculate current dt needed, new current time after update
     */

    // Find the maximum U and V, this is actually costly
    std::vector<double> maxUs;
    std::vector<double> maxVs;
    for (unsigned int i = 0 ; i < nl_; ++i){
        maxUs.push_back(*std::max_element(std::begin(uSolutionsMatrix_[i]), std::end(uSolutionsMatrix_[i])));
        maxVs.push_back(*std::max_element(std::begin(vSolutionsMatrix_[i]), std::end(vSolutionsMatrix_[i])));
    }
    double maxU = *std::max_element(std::begin(maxUs), std::end(maxUs));
    double maxV = *std::max_element(std::begin(maxVs), std::end(maxVs));

    double maxC = std::max(maxU, maxV); // Largest value in wave speed

    /*
     * Have Udt/dx + Udt/dy = SF -> dt = SF/(U/dx + U/dy)
     */
    double dt_adv;
    dt_adv = eps_/(maxC/dx_ + maxC/dy_); // eps_ is safety factor, typically 0.9

    // dt also needs to be calculated another means to ensure diffusive part is O.K.
    double dt_diff;
    double den = nu_ * (1./dx_/dx_ + 1./dy_/dy_);
    dt_diff = 0.5*eps_/den;

    // Take minimum dt from list
    dt_ = std::min(dt_diff,dt_adv);

    // start by looping over x-variable (space)

    // save previous column in holder
    std::vector<double> previous_colU (nh_, 0.0);
    std::vector<double> previous_colV (nh_, 0.0);
    for (unsigned int i = 0 ; i < nh_; ++i){
        previous_colU[i] = uSolutionsMatrix_[0][i];
        previous_colV[i] = vSolutionsMatrix_[0][i];
    }

    for (unsigned int  i = 1 ; i < nl_-1; ++i){
        double uSolnym1 = uSolutionsMatrix_[i][0];
        double vSolynm1 = vSolutionsMatrix_[i][0];
        for (unsigned int j = 1 ; j < nh_-1; ++j){
            // Solve U, V into temporary holders first
            double temp_u_at_ij;
            double temp_v_at_ij;

            // For easier reference
            double uij = uSolutionsMatrix_[i][j];
            double uip1j = uSolutionsMatrix_[i+1][j];
            double uim1j = previous_colU[j];
            double uijp1 = uSolutionsMatrix_[i][j+1];
            double uijm1 = uSolnym1;
            double vij = vSolutionsMatrix_[i][j];
            double vip1j = vSolutionsMatrix_[i+1][j];
            double vim1j = previous_colV[j];
            double vijp1 = vSolutionsMatrix_[i][j+1];
            double vijm1 = vSolynm1;
            double fx = 0.5*dt_/dx_;
            double fy = 0.5*dt_/dy_;
            double fxdiff = nu_*dt_/(dx_*dx_);
            double fydiff = nu_*dt_/(dy_*dy_);


            temp_u_at_ij = uij - uij * fx * (uip1j - uim1j) + fabs(uij) * fx * (uip1j - 2.*uij + uim1j)
                           - vij * fy * (uijp1 - uijm1) + fabs(vij) * fy * (uijp1 - 2.*uij + uijm1); // From burger

            temp_u_at_ij += fxdiff * (uip1j - 2. * uij + uim1j) + fydiff * (uijp1 - 2. * uij + uijm1); // From diffusion

            temp_v_at_ij = vij - uij * fx * (vip1j - vim1j) + fabs(uij) * fx * (vip1j - 2.*vij + vim1j)
                           - vij * fy * (vijp1 - vijm1) + fabs(vij) * fy * (vijp1 - 2.*vij + vijm1); // From burger

            temp_v_at_ij += fxdiff * (vip1j - 2. * vij + vim1j) + fydiff * (vijp1 - 2. * vij + vijm1); // From diffusion

            // Merely add on separate portions from before.

            // Write U, V into new form
            uSolutionsMatrix_[i][j] = temp_u_at_ij;
            vSolutionsMatrix_[i][j] = temp_v_at_ij;

            uSolnym1 = uSolutionsMatrix_[i][j];
            vSolynm1 = vSolutionsMatrix_[i][j];

        }

        //rewrite the previous_col
        for (unsigned int j = 0 ; j < nh_; ++j){
            previous_colU[j] = uSolutionsMatrix_[i][j];
            previous_colV[j] = vSolutionsMatrix_[i][j];
        }

    }

    // Add to extra time
    T_ += dt_;
}

//------------------------------------------------------------------------------------------------------------
void MultiDimBurger::write_to_file(std::string& template_file_name, unsigned int currentStep){
    // Copy that from the non linear advective equation -- same form
    //Write current time step to the file
    std::ofstream outFile(template_file_name + "U" + std::to_string(currentStep));
    for (unsigned int i = 0; i < nl_; ++i) {
        double pos_x;
        convert_idx_to_pos(i, pos_x);
        for (unsigned int j = 0; j < nh_; ++j){
            double pos_y;
            convert_idx_to_pos_y(j,pos_y);
            outFile << pos_x << "     ";
            outFile << pos_y << "     ";
            outFile << uSolutionsMatrix_[i][j] << std::endl;
        }
    } // writes file in two column format: x and u(x)
    outFile.close();

    //Write current time step to the file
    std::ofstream outVFile(template_file_name + "V" + std::to_string(currentStep));
    for (unsigned int i = 0; i < nl_; ++i) {
        double pos_x;
        convert_idx_to_pos(i, pos_x);
        for (unsigned int j = 0; j < nh_; ++j){
            double pos_y;
            convert_idx_to_pos_y(j,pos_y);
            outVFile << pos_x << "     ";
            outVFile << pos_y << "     ";
            outVFile << vSolutionsMatrix_[i][j] << std::endl;
        }
    } // writes file in two column format: x and u(x)
    outVFile.close();
}
