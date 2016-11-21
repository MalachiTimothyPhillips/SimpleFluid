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
FluidEquation *FluidEquation::make_fluid_equation(std::string &equationType, std::vector<double> &args) {

    if (equationType == "UpwindLinWave") {
        return new UpwindLinWave(args);
    } // FTBS/FTFS first order time, space upwind/downwind shifter method
    if (equationType == "DiffusionEquationFTCS") {
        return new DiffusionEquationFTCS(args);
    } // FTCS first order time, second order space
    if (equationType == "BurgerEquationFTCS") {
        return new BurgerEquationFTCS(args);
    } // FTBS/FTFS first order time, space upwind/downwind shifter with adapative timestep
    if (equationType == "TwoDimLinAdvectionEquation") {
        return new TwoDimLinAdvectionEquation(args);
    }
    if (equationType == "MultiDimDiffusion") {
        return new MultiDimDiffusion(args);
    }
    if (equationType == "MultiDimBurger") {
        return new MultiDimBurger(args);
    }
    if (equationType == "Laplacian") {
        return new Laplacian(args);
    }
    if (equationType == "LidDriven") {
        return new LidDriven(args);
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

//-----------------------------------------------------------------------------------------------------------
void FluidEquation::convert_idx_to_pos(unsigned int idx, double &pos) {
    pos = (double) idx * dx_ + lo_;
}

//============================================================================================================
/*
 * Linear Wave Equation, solved using FTBS/FTFS upwind/downwind shifter
 */
//============================================================================================================

UpwindLinWave::UpwindLinWave(std::vector<double> &args) : FluidEquation(args) {
    // Resize uSolutions_ to appropriate size
    uSolutions_.assign(nl_, 0.0);
    c_ = args[DOF_IDS::c];
    tf_ = args[DOF_IDS::tf];

    // Calculate the dt, CFL
    dt_ = tf_ / (double) nt_;
    CFL_ = c_ * dt_ / dx_;
}

//------------------------------------------------------------------------------------------------------------
void UpwindLinWave::apply_step() {
    // Apply single step
    // Do not include boundary points in solution
    double uSolnm1 = uSolutions_[0]; // Store solution of u_(i-1) at time n

    for (unsigned int i = 1; i < nl_ - 1; i++) { // Avoid overwriting the boundary value
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
DiffusionEquationFTCS::DiffusionEquationFTCS(std::vector<double> &args) : FluidEquation(args) {
    // Resize uSolutions_ to appropriate size
    uSolutions_.assign(nl_, 0.0);
    nu_ = args[DOF_IDS::c];
    tf_ = args[DOF_IDS::tf];
    dt_ = tf_ / (double) nt_;
    // Compute alpha
    alpha_ = nu_ * dt_ / dx_ / dx_; // For stability, alpha <= 0.5
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

BurgerEquationFTCS::BurgerEquationFTCS(std::vector<double> &args) : FluidEquation(args) {
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

    for (unsigned int i = 1; i < nl_ - 1; i++) { // May be forward, do not mess with boundary value
        // foward in time, backward/forward in space
        double CFL = uSolutions_[i] * dt_ / dx_;
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

ViscousBurgerEquationFTCS::ViscousBurgerEquationFTCS(std::vector<double> &args) : FluidEquation(args) {
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
    double dt_adv = dx_ * eps_ / Cmax; // dt as calculated from the advective portion
    double dt_diff = 0.5 * eps_ / nu_ * dx_ * dx_; // dt as calculated from the diffusive portion
    dt_ = std::min(dt_adv, dt_diff); // The correct dt is the minimum - depends on parameters if
    // the diffusive of the advective portion dominates the stability criterion

    for (unsigned int i = 1; i < nl_ - 1; i++) { // Central, so do not overwrite boundary value
        // foward in time, backward/forward in space
        double CFL = uSolutions_[i] * dt_ / dx_; // CFL number
        uSolutions_[i] = uSolutions_[i] - CFL / 2. * (uSolutions_[i + 1] - uSolnm1)
                         + std::abs(CFL / 2.) * (uSolutions_[i + 1] - 2 * uSolutions_[i] + uSolnm1); // Advective
        // Portion
        double alpha = nu_ * dt_ / dx_ / dx_; // Diffusive von Neumann number
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
        for (unsigned int j = 0; j < nh_; ++j) {
            double pos_y;
            convert_idx_to_pos_y(j, pos_y);
            outFile << pos_x << "     ";
            outFile << pos_y << "     ";
            outFile << uSolutionsMatrix_[i][j] << std::endl;
        }
        // everytime x changes, do a new line
        outFile << std::endl;
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

TwoDimLinAdvectionEquation::TwoDimLinAdvectionEquation(std::vector<double> &args) : TwoDimFluidEquation(args) {
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
    std::vector<double> previous_col(nh_, 0.0);
    for (unsigned int i = 0; i < nh_; ++i) {
        previous_col[i] = uSolutionsMatrix_[0][i];
    }

    for (unsigned int i = 1; i < nl_ - 1; ++i) {
        double uSolnym1 = uSolutionsMatrix_[i][0];
        for (unsigned int j = 1; j < nh_ - 1; ++j) {
            double uij = uSolutionsMatrix_[i][j];
            double uip1j = uSolutionsMatrix_[i + 1][j];
            double uim1j = previous_col[j];
            double uijp1 = uSolutionsMatrix_[i][j + 1];
            double uijm1 = uSolnym1;
            double fx = 0.5 * dt_ / dx_;
            double fy = 0.5 * dt_ / dy_;

            // Change to upshift/downshift equation form
            uSolutionsMatrix_[i][j] = uij - c_ * fx * (uip1j - uim1j) + fabs(c_) * fx * (uip1j - 2. * uij + uim1j)
                                      - c_ * fy * (uijp1 - uijm1) + fabs(c_) * fy * (uijp1 - 2. * uij + uijm1);
            uSolnym1 = uSolutionsMatrix_[i][j];
        }

        //rewrite the previous_col
        for (unsigned int j = 0; j < nh_; ++j) {
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
MultiDimNonLinAdvEqn::MultiDimNonLinAdvEqn(std::vector<double> &args) : TwoDimFluidEquation(args) {
    eps_ = args[DOF_IDS::eps];

    // set size of vSolutionsMatrix_
    vSolutionsMatrix_.resize(nl_, std::vector<double>(nh_));
}

//------------------------------------------------------------------------------------------------------------
void MultiDimNonLinAdvEqn::apply_step() {

    /*
     * calculate current dt needed, new current time after update
     */

    // Find the maximum U and V, this is actually costly
    std::vector<double> maxUs;
    std::vector<double> maxVs;
    for (unsigned int i = 0; i < nl_; ++i) {
        maxUs.push_back(*std::max_element(std::begin(uSolutionsMatrix_[i]), std::end(uSolutionsMatrix_[i])));
        maxVs.push_back(*std::max_element(std::begin(vSolutionsMatrix_[i]), std::end(vSolutionsMatrix_[i])));
    }
    double maxU = *std::max_element(std::begin(maxUs), std::end(maxUs));
    double maxV = *std::max_element(std::begin(maxVs), std::end(maxVs));

    double maxC = std::max(maxU, maxV); // Largest value in wave speed

    /*
     * Have Udt/dx + Udt/dy = SF -> dt = SF/(U/dx + U/dy)
     */

    dt_ = eps_ / (maxC / dx_ + maxC / dy_); // Eps is safety factor, typically 0.9

    // start by looping over x-variable (space)

    // save previous column in holder
    std::vector<double> previous_colU(nh_, 0.0);
    std::vector<double> previous_colV(nh_, 0.0);
    for (unsigned int i = 0; i < nh_; ++i) {
        previous_colU[i] = uSolutionsMatrix_[0][i];
        previous_colV[i] = vSolutionsMatrix_[0][i];
    }

    for (unsigned int i = 1; i < nl_ - 1; ++i) {
        double uSolnym1 = uSolutionsMatrix_[i][0];
        double vSolynm1 = vSolutionsMatrix_[i][0];
        for (unsigned int j = 1; j < nh_ - 1; ++j) {
            // Solve U, V into temporary holders first
            double temp_u_at_ij;
            double temp_v_at_ij;

            // For easier reference
            double uij = uSolutionsMatrix_[i][j];
            double uip1j = uSolutionsMatrix_[i + 1][j];
            double uim1j = previous_colU[j];
            double uijp1 = uSolutionsMatrix_[i][j + 1];
            double uijm1 = uSolnym1;
            double vij = vSolutionsMatrix_[i][j];
            double vip1j = vSolutionsMatrix_[i + 1][j];
            double vim1j = previous_colV[j];
            double vijp1 = vSolutionsMatrix_[i][j + 1];
            double vijm1 = vSolynm1;
            double fx = 0.5 * dt_ / dx_;
            double fy = 0.5 * dt_ / dy_;

            temp_u_at_ij = uij - uij * fx * (uip1j - uim1j) + fabs(uij) * fx * (uip1j - 2. * uij + uim1j)
                           - vij * fy * (uijp1 - uijm1) + fabs(vij) * fy * (uijp1 - 2. * uij + uijm1);

            temp_v_at_ij = vij - uij * fx * (vip1j - vim1j) + fabs(uij) * fx * (vip1j - 2. * vij + vim1j)
                           - vij * fy * (vijp1 - vijm1) + fabs(vij) * fy * (vijp1 - 2. * vij + vijm1);

            // Write U, V into new form
            uSolutionsMatrix_[i][j] = temp_u_at_ij;
            vSolutionsMatrix_[i][j] = temp_v_at_ij;

            uSolnym1 = uSolutionsMatrix_[i][j];
            vSolynm1 = vSolutionsMatrix_[i][j];

        }

        //rewrite the previous_col
        for (unsigned int j = 0; j < nh_; ++j) {
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
        for (unsigned int j = 0; j < nh_; ++j) {
            double pos_y;
            convert_idx_to_pos_y(j, pos_y);
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
        for (unsigned int j = 0; j < nh_; ++j) {
            double pos_y;
            convert_idx_to_pos_y(j, pos_y);
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
MultiDimDiffusion::MultiDimDiffusion(std::vector<double> &args) : TwoDimFluidEquation(args) {
    nu_ = args[DOF_IDS::c];
    tf_ = args[DOF_IDS::tf];
    dt_ = tf_ / (double) nt_;
} // constructor

//------------------------------------------------------------------------------------------------------------
void MultiDimDiffusion::apply_step() {
    // save previous column in holder
    std::vector<double> previous_col(nh_, 0.0);
    for (unsigned int i = 0; i < nh_; ++i) {
        previous_col[i] = uSolutionsMatrix_[0][i];
    }

    for (unsigned int i = 1; i < nl_ - 1; ++i) {
        double uSolnym1 = uSolutionsMatrix_[i][0];
        for (unsigned int j = 1; j < nh_ - 1; ++j) {
            double uij = uSolutionsMatrix_[i][j];
            double uip1j = uSolutionsMatrix_[i + 1][j];
            double uim1j = previous_col[j];
            double uijp1 = uSolutionsMatrix_[i][j + 1];
            double uijm1 = uSolnym1;
            double fx = nu_ * dt_ / (dx_ * dx_);
            double fy = nu_ * dt_ / (dy_ * dy_);

            uSolutionsMatrix_[i][j] = uij + fx * (uip1j - 2. * uij + uim1j)
                                      + fy * (uijp1 - 2. * uij + uijm1);


            uSolnym1 = uSolutionsMatrix_[i][j];
        }

        //rewrite the previous_col
        for (unsigned int j = 0; j < nh_; ++j) {
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
MultiDimBurger::MultiDimBurger(std::vector<double> &args) : TwoDimFluidEquation(args) {
    eps_ = args[DOF_IDS::eps];
    nu_ = args[DOF_IDS::c];
} // constructor

//------------------------------------------------------------------------------------------------------------
void MultiDimBurger::apply_step() {

    /*
     * calculate current dt needed, new current time after update
     */

    // Find the maximum U and V, this is actually costly
    std::vector<double> maxUs;
    std::vector<double> maxVs;
    for (unsigned int i = 0; i < nl_; ++i) {
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
    dt_adv = eps_ / (maxC / dx_ + maxC / dy_); // eps_ is safety factor, typically 0.9

    // dt also needs to be calculated another means to ensure diffusive part is O.K.
    double dt_diff;
    double den = nu_ * (1. / dx_ / dx_ + 1. / dy_ / dy_);
    dt_diff = 0.5 * eps_ / den;

    // Take minimum dt from list
    dt_ = std::min(dt_diff, dt_adv);

    // start by looping over x-variable (space)

    // save previous column in holder
    std::vector<double> previous_colU(nh_, 0.0);
    std::vector<double> previous_colV(nh_, 0.0);
    for (unsigned int i = 0; i < nh_; ++i) {
        previous_colU[i] = uSolutionsMatrix_[0][i];
        previous_colV[i] = vSolutionsMatrix_[0][i];
    }

    for (unsigned int i = 1; i < nl_ - 1; ++i) {
        double uSolnym1 = uSolutionsMatrix_[i][0];
        double vSolynm1 = vSolutionsMatrix_[i][0];
        for (unsigned int j = 1; j < nh_ - 1; ++j) {
            // Solve U, V into temporary holders first
            double temp_u_at_ij;
            double temp_v_at_ij;

            // For easier reference
            double uij = uSolutionsMatrix_[i][j];
            double uip1j = uSolutionsMatrix_[i + 1][j];
            double uim1j = previous_colU[j];
            double uijp1 = uSolutionsMatrix_[i][j + 1];
            double uijm1 = uSolnym1;
            double vij = vSolutionsMatrix_[i][j];
            double vip1j = vSolutionsMatrix_[i + 1][j];
            double vim1j = previous_colV[j];
            double vijp1 = vSolutionsMatrix_[i][j + 1];
            double vijm1 = vSolynm1;
            double fx = 0.5 * dt_ / dx_;
            double fy = 0.5 * dt_ / dy_;
            double fxdiff = nu_ * dt_ / (dx_ * dx_);
            double fydiff = nu_ * dt_ / (dy_ * dy_);


            temp_u_at_ij = uij - uij * fx * (uip1j - uim1j) + fabs(uij) * fx * (uip1j - 2. * uij + uim1j)
                           - vij * fy * (uijp1 - uijm1) + fabs(vij) * fy * (uijp1 - 2. * uij + uijm1); // From burger

            temp_u_at_ij += fxdiff * (uip1j - 2. * uij + uim1j) + fydiff * (uijp1 - 2. * uij + uijm1); // From diffusion

            temp_v_at_ij = vij - uij * fx * (vip1j - vim1j) + fabs(uij) * fx * (vip1j - 2. * vij + vim1j)
                           - vij * fy * (vijp1 - vijm1) + fabs(vij) * fy * (vijp1 - 2. * vij + vijm1); // From burger

            temp_v_at_ij += fxdiff * (vip1j - 2. * vij + vim1j) + fydiff * (vijp1 - 2. * vij + vijm1); // From diffusion

            // Merely add on separate portions from before.

            // Write U, V into new form
            uSolutionsMatrix_[i][j] = temp_u_at_ij;
            vSolutionsMatrix_[i][j] = temp_v_at_ij;

            uSolnym1 = uSolutionsMatrix_[i][j];
            vSolynm1 = vSolutionsMatrix_[i][j];

        }

        //rewrite the previous_col
        for (unsigned int j = 0; j < nh_; ++j) {
            previous_colU[j] = uSolutionsMatrix_[i][j];
            previous_colV[j] = vSolutionsMatrix_[i][j];
        }

    }

    // Add to extra time
    T_ += dt_;
}

//------------------------------------------------------------------------------------------------------------
void MultiDimBurger::write_to_file(std::string &template_file_name, unsigned int currentStep) {
    // Copy that from the non linear advective equation -- same form
    //Write current time step to the file
    std::ofstream outFile(template_file_name + "U" + std::to_string(currentStep));
    for (unsigned int i = 0; i < nl_; ++i) {
        double pos_x;
        convert_idx_to_pos(i, pos_x);
        for (unsigned int j = 0; j < nh_; ++j) {
            double pos_y;
            convert_idx_to_pos_y(j, pos_y);
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
        for (unsigned int j = 0; j < nh_; ++j) {
            double pos_y;
            convert_idx_to_pos_y(j, pos_y);
            outVFile << pos_x << "     ";
            outVFile << pos_y << "     ";
            outVFile << vSolutionsMatrix_[i][j] << std::endl;
        }
    } // writes file in two column format: x and u(x)
    outVFile.close();
}

//============================================================================================================
/*
 * 2D Laplacian Solver
 */
//============================================================================================================
Laplacian::Laplacian(std::vector<double> &args) : TwoDimFluidEquation(args) {
    // set eps_, beta_
    eps_ = args[DOF_IDS::eps];
    beta_ = dx_ / dy_;
}

//------------------------------------------------------------------------------------------------------------
void Laplacian::apply_step() {

    // check if the step needs to be applied
    double err = -100.; //  error term, initialized to -100
    // Temporary Solutions Vector
    matrix tempUSolutions_;
    tempUSolutions_.resize(nl_, std::vector<double>(nh_));
    if (!isConverged_) {// apply step: else do nothing
        // Loop over i, j index first
        // outside loop, calculate constant factor
        double cValue = 0.5 / (1 + beta_ * beta_);
        for (unsigned int i = 1; i < nl_ - 1; ++i) {
            for (unsigned int j = 1; j < nh_ - 1; ++j) {
                double uim1j = uSolutionsMatrix_[i - 1][j];
                double uip1j = uSolutionsMatrix_[i + 1][j];
                double uijp1 = uSolutionsMatrix_[i][j + 1];
                double uijm1 = uSolutionsMatrix_[i][j - 1];
                double uijnp1;
                uijnp1 = uim1j + uip1j + beta_ * beta_ * (uijm1 + uijp1);
                uijnp1 *= cValue; // multiply by constant
                // Write solution
                tempUSolutions_[i][j] = uijnp1;
                // Find err (difference between current time step and previous)
                double curr_err = fabs(tempUSolutions_[i][j] - uSolutionsMatrix_[i][j]);
                // if current error is larger than the error, re-write it
                if (curr_err > err) {
                    err = curr_err; // Ensures that only the greatest value
                }
            }
        }
    }

    // Once step is completed, re-write it to actual solutions vector
    for (unsigned int i = 0; i < nl_; ++i) {
        for (unsigned int j = 0; j < nh_; ++j) {
            uSolutionsMatrix_[i][j] = tempUSolutions_[i][j];
        }
    }

    // Once after step, check if the maximum error term is LESS THAN the error tolerance
    if (!isConverged_ && fabs(eps_) > fabs(err)) { // In case fabs(err) is not properly set
        isConverged_ = true; // solution converged
    }

}

//------------------------------------------------------------------------------------------------------------
void Laplacian::write_to_file(std::string &template_file_name, unsigned int currentStep) {
    TwoDimFluidEquation::write_to_file(template_file_name, currentStep); // Only truly care about last time step
}

//============================================================================================================
/*
 * 2D Lid Driven cavity problem
 */
//============================================================================================================
LidDriven::LidDriven(std::vector<double> &args) : TwoDimFluidEquation(args) {
    // set eps_, beta_
    eps_ = args[DOF_IDS::eps];
    beta_ = 0.0; // unused
    Re_ = args[DOF_IDS::c]; // use c as reynolds number

    dy_ = dx_; // force timestep to be same

    // set space as the same
    nh_ = nl_;
    ho_ = lo_;
    hf_ = lf_;


    // allocate sizes for matrices
    vSolutionsMatrix_.resize(nl_, std::vector<double>(nh_));
    phi_.resize(nl_, std::vector<double>(nh_));
    omega_.resize(nl_, std::vector<double>(nh_));
    p_.resize(nl_, std::vector<double>(nh_));
    w_.resize(nl_, std::vector<double>(nh_));
    pTemp_.resize(nl_, std::vector<double>(nh_));

    // initialize all as 0's
    for (unsigned int i = 0; i < nl_; ++i) {
        for (unsigned int j = 0; j < nh_; ++j) {
            uSolutionsMatrix_[i][j] = 0.;
            vSolutionsMatrix_[i][j] = 0.;
            phi_[i][j] = 0.;
            omega_[i][j] = 0.;
            p_[i][j] = 0.;
            w_[i][j] = 0.;
            pTemp_[i][j] = 0.;
        }
    }

}

//------------------------------------------------------------------------------------------------------------
void LidDriven::apply_step() {

    // Within each sub step, need to complete the following tasks:
    // Stream function calculations
    // Application of Boundary Values for vorticity (from phi's, handled internally)
    // RHS calculations
    // Vorticity update

    // Outside of iterations, calculate U,V fields and P.
    // For now, print out stream function and vorticity

    apply_stream_func();
    apply_vorticity_boundary();
    apply_rhs();
    update_vorticity();
    update_velocity(); // Not needed, except for noting solution progress
    update_pressure(); // Not needed, except for noting solution progress

}

//------------------------------------------------------------------------------------------------------------
void LidDriven::apply_stream_func() {

    double beta = 1.5; // for now, use 1.5

    // Can implement SOR/SUR if beta is changed, but for now, keep the same

    // Loop over total iterations
    double err = 0.0;
    for (unsigned int t = 0; t < nt_; ++t) {

        for (unsigned int i = 1; i < nl_ - 1; ++i) {
            for (unsigned int j = 1; j < nh_ - 1; ++j) {
                double phi_ij = phi_[i][j];
                double phi_ip1j = phi_[i + 1][j];
                double phi_im1j = phi_[i - 1][j];
                double phi_ijp1 = phi_[i][j + 1];
                double phi_ijm1 = phi_[i][j - 1];

                double omega_ij = omega_[i][j];

                phi_[i][j] = 0.25 * beta * (phi_ip1j + phi_im1j + phi_ijp1 + phi_ijm1 + dx_ * dx_ * omega_ij)
                             + (1.0 - beta) * phi_ij;
            }
        }

        // At end, find total errors
        for (unsigned int i = 0; i < nl_; ++i) {
            for (unsigned int j = 0; j < nh_; ++j) {
                err += err + fabs(w_[i][j] - phi_[i][j]);
            }
        }
        if (err < eps_) { // converged
            isConverged_ = true;
            break; // can stop iterations
        }
    }
}

//------------------------------------------------------------------------------------------------------------
void LidDriven::apply_vorticity_boundary() {
    for (unsigned int i = 0; i < nl_; ++i) {
        omega_[i][0] = -2.0 * phi_[i][1] / dx_ / dx_;
        omega_[i][nh_ - 1] = -2.0 * phi_[i][nh_ - 2] / dx_ / dx_ - 2.0 / dx_;
    }
    for (unsigned int j = 0; j < nh_; ++j) {
        omega_[0][j] = -2.0 * phi_[1][j] / dx_ / dx_;
        omega_[nl_ - 1][j] = -2.0 * phi_[nl_ - 2][j] / dx_ / dx_;
    }
}

//------------------------------------------------------------------------------------------------------------
void LidDriven::apply_rhs() {
    for (unsigned int i = 1; i < nl_ - 1; ++i) {
        for (unsigned int j = 1; j < nh_ - 1; ++j) {
            double phi_ip1j = phi_[i + 1][j];
            double phi_im1j = phi_[i - 1][j];
            double phi_ijp1 = phi_[i][j + 1];
            double phi_ijm1 = phi_[i][j - 1];
            double omega_ip1j = omega_[i + 1][j];
            double omega_im1j = omega_[i - 1][j];
            double omega_ijp1 = omega_[i][j + 1];
            double omega_ijm1 = omega_[i][j - 1];
            double omega_ij = omega_[i][j];

            w_[i][j] = -0.25 * ((phi_ijp1 - phi_ijm1) * (omega_ip1j - omega_im1j)
                                - (phi_ip1j - phi_im1j) * (omega_ijp1 - omega_ijm1)) / (dx_ * dx_)
                       + (1 / Re_) * (omega_ip1j + omega_im1j + omega_ijp1 + omega_ijm1 -
                                      4.0 * omega_ij) / (dx_ * dx_);
        }
    }
}

//------------------------------------------------------------------------------------------------------------
void LidDriven::update_vorticity() {
    for (unsigned int i = 1; i < nl_ - 1; ++i) {
        for (unsigned int j = 1; j < nh_ - 1; ++j) {
            omega_[i][j] = omega_[i][j] + 0.01 * w_[i][j]; // Force dt = 0.01
            // Can change later for stability - this time is artificial
        }
    }
}

//------------------------------------------------------------------------------------------------------------
void LidDriven::update_velocity() {
    // Calculate the velocities
    // Edges:
    for (unsigned int i = 0; i < nl_; ++i) {
        uSolutionsMatrix_[i][nh_ - 1] = 1.;
    }
    for (unsigned int j = 0; j < nh_; ++j) {
        vSolutionsMatrix_[nl_ - 1][j] = 0.02;
    }
    for (unsigned int i = 1; i < nl_ - 1; ++i) {
        for (unsigned int j = 0; j < nh_ - 1; ++j) {
            uSolutionsMatrix_[i][j] = (phi_[i][j + 1] - phi_[i][j]) / (2 * dx_);
            vSolutionsMatrix_[i][j] = (phi_[i + 1][j] - phi_[i][j]) / (2 * dx_);
        }
    }

    // Not needed for solutions, but may be plotted with time for solution progresion
}

//------------------------------------------------------------------------------------------------------------
void LidDriven::update_pressure() {
    matrix rhs;
    rhs.resize(nl_, std::vector<double>(nh_));
    for (unsigned int i = 1; i < nl_ - 1; ++i) {
        for (unsigned int j = 1; j < nh_ - 1; ++j) {
            double phi_ij = phi_[i][j];
            double phi_ip1j = phi_[i + 1][j];
            double phi_im1j = phi_[i - 1][j];
            double phi_ijp1 = phi_[i][j + 1];
            double phi_ijm1 = phi_[i][j - 1];
            double phi_ip1jp1 = phi_[i + 1][j + 1];
            double phi_ip1jm1 = phi_[i + 1][j - 1];
            double phi_im1jp1 = phi_[i - 1][j + 1];
            double phi_im1jm1 = phi_[i - 1][j - 1];
            rhs[i][j] = ((phi_im1j - 2 * phi_ij + phi_ip1j) / (dx_ * dx_) *
                         ((phi_ijm1 - 2 * phi_ij + phi_ijp1) / (dx_ * dx_))) -
                        (phi_ip1jp1 - phi_ip1jm1 - phi_im1jp1 + phi_im1jm1) / (4 * (dx_ * dx_));

            double p_ip1j = p_[i + 1][j];
            double p_im1j = p_[i - 1][j];
            double p_ijp1 = p_[i][j + 1];
            double p_ijm1 = p_[i][j - 1];

            pTemp_[i][j] = (0.25 * (p_ip1j + p_im1j + p_ijp1 + p_ijm1)
                            - 0.5 * ((rhs[i][j] * dx_ * dx_ * dx_ * dx_)));
        }
    }

    // rewrite pressure solutions
    for (unsigned int i = 0; i < nl_; ++i) {
        for (unsigned int j = 0; j < nh_; ++j) {
            p_[i][j] = pTemp_[i][j];
        }
    }
}

//------------------------------------------------------------------------------------------------------------
void LidDriven::write_to_file(std::string &template_file_name, unsigned int currentStep) {
    // only output on LAST STEP
    if (currentStep > nt_ - 3) {
        std::ofstream outFile(template_file_name + "PSI" + std::to_string(currentStep));
        for (unsigned int i = 0; i < nl_; ++i) {
            double pos_x;
            convert_idx_to_pos(i, pos_x);
            for (unsigned int j = 0; j < nh_; ++j) {
                double pos_y;
                convert_idx_to_pos_y(j, pos_y);
                outFile << pos_x << "     ";
                outFile << pos_y << "     ";
                outFile << phi_[i][j] << std::endl; // Try to print out stream function
            }
            outFile << std::endl;
        } // writes file in two column format: x and u(x)
        outFile.close();

        std::ofstream outOmFile(template_file_name + "OMEGA" + std::to_string(currentStep));
        for (unsigned int i = 0; i < nl_; ++i) {
            double pos_x;
            convert_idx_to_pos(i, pos_x);
            for (unsigned int j = 0; j < nh_; ++j) {
                double pos_y;
                convert_idx_to_pos_y(j, pos_y);
                outOmFile << pos_x << "     ";
                outOmFile << pos_y << "     ";
                outOmFile << omega_[i][j] << std::endl; // Try to print out stream function
            }
            outOmFile << std::endl;
        } // writes file in two column format: x and u(x)
        outOmFile.close();

        std::ofstream outUFile(template_file_name + "U" + std::to_string(currentStep));
        for (unsigned int i = 0; i < nl_; ++i) {
            double pos_x;
            convert_idx_to_pos(i, pos_x);
            for (unsigned int j = 0; j < nh_; ++j) {
                double pos_y;
                convert_idx_to_pos_y(j, pos_y);
                outUFile << pos_x << "     ";
                outUFile << pos_y << "     ";
                outUFile << uSolutionsMatrix_[i][j] << std::endl; // Try to print out stream function
            }
            outUFile << std::endl;
        } // writes file in two column format: x and u(x)
        outUFile.close();

        std::ofstream outVFile(template_file_name + "V" + std::to_string(currentStep));
        for (unsigned int i = 0; i < nl_; ++i) {
            double pos_x;
            convert_idx_to_pos(i, pos_x);
            for (unsigned int j = 0; j < nh_; ++j) {
                double pos_y;
                convert_idx_to_pos_y(j, pos_y);
                outVFile << pos_x << "     ";
                outVFile << pos_y << "     ";
                outVFile << vSolutionsMatrix_[i][j] << std::endl; // Try to print out stream function
            }
            outVFile << std::endl;
        } // writes file in two column format: x and u(x)
        outVFile.close();

        std::ofstream outPFile(template_file_name + "P" + std::to_string(currentStep));
        for (unsigned int i = 0; i < nl_; ++i) {
            double pos_x;
            convert_idx_to_pos(i, pos_x);
            for (unsigned int j = 0; j < nh_; ++j) {
                double pos_y;
                convert_idx_to_pos_y(j, pos_y);
                outPFile << pos_x << "     ";
                outPFile << pos_y << "     ";
                outPFile << p_[i][j] << std::endl; // Try to print out stream function
            }
            outPFile << std::endl;
        } // writes file in two column format: x and u(x)
        outPFile.close();
    }
}

