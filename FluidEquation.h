//
// Created by Malachi Phillips on 11/10/16.
//

#ifndef CFD_HW_FLUIDEQUATION_H
#define CFD_HW_FLUIDEQUATION_H
#include <vector>
#include "RuntimeParameters.h"
#include "BoundaryCond.h"
#include "IntCond.h"
#include <boost/multi_array.hpp>
#include <cassert>

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iostream>
#include <iomanip>

///*
// * Forwrd declarations
//*/
//
//class FluidEquation;

typedef std::vector<std::vector<double>> matrix;

//============================================================================================================
/*
 * Struct holding order of passed parameters, regardless if used
 */
//============================================================================================================
struct DOF_IDS{
    const unsigned static int lo = 0;
    const unsigned static int lf = 1;
    const unsigned static int nl = 2;
    const unsigned static int ho = 3;
    const unsigned static int hf = 4;
    const unsigned static int nh = 5;
    const unsigned static int tf = 6;
    const unsigned static int nt = 7;
    const unsigned static int c = 8;
    const unsigned static int eps = 9;
    const unsigned static int wfreq = 10;
};

//============================================================================================================
 /* Base class for equation
 *
 * An equation holds information about the PDE and the method of solution
 * ie. ut + cux = 0 is a PDE, but has several methods of being solved
 * Each method would have its own equation (although they can easily borrow off others)
 */
//============================================================================================================
class FluidEquation{
public:
    // Constructor
    FluidEquation(std::vector<double>& args);
    FluidEquation* make_fluid_equation(std::string &equationType, std::vector<double>& args);
    virtual void apply_step() = 0; // Pure virtual, must be implemented
    virtual void write_to_file(std::string& template_file_name, unsigned int currentStep) = 0;
    void convert_idx_to_pos(unsigned int idx, double &pos);

    // getters -- needed by the procedure
    unsigned int get_nt(){return nt_;};
    unsigned int get_nl(){return nl_;};
    double get_lo(){return lo_;};
    double get_dx(){return dx_;};
    double get_lf(){return lf_;};

    std::vector<double> uSolutions_; // For 1D case
    matrix uSolutionsMatrix_; // For 2D case
    matrix vSolutionsMatrix_; // For 2D paired case


    /*
     * Matrix solutions for lid-driven cavity problem
     */
    matrix phi_; // streams
    matrix omega_; // vorticity
    matrix p_; // pressure
    matrix w_;

protected:
    // Base class does not include solutions (difference between 1D and 2D equations)


    // Data shared by all solution procedures
    double lo_;
    double lf_;
    unsigned int nl_;
    unsigned int nt_;
    double dx_;

private:
};

//============================================================================================================
/*
 * One dimensional linear advection equation solved by upwind/downwind FTBS/FTFS method (first ordr)
 */
//============================================================================================================
class UpwindLinWave : public FluidEquation{
public:
    // constructor
    UpwindLinWave(std::vector<double>& args);

    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);
protected:

    // Data needed for linear wave equation
    double c_;
    double tf_;
    double dt_;
    double CFL_;
private:
};

//============================================================================================================
/*
 * DiffusionEquation, solve by FTCS
 */
//============================================================================================================
class DiffusionEquationFTCS : public FluidEquation{
public:
    // constructor
    DiffusionEquationFTCS(std::vector<double>&args);// : DiffusionEquation(args) {}
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);

protected:
    // Data needed for linear wave equation
    double nu_;
    double tf_;
    double dt_;
    double alpha_;

private:
};

//============================================================================================================
/*
 * One dimensional burger equation, FTCS
 */
//============================================================================================================
class BurgerEquationFTCS : public FluidEquation{
public:
    // constructor
    BurgerEquationFTCS(std::vector<double>&args);
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);

protected:
    // Data needed for linear wave equation
    double dt_;
    double eps_;
    double T_ = 0;
private:
};

//============================================================================================================
/*
 * One dimensional viscous burger equation, FTCS
 */
//============================================================================================================
class ViscousBurgerEquationFTCS : public FluidEquation{
public:
    // constructor
    ViscousBurgerEquationFTCS(std::vector<double>&args);
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);

protected:
    double nu_;
    double dt_;
    double eps_;
    double T_ = 0;
private:
};


//============================================================================================================
/*
 * Base class for two dimensional problems
 */
//============================================================================================================
class TwoDimFluidEquation : public FluidEquation{
public:
    // constructor
    TwoDimFluidEquation(std::vector<double>& args);
    virtual void apply_step() = 0;
    void write_to_file(std::string& template_file_name, unsigned int currentStep);
    void convert_idx_to_pos_y(unsigned int idx, double &pos);

    // getters
    unsigned int get_nh(){return nh_;};
    double get_ho(){return ho_;};
    double get_dy(){return dy_;};
    double get_hf(){return hf_;};

protected:

    // Data needed for all 2D test problems
    double ho_;
    double hf_;
    unsigned int nh_;
    double dy_;

private:
};

//============================================================================================================
/*
 * 2D linear advection equation FTFS/FTBS
 */
//============================================================================================================
class TwoDimLinAdvectionEquation : public TwoDimFluidEquation{
public:
    // constructor
    TwoDimLinAdvectionEquation(std::vector<double>& args);
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);

protected:
    // Data needed for linear wave equation
    double c_;
    double tf_;
    double dt_;
private:

};

//============================================================================================================
/*
 *  2D Non-linear coupled advective equation
 */
//============================================================================================================
class MultiDimNonLinAdvEqn : public TwoDimFluidEquation{
public:
    // constructor
    MultiDimNonLinAdvEqn(std::vector<double>& args);
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);

protected:

    double eps_;
    double dt_;
    double T_ = 0;

private:

};

//============================================================================================================
/*
 *  2D Diffusion Equation
 */
//============================================================================================================
class MultiDimDiffusion : public TwoDimFluidEquation{
public:
    //constructor
    MultiDimDiffusion(std::vector<double>& args);
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);

protected:

    double nu_;
    double dt_;
    double tf_;

private:

};


class MultiDimBurger : public TwoDimFluidEquation{
public:
    MultiDimBurger(std::vector<double>& args);
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);
protected:

    double eps_;
    double dt_;
    double T_ = 0;

    double nu_; // for diffusive portion

private:
};

//============================================================================================================
/*
 * 2D Laplacian Solver
 */
//============================================================================================================
class Laplacian : public TwoDimFluidEquation{
public:
    // constructor
    Laplacian(std::vector<double>& args);
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);

protected:
    // Data needed for linear wave equation
    double eps_;
    double beta_;

    bool isConverged_ = false;
private:

};

//============================================================================================================
/*
 * 2D Lid Driven Cavity
 */
//============================================================================================================
class LidDriven : public TwoDimFluidEquation{
public:
    // constructor
    LidDriven(std::vector<double>& args);
    void apply_step();
    void write_to_file(std::string& template_file_name, unsigned int currentStep);

    void apply_stream_func();
    void apply_vorticity_boundary();
    void apply_rhs();
    void update_vorticity();
    void update_velocity();
    void update_pressure();

    matrix pTemp_; //temporary p matrix

protected:
    // Data needed for linear wave equation
    double eps_;
    double beta_;

    double Re_; // Use the time specification for the reynolds number

    // For now, force sub iterations for stream function solution to be the same as overall
    // Force user to pick delta X, delta Y to be the same
    // Force numerical epsilon to be the same for the stream function solution and others

    bool isConverged_ = false;
private:

};
#endif //CFD_HW_FLUIDEQUATION_H
