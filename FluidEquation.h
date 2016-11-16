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
    void write_to_file(std::string& template_file_name, unsigned int currentStep);
    void convert_idx_to_pos(unsigned int idx, double &pos);

    // getters -- needed by the procedure
    unsigned int get_nt(){return nt_;};
    unsigned int get_nl(){return nl_;};
    double get_lo(){return lo_;};
    double get_dx(){return dx_;};
    double get_lf(){return lf_;};

    std::vector<double> uSolutions_;

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

////============================================================================================================
///*
// * One dimensional equations base class
// */
////============================================================================================================
//class OneDimFluidEquation : public FluidEquation{
//public:
//    // Constructor
//    OneDimFluidEquation(std::vector<double>& args);
//    void write_to_file(std::string& template_file_name, unsigned int currentStep); // write to file, perhaps plot?
//    void convert_idx_to_pos(unsigned int idx, double& pos);
//
//protected:
//
////    std::vector<double> uSolutions_;
//
//private:
//};
//
////============================================================================================================
///*
// * One dimensional linear wave equation -- base class
// */
////============================================================================================================
//class LinearWaveEquation : public OneDimFluidEquation{
//public:
//    // constructor
//    LinearWaveEquation(std::vector<double>&args);
//protected:
//
//    // Data needed for linear wave equation
//    double c_;
//    double tf_;
//    double dt_;
//    double CFL_;
//
//private:
//};

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

////============================================================================================================
///*
// * One dimensional diffusion equation - base case
// */
////============================================================================================================
//class DiffusionEquation : public OneDimFluidEquation{
//public:
//    // constructor
//    DiffusionEquation(std::vector<double>&args);// : OneDimFluidEquation(args) {
////        nu_ = args[DOF_IDS::c];
////        tf_ = args[DOF_IDS::tf];
////        dt_ = tf_/(double) nt_;
////        // Compute alpha
////        alpha_ = nu_*dt_/dx_/dx_;
////    }
//
//protected:
//
//    // Data needed for linear wave equation
//    double nu_;
//    double tf_;
//    double dt_;
//    double alpha_;
//
//private:
//};
//
////============================================================================================================
///*
// * DiffusionEquation, solve by FTCS
// */
////============================================================================================================
//class DiffusionEquationFTCS : public DiffusionEquation{
//public:
//    // constructor
//    DiffusionEquationFTCS(std::vector<double>&args);// : DiffusionEquation(args) {}
//    void apply_step();
//    void write_to_file(std::string& template_file_name, unsigned int currentStep);
//
//protected:
//private:
//};
//
////============================================================================================================
///*
// * One dimensional burger equation, base class
// */
////============================================================================================================
//class BurgerEquation : public OneDimFluidEquation{
//public:
//    // constructor
//    BurgerEquation(std::vector<double>&args);// : OneDimFluidEquation(args) {
////        eps_ = args[DOF_IDS::eps];
////    }
//
//protected:
//
//    // Data needed for linear wave equation
//    double dt_;
//    double eps_;
//    double T_ = 0;
//
//private:
//};
//
////============================================================================================================
///*
// * One dimensional burger equation, FTCS
// */
////============================================================================================================
//class BurgerEquationFTCS : public BurgerEquation{
//public:
//    // constructor
//    BurgerEquationFTCS(std::vector<double>&args);// : BurgerEquation(args) {}
//    void apply_step();
//    void write_to_file(std::string& template_file_name, unsigned int currentStep);
//
//protected:
//
//private:
//};
//
////============================================================================================================
///*
// * One dimensional viscous burger equation, base class
// */
////============================================================================================================
//class ViscousBurgerEquation : public OneDimFluidEquation{
//public:
//    // constructor
//    ViscousBurgerEquation(std::vector<double>&args);// : OneDimFluidEquation(args) {}
//
//protected:
//
//    // Data needed for linear wave equation
//    double nu_;
//    double dt_;
//    double eps_;
//
//private:
//};
//
////============================================================================================================
///*
// * One dimensional viscous burger equation, FTCS
// */
////============================================================================================================
//class ViscousBurgerEquationFTCS : public ViscousBurgerEquation{
//public:
//    // constructor
//    ViscousBurgerEquationFTCS(std::vector<double>&args);// : ViscousBurgerEquation(args) {}
//    void apply_step();
//    void write_to_file(std::string& template_file_name, unsigned int currentStep);
//
//protected:
//private:
//};


////============================================================================================================
///*
// * Base class for two dimensional problems
// */
////============================================================================================================
//class TwoDimFluidEquation : public FluidEquation{
//public:
//    // constructor
//    TwoDimFluidEquation(std::vector<double>& args);
//
//protected:
//    matrix uSolutions_; // 2D solutions matrix
//
//    // Data needed for all 2D test problems
//    double ho_;
//    double hf_;
//    unsigned int nh_;
//    double dy_;
//
//private:
//};
//
////============================================================================================================
///*
// * Base class for 2D linear advection equation
// */
////============================================================================================================
//class TwoDimLinAdvectionEquation : public TwoDimFluidEquation{
//public:
//    // constructor
//    TwoDimLinAdvectionEquation(std::vector<double>& args) : TwoDimFluidEquation(args){}
//
//protected:
//private:
//
//};

#endif //CFD_HW_FLUIDEQUATION_H
