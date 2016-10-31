//
// Created by Malachi Phillips on 9/25/16.
//

#ifndef CFD_HW_RUNTIMEPARAMETERS_H
#define CFD_HW_RUNTIMEPARAMETERS_H

// Controls solution parameters, from user
/*
 * User must input grid spacing (dx), length, and total time.
 * CFL, dt are inferred from these settings to ensure a reasonable PDE solution
 */

#include <string>
#include <vector>


//============================================================================================================
/*
 * Base class for run time parameters
 */
//============================================================================================================
//class RuntimeParameters{
//public:
//    static RuntimeParameters* create_runtime_param(std::string& whichPDE);
//    unsigned int get_space_iterations(){return space_iterations_;};
//    unsigned int get_time_iterations(){return time_iterations_;};
//    double get_wall_value(){return wallValue_;};
//    double get_xo(){return xo_;};
//    double get_xf(){return xf_;};
//    double get_dx(){return dx_;};
//    double get_to(){return to_;};
//    double get_CFL(){return CFL_;};
//    double get_alpha(){return alpha_;};
//    std::vector<double>& get_wall_vals(){return wallValues_;};
//
//    void read_parameters_from_file(std::string& runtime_args) = 0;
//
//protected:
//
//    unsigned int space_iterations_;
//    unsigned int time_iterations_;
//    double wallValue_;
//    double xo_;
//    double xf_;
//    double dx_;
//    double to_;
//    double CFL_;
//    double alpha_;
//    std::vector<double> wallValues_;
//    void set_space_iterations();
//    void set_wall_value();
//    void set_wall_boundary(std::vector<double>& wallVals);
//};

//============================================================================================================
/*
 * Linear wave runtime parameter handler
 */
//============================================================================================================
class RuntimeParametersLinWave{

public:

    // Read in parameters from input file
    void read_parameters_from_file(std::string& runtime_args);

    double get_to(){return to_;};
    double get_tf(){return tf_;};
    double get_dx(){return dx_;};
    double get_dt(){return dt_;};
    double get_CFL(){return CFL_;};
    double get_xo(){return xo_;};
    double get_xf(){return xf_;};
    double get_c(){return c_;};
    double get_wall_value() {return wallVal_;};
    unsigned int get_time_iterations(){return time_iterations_;};
    unsigned int get_space_iterations(){return space_iterations_;};


protected:
    double to_;
    double tf_;
    double dx_;
    double dt_;
    double CFL_;
    double xo_;
    double xf_;
    double wallVal_;
    unsigned int time_iterations_;
    unsigned int space_iterations_;
    double c_; // constant parameter for wave equation
    double safetyFactor_; // Default:: force CFL to 1 for constant coefficient case

    // setters not needed - from runtime parameters
    void set_to(double to){to_ = to;};
    void set_tf(double tf){tf_ = tf;};
    void set_dt(double dt){dt_ = dt;};
    void set_dx(double dx){dx_ = dx;};
    void set_xo(double xo){xo_ = xo;};
    void set_xf(double xf){xf_ = xf;};
    void set_c(double c){c_ = c;};
    void set_safety_factor(double safety_factor);
    void set_wall_boundary(double wallVal){wallVal_ = wallVal;};
    void set_CFL();
    void set_CFL(double CFL) {CFL_ = CFL;};
    void set_time_iterations(unsigned int timeIterations){ time_iterations_ = timeIterations;};
    void set_space_iterations(); // inferred from CFD number

private:

};

//============================================================================================================
/*
 * Diffusion PDE runtime parameter handler
 */
//============================================================================================================

class RuntimeParametersDiffusion{
public:
    // Read in parameters from input file
    void read_parameters_from_file(std::string& runtime_args);

    double get_to(){return to_;};
    double get_tf(){return tf_;};
    double get_dx(){return dx_;};
    double get_dt(){return dt_;};
    double get_CFL(){return CFL_;};
    double get_xo(){return xo_;};
    double get_xf(){return xf_;};
    double get_nu(){return nu_;};
    double get_wall_value() {return wallVal_;};
    unsigned int get_time_iterations(){return time_iterations_;};
    unsigned int get_space_iterations(){return space_iterations_;};
    double get_alpha(){return alpha_;};

    void add_wall_value(double wallValue){wallValues_.push_back(wallValue);};
    std::vector<double>& get_wall_vals(){return wallValues_;};

protected:
    double to_;
    double tf_;
    double dx_;
    double dt_;
    double CFL_;
    double xo_;
    double xf_;
    double wallVal_;
    unsigned int time_iterations_;
    unsigned int space_iterations_;
    double nu_; // constant parameter for wave equation
    double safetyFactor_; // Default:: force CFL to 1 for constant coefficient case
    double alpha_;
    std::vector<double> wallValues_;

    // setters not needed - from runtime parameters
    void set_to(double to){to_ = to;};
    void set_tf(double tf){tf_ = tf;};
    void set_dt();
    void set_dx(double dx){dx_ = dx;};
    void set_xo(double xo){xo_ = xo;};
    void set_xf(double xf){xf_ = xf;};
    void set_nu(double nu){nu_ = nu;};
    void set_safety_factor(double safety_factor){safetyFactor_ = safety_factor;};
    void set_wall_boundary(std::vector<double>& wallVals);
    void set_CFL();
    void set_CFL(double CFL) {CFL_ = CFL;};
    void set_time_iterations(unsigned int time){time_iterations_ = time;};
    void set_time_iterations();
    void set_space_iterations(); // inferred from CFD number
    void set_alpha();
private:

};

//============================================================================================================
/*
 * Burger PDE runtime parameter handler
 */
//============================================================================================================

class RuntimeParametersBurger{
public:
    // Read in parameters from input file
    void read_parameters_from_file(std::string& runtime_args);

    double get_to(){return to_;};
    double get_tf(){return tf_;};
    double get_dx(){return dx_;};
    double get_dt(){return dt_;};
    double get_CFL(){return CFL_;};
    double get_xo(){return xo_;};
    double get_xf(){return xf_;};
    double get_nu(){return nu_;};
    double get_wall_value() {return wallVal_;};
    unsigned int get_time_iterations(){return time_iterations_;};
    unsigned int get_space_iterations(){return space_iterations_;};
    double get_alpha(){return alpha_;};

    void add_wall_value(double wallValue){wallValues_.push_back(wallValue);};
    std::vector<double>& get_wall_vals(){return wallValues_;};

protected:
    double to_;
    double tf_;
    double dx_;
    double dt_;
    double CFL_;
    double xo_;
    double xf_;
    double wallVal_;
    unsigned int time_iterations_;
    unsigned int space_iterations_;
    double nu_; // constant parameter for wave equation
    double safetyFactor_; // Default:: force CFL to 1 for constant coefficient case
    double alpha_;
    double cmax_;
    std::vector<double> wallValues_;

    // setters not needed - from runtime parameters
    void set_to(double to){to_ = to;};
    void set_tf(double tf){tf_ = tf;};
    void set_dt();
    void set_dx(double dx){dx_ = dx;};
    void set_xo(double xo){xo_ = xo;};
    void set_xf(double xf){xf_ = xf;};
    void set_nu(double nu){nu_ = nu;};
    void set_safety_factor(double safety_factor){safetyFactor_ = safety_factor;};
    void set_wall_boundary(std::vector<double>& wallVals);
    void set_CFL();
    void set_CFL(double CFL) {CFL_ = CFL;};
    void set_time_iterations(unsigned int time){time_iterations_ = time;};
    void set_time_iterations();
    void set_space_iterations(); // inferred from CFD number
    void set_alpha();
    void set_cmax(double c){cmax_ = c;};
private:

};

#endif //CFD_HW_RUNTIMEPARAMETERS_H
