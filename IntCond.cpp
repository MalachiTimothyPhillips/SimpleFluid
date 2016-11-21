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
InitCond* InitCond::make_initial_condition(std::string& init_cond, std::vector<double>& args)
{
    // Handle initial condition
    if(init_cond == "SinWave"){
        return new SinWave(args);
    }

    if(init_cond == "Curvilinear"){
        return new Curvilinear(args);
    }

    if(init_cond == "ExtendedCurvilinear"){
        return new ExtendedCurvilinear(args);
    }

    if(init_cond == "LaplacianBoundary"){
        return new LaplacianBoundary(args);
    }

}

// Constructor
InitCond::InitCond(std::vector<double>& args){
    lo_ = args[DOF_IDS::lo];
    nl_ = (unsigned int)args[DOF_IDS::nl];
    lf_ = args[DOF_IDS::lf];
    dx_ = (lf_-lo_)/(double) nl_;
}

void InitCond::convert_idx_to_pos(unsigned int idx, double &pos) {
    pos = (double)idx * dx_ + lo_;
}

// SinWave constructor
SinWave::SinWave(std::vector<double>& args) : InitCond(args){
    // Nothing to do here
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

Curvilinear::Curvilinear(std::vector<double>& args) : InitCond(args){
    ho_ = args[DOF_IDS::ho];
    nh_ = (unsigned int)args[DOF_IDS::nh];
    hf_ = args[DOF_IDS::hf];
    dy_ = (hf_ - ho_)/(double)nh_;

}

void Curvilinear::convert_idx_to_pos_y(unsigned int idx, double &pos) {
    pos = (double)idx * dy_ + ho_;
}
void Curvilinear::pos_func(double x, double y, double& u) {
    if ( x < 1.0 && y < 1.0){
        u = 5.0;
    }
    else{
        u = 1.0;
    }
}

void Curvilinear::apply_initial_cond(){
    unsigned int sizeX = nl_;
    unsigned int sizeY = nh_;
    for (unsigned int i = 0 ; i < sizeX; ++i){
        for (unsigned int j = 0 ; j < sizeY; ++j){
            double x_curr;
            convert_idx_to_pos(i, x_curr);
            double y_curr;
            convert_idx_to_pos_y(j, y_curr);

            double waveVal;
            pos_func(x_curr, y_curr, waveVal);
            fluidEquation_->uSolutionsMatrix_[i][j] = waveVal;
        }
    }
}

void Curvilinear::enforce_boundary(){ // By definition, initial condition function must satisfy boundary
                                      // These are the same evaluation
    // Enforce boundary conditions

    // on LHS wall, x = lo_ always, y may vary
    for (unsigned int j = 0 ; j < nh_; j++){
        double curr_x = lo_;
        double curr_y;
        convert_idx_to_pos_y(j,curr_y);
        double waveValue;
        pos_func(curr_x, curr_y, waveValue); // Calculates waveValue
        fluidEquation_->uSolutionsMatrix_[0][j] = waveValue;
    }

    // on RHS wall, x=lf_ always, y may vary
    for (unsigned int j = 0 ; j < nh_; j++){
        double curr_x = lf_;
        double curr_y;
        convert_idx_to_pos_y(j,curr_y);
        double waveValue;
        pos_func(curr_x, curr_y, waveValue); // Calculates waveValue
        fluidEquation_->uSolutionsMatrix_[nl_-1][j] = waveValue;
    }

    // on bottom wall, y=ho_ always, x may vary
    for (unsigned int i = 0 ; i < nl_; ++i){
        double curr_x;
        convert_idx_to_pos(i,curr_x);
        double curr_y = ho_;
        double waveValue;
        pos_func(curr_x, curr_y, waveValue); // calculates the waveValue
        fluidEquation_->uSolutionsMatrix_[i][0] = waveValue;
    }

    // on top wall, y=hf_ always, x may vary
    for (unsigned int i = 0 ; i < nl_; ++i){
        double curr_x;
        convert_idx_to_pos(i,curr_x);
        double curr_y = hf_;
        double waveValue;
        pos_func(curr_x, curr_y, waveValue); // calculates the waveValue
        fluidEquation_->uSolutionsMatrix_[i][nh_-1] = waveValue;
    }
}

ExtendedCurvilinear::ExtendedCurvilinear(std::vector<double>& args) : InitCond(args){
    ho_ = args[DOF_IDS::ho];
    nh_ = (unsigned int)args[DOF_IDS::nh];
    hf_ = args[DOF_IDS::hf];
    dy_ = (hf_ - ho_)/(double)nh_;
}

void ExtendedCurvilinear::enforce_boundary(){ // By definition, initial condition function must satisfy boundary
    // These are the same evaluation
    // Enforce boundary conditions

    // on LHS wall, x = lo_ always, y may vary
    for (unsigned int j = 0 ; j < nh_; j++){
        double curr_x = lo_;
        double curr_y;
        convert_idx_to_pos_y(j,curr_y);
        double waveValue;
        pos_func(curr_x, curr_y, waveValue); // Calculates waveValue
        fluidEquation_->uSolutionsMatrix_[0][j] = waveValue;
        fluidEquation_->vSolutionsMatrix_[0][j] = waveValue;
    }

    // on RHS wall, x=lf_ always, y may vary
    for (unsigned int j = 0 ; j < nh_; j++){
        double curr_x = lf_;
        double curr_y;
        convert_idx_to_pos_y(j,curr_y);
        double waveValue;
        pos_func(curr_x, curr_y, waveValue); // Calculates waveValue
        fluidEquation_->uSolutionsMatrix_[nl_-1][j] = waveValue;
        fluidEquation_->vSolutionsMatrix_[nl_-1][j] = waveValue;
    }

    // on bottom wall, y=ho_ always, x may vary
    for (unsigned int i = 0 ; i < nl_; ++i){
        double curr_x;
        convert_idx_to_pos(i,curr_x);
        double curr_y = ho_;
        double waveValue;
        pos_func(curr_x, curr_y, waveValue); // calculates the waveValue
        fluidEquation_->uSolutionsMatrix_[i][0] = waveValue;
        fluidEquation_->vSolutionsMatrix_[i][0] = waveValue;
    }

    // on top wall, y=hf_ always, x may vary
    for (unsigned int i = 0 ; i < nl_; ++i){
        double curr_x;
        convert_idx_to_pos(i,curr_x);
        double curr_y = hf_;
        double waveValue;
        pos_func(curr_x, curr_y, waveValue); // calculates the waveValue
        fluidEquation_->uSolutionsMatrix_[i][nh_-1] = waveValue;
        fluidEquation_->vSolutionsMatrix_[i][nh_-1] = waveValue;
    }
}

void ExtendedCurvilinear::convert_idx_to_pos_y(unsigned int idx, double &pos) {
    pos = (double)idx * dy_ + ho_;
}

void ExtendedCurvilinear::pos_func(double x, double y, double& u) {
    if ( x < 1.0 && y < 1.0){
        u = 5.0;
    }
    else{
        u = 1.0;
    }
}

void ExtendedCurvilinear::apply_initial_cond(){
    unsigned int sizeX = nl_;
    unsigned int sizeY = nh_;
    for (unsigned int i = 0 ; i < sizeX; ++i){
        for (unsigned int j = 0 ; j < sizeY; ++j){
            double x_curr;
            convert_idx_to_pos(i, x_curr);
            double y_curr;
            convert_idx_to_pos_y(j, y_curr);

            double waveVal;
            pos_func(x_curr, y_curr, waveVal);
            fluidEquation_->uSolutionsMatrix_[i][j] = waveVal;
            fluidEquation_->vSolutionsMatrix_[i][j] = waveVal;
        }
    }
}

LaplacianBoundary::LaplacianBoundary(std::vector<double>& args) : InitCond(args){
    ho_ = args[DOF_IDS::ho];
    nh_ = (unsigned int)args[DOF_IDS::nh];
    hf_ = args[DOF_IDS::hf];
    dy_ = (hf_ - ho_)/(double)nh_;
}

void LaplacianBoundary::apply_initial_cond(){
    // initial condition is just zeros
    unsigned int XSIZE = nl_;
    unsigned int YSIZE = nh_;
    for (unsigned int i = 0 ; i < XSIZE; ++i){
        for (unsigned int j = 0 ; j < YSIZE; ++j){
            fluidEquation_->uSolutionsMatrix_[i][j] = 0.0;
        }
    }
}

void LaplacianBoundary::enforce_boundary(){
    unsigned int XSIZE = nl_;
    unsigned int YSIZE = nh_;

    // fixed, x = 0
    for (unsigned int j = 0 ; j < YSIZE; ++j){
        fluidEquation_->uSolutionsMatrix_[0][j] = 10.;
    }

    // fixed, x = L
    for (unsigned int j = 0 ; j < YSIZE; ++j){
        fluidEquation_->uSolutionsMatrix_[nl_-1][j] = 5.0;
    }

    // fixed, y = 0
    for (unsigned int i = 0 ; i < XSIZE; ++i){
        fluidEquation_->uSolutionsMatrix_[i][0] = 4.0;
    }

    // fixed, y = H
    for (unsigned int i = 0 ; i < XSIZE; ++i){
        fluidEquation_->uSolutionsMatrix_[i][nh_-1] = 2.0;
    }

}

void LaplacianBoundary::convert_idx_to_pos_y(unsigned int idx, double& pos){
    pos = (double)idx * dy_ + ho_;
}


