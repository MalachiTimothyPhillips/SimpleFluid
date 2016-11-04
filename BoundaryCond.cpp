//
// Created by Malachi Phillips on 9/25/16.
//

#include "BoundaryCond.h"

//Factory
BoundaryCondition *BoundaryCondition::make_boundary_condition(std::string& conditions)
{
    return new LHSWall; //Hard code in single boundary condition, for now


}

void LHSWall::enforce_boundary_conditions(std::vector<double>& u_solutions, std::vector<double>& wallValues ){
    /*
     * enforce boundary conditions at wall
     * ie. at x = 0, u = 3, for example (some constant)
     */

    // Simple case: index is always 0 for LHS, at current implementation (needs changing)
    u_solutions[0] = wallValues[0]; // hard code in value at wall, currently 0

}

void TwoWall::enforce_boundary_conditions(std::vector<double>& u_solutions, std::vector<double>& wallValues ){
    /*
     * enforce boundary conditions at wall
     * ie. at x = 0, u = 3, for example (some constant)
     */

    // Simple case: index is always 0 for LHS, at current implementation (needs changing)
    u_solutions[0] = 2; // hard code in value at wall, currently 0

    u_solutions[u_solutions.size()-1] = wallValues[1];

}

void BoxBoundaryCondition::enforce_boundary_conditions(std::vector<std::vector<double>>& myArray, RuntimeParamMultiDim runtime){

    /*
     * Enforce boundary condition along the edges
     */

    // size of myArray has already been determined at this point

    // Draw values along top and bottom boundaries
    for (unsigned int i = 0 ; i < runtime.get_x_iterations(); ++i){
        // convert index along i to an actual x value
        double curr_x = runtime.get_xo() + static_cast<double>(i) * runtime.get_dx();
        myArray[i][0] = TopWall(curr_x, runtime);
        myArray[i][runtime.get_y_iterations()-1] = BottomWall(curr_x, runtime);
    }

    // Draw values along the LHS and RHS of the boundaries
    for (unsigned int i = 0 ; i < runtime.get_y_iterations(); ++i){
        // convert index along i to an actual y value
        double curr_y = runtime.get_yo() + static_cast<double>(i) + runtime.get_dy();
        myArray[0][i] = LHSWall(curr_y, runtime);
        myArray[runtime.get_x_iterations()-1][i] = RHSWall(curr_y, runtime);
    }
}

// For now, fix the boundary condition to at least match
double BoxBoundaryCondition::LHSWall(double y, RuntimeParamMultiDim& runtime){
    //return sin(y);
    return 1.0;
    //return -25.*(y-5.)*(y-5.);
}

double BoxBoundaryCondition::RHSWall(double y, RuntimeParamMultiDim& runtime){
    double L = runtime.get_xf()-runtime.get_xo();
    //return L*L*cos(y) + sin(y);
    //return 2.0;
    //return L*cos(L*y) - (L-5.)*(L-5.)*(y-5.)*(y-5.);
    //return sin(L) * sin(y);
    return 1.0;
}

double BoxBoundaryCondition::BottomWall(double x, RuntimeParamMultiDim& runtime){
    //return x*x;
    //return 3.0;
    //return x-25*(x-5.)*(x-5.);
    return 1.0;
}

double BoxBoundaryCondition::TopWall(double x, RuntimeParamMultiDim& runtime){
    double H = runtime.get_yf() - runtime.get_yo();
    //return x*x*cos(H) + sin(H);
    //return 4.0;
    //return x*cos(x*H) - (x-5.)*(x-5.)*(H-5.)*(H-5.);
    //return sin(H) * sin(x);
    return 1.0;
}