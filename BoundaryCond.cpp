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