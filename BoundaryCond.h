//
// Created by Malachi Phillips on 9/25/16.
//

/*
 * Boundary Conditions
 */

#ifndef CFD_HW_BOUNDARYCOND_H
#define CFD_HW_BOUNDARYCOND_H

#include "RuntimeParameters.h"
#include <string>
#include <vector>
#include <boost/multi_array.hpp>

// Base class
class BoundaryCondition{
public:
    // Factory method -- pass string containing boundary information from run-time commands
    static BoundaryCondition *make_boundary_condition(std::string& conditions);
    virtual void enforce_boundary_conditions(std::vector<double>& u_solutions, std::vector<double>& wallValues) = 0;
protected:
private:
    //virtual ~BoundaryCondition(){};
};

// Derived base class for 2D case
class MultiDimensionBoundaryCondition{
public:
    virtual void enforce_boundary_conditions(std::vector<std::vector<double>>& myArray, RuntimeParamMultiDim runtime ) = 0; // pure virtual: requires definition
protected:
private:
    //virtual ~BoundaryCondition(){};
};

// LHS Constant value wall boundary
class LHSWall : public BoundaryCondition{
public:

    double wallValue_ = 0.0;
    void enforce_boundary_conditions(std::vector<double>& u_solutions, std::vector<double>& wallValues);

protected:

private:
    //virtual ~LHSWall(){};

};

// LHS/RHS constant wall value
class TwoWall : public BoundaryCondition{
public:

    double LHSWallValue_ = 0.0;
    double RHSWallValue_ = 0.0;
    void enforce_boundary_conditions(std::vector<double>& u_solutions, std::vector<double>& wallValues);

protected:

private:
    //virtual ~LHSWall(){};
};

// 2D Boundary Condition
class BoxBoundaryCondition : public MultiDimensionBoundaryCondition{
    void enforce_boundary_conditions(std::vector<std::vector<double>>& myArray, RuntimeParamMultiDim runtime);
    double LHSWall(double y, RuntimeParamMultiDim& runtime);
    double RHSWall(double y, RuntimeParamMultiDim& runtime);

    double TopWall(double x, RuntimeParamMultiDim& runtime);
    double BottomWall(double x, RuntimeParamMultiDim& runtime);

};


// TODO: Generalized boundary conditions
#endif //CFD_HW_BOUNDARYCOND_H
