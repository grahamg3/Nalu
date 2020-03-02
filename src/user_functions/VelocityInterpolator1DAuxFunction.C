/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/VelocityInterpolator1DAuxFunction.h>
#include <algorithm>

//basic c++
#include <cmath>
#include <vector>
#include <stdexcept>
#include <cstdlib>

namespace sierra{
namespace nalu{

VelocityInterpolator1DAuxFunction::VelocityInterpolator1DAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos,endPos),
  // coordinates (this defines which coordinate is used to locate the interpolation)
  x_(0),
  y_(0),
  z_(0),
  r_(0),
  // positions vector
  minpos_(0.0),
  maxpos_(0.0),
  pos_({0.0, 0.0, 0.0}),
  // velocities array
  vel_({{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}})
{
  //extract the parameters
  if ( params.size() != 6 || params.empty() )
    throw std::runtime_error("VelocityInterpolator1DAuxFunction: requires direction, default velocity components, table"); 
    // First input is a number from 1 to 4 defining direction (1=x, 2=y, 3=z, 4=r)
    // Afterward, there is a table formed by alternating the position and the three velocity components
    // e.g. position, component1, component2, component3, position, component1, component2, component3 etc.
    // Second, third, fourth inputs are the values used outside of the range of data provided
    // If interpolation direction is x, y, or z, cartesian components are assumed for velocity (vx, vy, vz)
    // If interpolation direction is r, cylindrical components are assumed for velocity (vr, vtheta, vz)
  x_ = (params[0] == 1);
  y_ = (params[0] == 2);
  z_ = (params[0] == 3);
  r_ = (params[0] == 4);
  minpos_ = params[4];
  maxpos_ = params[params.size() - 4];
  for(unsigned n=0; n < params.size(); ++n) {
  }
}

void
VelocityInterpolator1DAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    // seed the random functions for this node at this time
    std::srand(time*1e8);
    std::srand(std::rand()*(p+1)*time*1e4);

    // calculate velocity components
    fieldPtr[0] = ux_ + 2*((std::rand()-0.1)/RAND_MAX - 0.5)*ax_;
    fieldPtr[1] = uy_ + 2*((std::rand()-0.1)/RAND_MAX - 0.5)*ay_;
    fieldPtr[2] = uz_ + 2*((std::rand()-0.1)/RAND_MAX - 0.5)*az_;

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
