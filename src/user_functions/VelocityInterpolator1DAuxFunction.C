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
  // position limits
  minpos_(0.0),
  maxpos_(0.0),

{
  //extract the parameters
  if ( params.empty() )
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
  const unsigned /*endPos*/,
  const std::vector<double> &params)
{
  for(unsigned p=0; p < numPoints; ++p) {
    
    // initialize variables at point
    double v0 = params[1];
    double v1 = params[2];
    double v2 = params[3];
    const double xp = coords[0];
    const double yp = coords[1];
    const double zp = coords[2];
    const double radius = std::sqrt(xp*xp + yp*yp);
    const double theta = std::atan2(yp, xp);

    // find the value of the direction used for index
    const double index = x_*xp + y_*yp + z_*zp + r_*(radius);
    
    // interpolate velocities using index
    if ( index < minpos_ || index > maxpos_ ) {
      v0 = params[1];
      v1 = params[2];
      v2 = params[3];
    } else {
      bool passed = false;
      for(unsigned n=4; n < params.size(); n += 4) {
        if(index == params[n]) {
          v0 = params[n+1];
          v1 = params[n+2];
          v2 = params[n+3];
        } else if(index > params[n] && passed == false) {
          passed = true;
          double ratio = (index - params[n])/(params[n+4] - params[n]);
          v0 = params[n+1] + ratio*(params[n+5] - params[n+1]);
          v1 = params[n+2] + ratio*(params[n+6] - params[n+2]);
          v2 = params[n+3] + ratio*(params[n+7] - params[n+3]);
        }
      }
    }
    
    // transform results into correct coordinates and apply to field
    if ( r_ == 1 ) {
      fieldPtr[0] = v0*cos(theta) - v1*radius*std::sin(theta);
      fieldPtr[1] = v0*sin(theta) + v1*radius*std::cos(theta);
      fieldPtr[2] = v2;
    } else {
      fieldPtr[0] = v0;
      fieldPtr[1] = v1;
      fieldPtr[2] = v2;
    }

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
