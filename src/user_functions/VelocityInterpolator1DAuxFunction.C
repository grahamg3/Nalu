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
  vel_({{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}),
  defaultvel_({0.0,0.0,0.0})
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
  for(unsigned n=4; n < params.size(); n += 4) {
    pos_[n/4-1] = params[n];
    vel_[n/4-1][0] = params[n+1];
    vel_[n/4-1][1] = params[n+2];
    vel_[n/4-1][2] = params[n+3];
  }
  defaultvel_ = {params[1], params[2], params[3]};
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
    
    // initialize variables at point
    double pointvelocity_[3] = defaultvel_;
    const double xp = coords[0];
    const double yp = coords[1];
    const double zp = coords[2];
    const double radius = std::sqrt(xp*xp + yp*yp);
    const double theta = std::atan2(yp, xp);

    // find the value of the direction used for index
    const double index_ = x_*xp + y_*yp + z_*zp + r_*(radius);
    
    // interpolate velocities using index
    if ( index_ < minpos || index > maxpos ) {
      pointvelocity_ = defaultvel_;
    } else {
      bool passed_ = false;
      for(unsigned m=0; m < pos_.size(); ++m) {
        if(index_ == pos_[m]) {
          pointvelocity_ = {vel_[m][0], vel_[m][1], vel_[m][2]};
        } else if(index_ > pos_[m] && passed_ == false) {
          passed_ = true;
          double ratio_ = (index_ - pos_[m])/(pos_[m+1] - pos[m]);
          double comp1_ = vel_[m][0] + ratio_*(vel_[m+1][0] - vel_[m][0]);
          double comp2_ = vel_[m][1] + ratio_*(vel_[m+1][1] - vel_[m][1]);
          double comp3_ = vel_[m][2] + ratio_*(vel_[m+1][2] - vel_[m][2]);
          pointvelocity_ = {comp1_, comp2_, comp3_};
        }
      }
    }
    
    // transform results into correct coordinates and apply to field
    if ( r_ == 1 ) {
      fieldPtr[0] = pointvelocity_[0]*cos(theta) - pointvelocity_[1]*radius*std::sin(theta);
      fieldPtr[1] = pointvelocity_[0]*sin(theta) + pointvelocity_[1]*radius*std::cos(theta);
      fieldPtr[2] = pointvelocity_[2];
    } else {
      fieldPtr[0] = pointvelocity_[0]
      fieldPtr[1] = pointvelocity_[1];
      fieldPtr[2] = pointvelocity_[2];
    }

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
