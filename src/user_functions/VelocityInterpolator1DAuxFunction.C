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
  // ramp
  ramp_(0.0),
  t_ramp_(0.0),
  t_frac_(0.0),
  // params vector
  params_(params)
{
  //extract the parameters
  if ( params.empty() )
    throw std::runtime_error("VelocityInterpolator1DAuxFunction: requires direction, default velocity components, table"); 
    // First input is a number from 1 to 4 defining direction (1=x, 2=y, 3=z, 4=r)
    // Afterward, there is a table formed by alternating the position and the three velocity components
    // e.g. position, component1, component2, component3, position, component1, component2, component3 etc.
    // Second, third, fourth inputs are the values used outside of the range of data provided
    // Fifth, sixth, seventh inputs define the ramp function
    // If interpolation direction is x, y, or z, cartesian components are assumed for velocity (vx, vy, vz)
    // If interpolation direction is r, cylindrical components are assumed for velocity (vr, vtheta, vz)
  x_ = (params[0] == 1);
  y_ = (params[0] == 2);
  z_ = (params[0] == 3);
  r_ = (params[0] == 4);
  minpos_ = params[7];
  maxpos_ = params[params.size() - 4];
  ramp_ = params[4];
  t_ramp_ = params[5];
  t_frac_ = params[6];
  params_ = params;
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
  double frac_ = 1.0
  if(ramp_ == 1.0) {
    if(time <= t_ramp_) {
      double A = 1/(2*t_frac_);
      double C = 2*std::atanh(t_frac_)/t_ramp_;
      double W = t_ramp_/2;
      frac_ = A*std::tanh(C*(time-W)) + 0.5;
    }
  }
  
  for(unsigned p=0; p < numPoints; ++p) {
    
    // initialize variables at point
    double v0 = params_[1];
    double v1 = params_[2];
    double v2 = params_[3];
    const double xp = coords[0];
    const double yp = coords[1];
    const double zp = coords[2];
    const double radius = std::sqrt(xp*xp + yp*yp);
    const double theta = std::atan2(yp, xp);

    // find the value of the direction used for index
    const double index = x_*xp + y_*yp + z_*zp + r_*(radius);
    
    // interpolate velocities using index
    if ( index < minpos_ || index > maxpos_ ) {
      v0 = params_[1];
      v1 = params_[2];
      v2 = params_[3];
    } else {
      for(unsigned n=7; n < params_.size(); n += 4) {
        if(index == params_[n]) {
          v0 = params_[n+1];
          v1 = params_[n+2];
          v2 = params_[n+3];
        } else if(index > params_[n] && index < params_[n+4]) {
          double ratio = (index - params_[n])/(params_[n+4] - params_[n]);
          v0 = params_[n+1] + ratio*(params_[n+5] - params_[n+1]);
          v1 = params_[n+2] + ratio*(params_[n+6] - params_[n+2]);
          v2 = params_[n+3] + ratio*(params_[n+7] - params_[n+3]);
        }
      }
    }
    
    // transform results into correct coordinates, scale, and apply to field
    if ( r_ == 1 ) {
      fieldPtr[0] = frac_*(v0*cos(theta) - v1*std::sin(theta));
      fieldPtr[1] = frac_*(v0*sin(theta) + v1*std::cos(theta));
      fieldPtr[2] = frac_*v2;
    } else {
      fieldPtr[0] = frac_*v0;
      fieldPtr[1] = frac_*v1;
      fieldPtr[2] = frac_*v2;
    }

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
