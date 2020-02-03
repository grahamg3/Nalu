/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/RotatingWallAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

RotatingWallAuxFunction::RotatingWallAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  xCentroid_(0.0),
  yCentroid_(0.0),
  angularVelocity_(10.0),
  zVelocity_(0.0)
{
  if ( params.size() != 4 || params.empty() )
    throw std::runtime_error("RotatingWallAuxFunction: requires xCentroid_, yCentroid_, angularVelocity_, and zVelocity_.");
  
  // extract them
  xCentroid_ = params[0];
  yCentroid_ = params[1];
  angularVelocity_ = params[2];
  zVelocity_ = params[3];
}

void
RotatingWallAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    const double x = coords[0] - xCentroid_;
    const double y = coords[1] - yCentroid_;

    const double radius = std::sqrt(x*x + y*y);
    const double theta = std::atan2(y, x);
    
    fieldPtr[0] = -angularVelocity_*radius*std::sin(theta);
    fieldPtr[1] = angularVelocity_*radius*std::cos(theta);
    fieldPtr[2] = zVelocity_;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
