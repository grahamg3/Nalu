/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <RotatingWallAuxFunction.h>
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
  MaxR_(8.31),
  MaxMagnitude_(98.21)
{
  if ( params.size() != 4 || params.empty() )
    throw std::runtime_error("RotatingWallAuxFunction: requires xCentroid_, yCentroid_, MaxR_ and MaxMagnitude_");
  
  // extract them
  xCentroid_ = params[0];
  yCentroid_ = params[1];
  MaxR_ = params[2];
  MaxMagnitude_ = params[3];
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

    const double rnormalized = std::sqrt(x*x + y*y)/MaxR_;
    
    fieldPtr[0] = rnormalized_*MaxMagnitude_*std::cos(0.8496);
    fieldPtr[1] = rnormalized_*MaxMagnitude_*std::sin(0.8496);
    fieldPtr[2] = MaxMagnitude_*0.5274;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
