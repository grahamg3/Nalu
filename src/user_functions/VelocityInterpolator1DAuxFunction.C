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
  // average velocities
  ux_(0.0),
  uy_(0.0),
  uz_(1.0),
  // max amplitudes
  ax_(0.1),
  ay_(0.1),
  az_(0.1)
{
  //extract the parameters
  if ( params.size() != 6 || params.empty() )
    throw std::runtime_error("VelocityInterpolator1DAuxFunction: requires u, v, w, ax, ay, az"); 
  ux_ = params[0];
  uy_ = params[1];
  uz_ = params[2];
  ax_ = params[3];
  ay_ = params[4];
  az_ = params[5];
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
