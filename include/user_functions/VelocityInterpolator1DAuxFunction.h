/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VelocityInterpolator1DAuxFunction_h
#define VelocityInterpolator1DAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class VelocityInterpolator1DAuxFunction : public AuxFunction
{
public:

  VelocityInterpolator1DAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    const std::vector<double> &params);

  virtual ~VelocityInterpolator1DAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
private:
  bool x_;
  bool y_;
  bool z_;
  bool r_;
  double minpos_;
  double maxpos_;
};

} // namespace nalu
} // namespace Sierra

#endif
