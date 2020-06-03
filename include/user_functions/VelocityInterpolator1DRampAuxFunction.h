/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VelocityInterpolator1DRampAuxFunction_h
#define VelocityInterpolator1DRampAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class VelocityInterpolator1DRampAuxFunction : public AuxFunction
{
public:

  VelocityInterpolator1DRampAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    const std::vector<double> &params);

  virtual ~VelocityInterpolator1DRampAuxFunction() {}
  
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
  double ramp_;
  double t_ramp_;
  double t_frac_;
  std::vector<double> params_;
};

} // namespace nalu
} // namespace Sierra

#endif
