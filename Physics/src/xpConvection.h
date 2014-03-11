/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPCONVECT_H
#define _XPCONVECT_H

#include "xpThermic.h"
using namespace xfem;

/// defines a steady State thermal problem with convection.
class xpConvection : public xpThermic
{
public :
  xpConvection  (xData * d, int deg=1) : xpThermic(d, deg), SUPG(true)
  { };
  xpConvection (xData * d, xLevelSet& lsmat, int deg=1) : xpThermic(d, lsmat, deg), SUPG(true)
  { };

  /// compute element_size and stores it as an internal variable.
  void computeElementLength();
 
  /// only adds the convection term to the residual computed in the parent class
  virtual void setFunction(lalg::xCSRVector& F);
  
  /// only adds the convection matrix to the parent computation
  virtual void setJacobian(lalg::xCSRMatrix&  J);

  virtual void  exportFields(int details, const std::string& extension, bool binary=false, bool sorted =false);

  //sets a velocity, if not already defined by an external fluid problem.
  void setVelocity(xEvalConstant<xVector> my_eval);
  void setVelocity(double vx, double vy, double vz)
  {
    xEvalConstant<xVector> veloc(xVector(vx,vy,vz));
    setVelocity(veloc);
    return;
  }

  void initialize();

  void setSUPG() {SUPG = true; return;}
  void unsetSUPG() {SUPG = false; return;}

protected:
  bool SUPG;
};
 
#endif
