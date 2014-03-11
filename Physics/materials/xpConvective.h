/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#ifndef __XPTHERMOCONVECT_
#define __XPTHERMOCONVECT_
#include "xpConductive.h"

namespace xfem
{

  /// This class only adds a velocity framework in the xpConductive material.  This material can be used with an xConductive formulation
  /*!
    It supposes that the velocity is known. It only computes the characteristic time needed in the SUPG formulation for convection as defined by [Brooks and Hughes 1982].
    This material need a material file that contains:
    THERMAL_SOURCE:s
    THERMIC_CONDUCTIVITY : lambda
    THERMIC_CAPACITY
   */
class xpConvective : virtual public xpConductive 
{
public:

  xpConvective(): xpConductive()
  { 
    variables_signature.register_scalar("characteristic_time");
    variables_signature.register_scalar("element_size");
    variables_signature.register_vector("velocity");
    properties_signature.register_scalar("THERMIC_CAPACITY");
    properties.setSignature(&properties_signature);
    properties.astring("MATERIAL_CLASS") = "MATERIAL_THERMO_CONVECT";
  }

  virtual void computeCurrentState()
  {
    xpConductive::computeCurrentState();

    xVector V = curr->vector("velocity");
    double k = properties.scalar("THERMIC_CONDUCTIVITY");
    double h = curr->scalar("element_size");
    double normV = V.mag();
    double Pe = h * normV / (2*k);//Peclet number
    curr->scalar("characteristic_time") = h/(2.*normV) * (1./tanh(Pe) - 1./Pe);
    // cf Brooks and Hughes 1982 and the SUPG method
  }
};

} 
#endif


