/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#ifndef __XPCAPCITIVE_H
#define __XPCAPCITIVE_H
#include "xpConductive.h"

namespace xfem
{
  /// This class defines a thermal transient material useful for solving an xmThermicTransient problem. 
  /*!It can be usefull for solving either a transient diffusion (properties may be improved to be temperature dependant)  or a transient convection (Using SUPG).
  This material will be used with an xpThermicTransient formulation. The thermal capacity dependance to temperature is only linear. The conductivity is constant but these may be improved.
  This material will need the following properties:
  THERMIC_CONDUCTIVITY
  REF_THERMIC_CAPACITY
  THERMIC_CAPACITY_SLOPE
  THERMAL_SOURCE
  */
  class xpCapacitive :virtual public xpConductive
  {
  public:
    xpCapacitive();
    
    virtual void checkProperties(); 
    virtual void computeCurrentState();
    virtual void computeCurrentState(string phys_to_update);
  
    void sensitivityTo(const std::string& phys_token, double& sensitivity);
    void sensitivityTo(const std::string& phys_token, xVector& sensitivity);
    void  sensitivityTo(const std::string& phys_token, xTensor2& sensitivity);
    //sensitivity to temp gradient (ie conductivity) is already in the parent class
  };
  
} // end of namespace

#endif





