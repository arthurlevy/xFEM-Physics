/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#ifndef __XPPOWERLAWFLUID_H
#define __XPPOWERLAWFLUID_H
#include <iostream>
#include <cassert>
#include <string>
#include "AOMDfwd.h"
#include "xTensor2.h"
#include "xTensor3.h"
#include "xTensor4.h"
#include "xTensors.h"
#include "xpNewtonianIncompFluid.h"

namespace xfem
{
 ///This class defines a powerlaw material which is a quasi newtonian material.It works with an xmFluidMechanics formulation
  /*! The carreau material follows the law sigma = eta * D where eta = eta0*(Deq)^(m-1)
This class works with an xmFluidMechanics formulation and needs a material file that contains : 
CARREAU_TIME:lambda
POWER_LAW_INDEX: m
VISCOSITY:eta0
REG_COEFF : constrain Deq^2 not to be smaller than this reg_coeff in order to  avoid infinite stress
  */
class xpPowerLawFluid : virtual public xpNewtonianIncompFluid {

public:
  xpPowerLawFluid();
  void checkProperties(); 
  void computeCurrentState();
  void sensitivityTo(const std::string& phys_token, xfem::xTensor4& sensitivity);  
};

} // end of namespace

#endif





