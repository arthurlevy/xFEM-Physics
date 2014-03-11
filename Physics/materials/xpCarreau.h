/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/


#ifndef __XPCARREAU_H
#define __XPCARREAU_H
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
  ///This class defines a carreau material which is a quasi newtonian material.It works with an xmFluidMechanics formulation
  /*! The carreau material follows the law sigma = 2 eta * D where eta = eta0*(1+(lambda*Deq)^2)^((m-1)/2) + eta_inf.
This class works with an xmFluidMechanics formulation and needs a material file that contains : 
CARREAU_TIME:lambda
CARREAU_POWER: m
VISCOSITY:eta0
INFINITE_VISCOSITY = eta_inf
VOLUMIC_FORCES: vector (see xpNewtonianIncompFluid)
  */

class xpCarreau : virtual public xpNewtonianIncompFluid {
public:
  xpCarreau();
  void checkProperties(); 
  void computeCurrentState();
  void sensitivityTo(const std::string& phys_token, xfem::xTensor4& sensitivity);  
};
} // end of namespace
#endif





