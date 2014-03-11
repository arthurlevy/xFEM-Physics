/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#ifndef __XPNEWTONIANIMCOMPFLUID_H
#define __XPNEWTONIANIMCOMPFLUID_H
#include <iostream>
#include <cassert>
#include <string>
#include "AOMDfwd.h"
#include "xTensor2.h"
#include "xTensor3.h"
#include "xTensor4.h"
#include "xTensors.h"
#include "xMaterial.h"

namespace xfem
{
  ///This class defines a Newtonian material which works with an xpFluidMechanics formulation
  /*! The newtonian material follows the law sigma = eta * D, sigma being the extra stress
This class works with an xpFluidMechanics formulation and needs a material file that contains : 
VISCOSITY:eta0
  */
class xpNewtonianIncompFluid : virtual public xMaterial {

public:
  xpNewtonianIncompFluid();
  void checkProperties(); 
   void computeCurrentState();
   void sensitivityTo(const std::string& phys_token, xfem::xTensor4& sensitivity);

protected:
  double delta(int i, int j) {return ((i == j)?1.0:0.0);}
  xfem::xTensor4  T4const;
  xfem::xTensor2  ident2;
  double viscosity;
};

}

#endif





