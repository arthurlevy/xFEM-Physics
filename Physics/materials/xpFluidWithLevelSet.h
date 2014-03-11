/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#ifndef __FLUIDWITHLS_H
#define __FLUIDWITHLS_H

#include "xpNewtonianIncompFluid.h"

namespace xfem
{
  ///This class defines a Newtonian material with a bimaterial handled with Levelset.It works with an xmFluidMechanics formulation and a xmLevelSetConvector formulation. It is usefull for convecting a levelset using SUPG.
  /*! The newtonian material follows the law sigma = eta * D
This class works with an xmFluidMechanics formulation and needs a material file that contains : 
VISCOSITY:eta0
  */

class  xpFluidWithLevelSet: virtual public xpNewtonianIncompFluid {

public:
	xpFluidWithLevelSet();
	virtual void computeCurrentState(string phys_to_update);
};


} // end of namespace

#endif





