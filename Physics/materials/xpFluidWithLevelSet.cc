/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include <cstdio>
#include "xpFluidWithLevelSet.h"

using namespace std;

namespace xfem
{
xpFluidWithLevelSet::xpFluidWithLevelSet() 
{
  variables_signature.register_scalar("levelset");
  variables_signature.register_vector("levelset_gradient");
  variables_signature.register_vector("modified_velocity");
  variables_signature.register_scalar("levelset_timederiv");
  variables_signature.register_scalar("element_size");
  variables_signature.register_scalar("characteristic_time");

  properties.setSignature(&properties_signature);
  properties.astring("MATERIAL_CLASS") = "MATERIAL_NEWTONIAN_FLUID_WITH_LS";
}


void xpFluidWithLevelSet::computeCurrentState(string phys_to_update)
{
	if (phys_to_update=="characteristic_time")
	{
		// constantes et variables utiles
		xVector V = curr->vector("velocity");
		double normV = V.mag();
		if (normV<1e-10) normV=1e-10;
		// in the case where there is no diffusion, ts = h/(2*V.norm) because peclet = inf.
		curr->scalar("characteristic_time")= curr->scalar("element_size")/(2*normV);//*(1/tanh(peclet)-1/peclet);
	}
	else
	{
		xpNewtonianIncompFluid::computeCurrentState();	
		curr->vector("modified_velocity") = curr->vector("levelset_gradient")*
								(curr->vector("velocity")*curr->vector("levelset_gradient"));
	}
}


}
