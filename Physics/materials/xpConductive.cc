/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/


#include "AOMDfwd.h"
#include "xTensor2.h"
#include "xTensor3.h"
#include "xTensor4.h"
#include "xTensors.h"

#include "xpConductive.h"

namespace xfem
{
xpConductive::xpConductive()
  {
    properties_signature.register_string("MATERIAL_CLASS");
    properties_signature.register_string("NAME");
    properties_signature.register_scalar("THERMAL_SOURCE");
    properties_signature.register_scalar("THERMIC_CONDUCTIVITY");
    variables_signature.register_vector("heat_flux");
    variables_signature.register_vector("temperature_gradient");
    variables_signature.register_scalar("thermal_source");
    properties.setSignature(&properties_signature);
    properties.astring("MATERIAL_CLASS") = "MATERIAL_CONDUCTIVE";
  }

void xpConductive::checkProperties() 
  {
    //check the properties
    double conductivity = properties.scalar("THERMIC_CONDUCTIVITY");
    assert(conductivity > 0.);
    
    //Set the derived informations
    SetThermicConductivityIsotropic (conductivity, thermic_conductivity);
  }  

void  xpConductive::sensitivityTo(const std::string& phys_token, 
		      xTensor2& sensitivity) 
  {
    if (phys_token == "temperature_gradient") sensitivity = -thermic_conductivity;
    else {      fprintf(stderr, "No sensitivity coded\n"); assert(0); }
    return;
  }

void xpConductive::computeCurrentState()
  {
    curr->vector("heat_flux") = - (thermic_conductivity * curr->vector("temperature_gradient"));
    curr->scalar("thermal_source") = properties.scalar("THERMAL_SOURCE");
  }
  
void xpConductive::SetThermicConductivityIsotropic(double k, xTensor2& conductivity){
    for (int i = 0; i < 3; ++i){
      for (int j = 0; j < 3; ++j){
	conductivity(i,j) = k * delta(i,j);
      }
    }
    return;
  }

}
