/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include <cstdio>
#include "xpNewtonianIncompFluid.h"

//using namespace std;

namespace xfem
{
xpNewtonianIncompFluid::xpNewtonianIncompFluid() 
{
  variables_signature.register_tensor2("stress_viscous");
  variables_signature.register_tensor2("stress_total");
  variables_signature.register_tensor2("strain_rate");
  variables_signature.register_scalar("pressure");
  variables_signature.register_scalar("divergence_V");
  variables_signature.register_vector("velocity");
  variables_signature.register_vector("volume_force");
  
  properties_signature.register_string("MATERIAL_CLASS");
  properties_signature.register_string("NAME");

  //fluid constants
  properties_signature.register_vector("VOLUMIC_FORCE");
  properties_signature.register_scalar("VISCOSITY");
 
  properties.setSignature(&properties_signature);
  properties.astring("MATERIAL_CLASS") = "MATERIAL_NEWTONIAN_FLUID";
}

 void xpNewtonianIncompFluid::checkProperties() 
{
  const bool debug = true;
  double eta = properties.scalar("VISCOSITY");
  assert(eta > 0);
  viscosity = eta;
  if (debug) cout << "viscosite                 : " << eta << endl;
  
  double m1_3 = -1./3.;
  for (int i = 0; i < 3; ++i)
  { 
    for (int j = 0; j < 3; ++j)
    {
      ident2(i,j) = delta(i,j);
      for (int k = 0; k < 3; ++k)
      {
		for (int l = 0; l < 3; ++l)
		{
	  		T4const(i,j,k,l) = delta(i,k) * delta(j,l);
		}
      }
    }
  }
}

//protected member functions
 void  xpNewtonianIncompFluid::sensitivityTo(const std::string& phys_token, 
			       xTensor4& sensitivity) {
  if (phys_token == "strain_rate") {
      
    for (int i = 0; i < 3; ++i)
    { 
      for (int j = 0; j < 3; ++j)
      {
      	for (int k = 0; k < 3; ++k)
      	{
      		for (int l = 0; l < 3; ++l)
      		{
	    	/* NB : cette expression s'applique que l'on prenne D complet ou 
	       	D deviatoire dans la loi de comportement.
	       	Si on veut prendre la loi de comportement sous sa forme deviatoire,
	       	il suffit de modifier T4const en prenant :
	       	- Identite4 -1/3(Identite2(x)Identite2 
	       	- au lieu de Identite4 */
	    	sensitivity(i,j,k,l) = 2*viscosity* T4const(i,j,k,l);
	  		}
		}
      }
    } 
  }
  else { fprintf(stderr, "No sensitivity coded\n"); assert(0); } 
  return; 
}

 void xpNewtonianIncompFluid::computeCurrentState()
{
  const bool mydbg=false;

  // constantes et variables utiles
  xTensor2 D     = curr->tensor2("strain_rate");
  double p       = curr->scalar("pressure");

  // divergence of the velocity field
  curr->scalar("divergence_V") = D.trace() ; 

  // contrainte visqueuse
  curr->tensor2("stress_viscous") = D * 2 * viscosity; 
  
// contrainte totale
  curr->tensor2("stress_total") = curr->tensor2("stress_viscous") - ident2 * p; 
  
  curr->vector("volume_force") = properties.vector("VOLUMIC_FORCE");

  return;
}
}
