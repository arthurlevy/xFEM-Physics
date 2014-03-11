/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include <cstdio>
#include "xpPowerLawFluid.h"
#include "xTensor2.h"
#include "xDebug.h"

using namespace std;

namespace xfem
{
xpPowerLawFluid::xpPowerLawFluid() 
{
  variables_signature.register_scalar("equ_strain_rate"); 
  properties_signature.register_scalar("POWER_LAW_INDEX");
  properties_signature.register_scalar("REG_COEFF");
  properties.setSignature(&properties_signature);
  properties.astring("MATERIAL_CLASS") = "MATERIAL_POWER_LAW_FLUID";
 }

void xpPowerLawFluid::checkProperties() 
{
  xpNewtonianIncompFluid::checkProperties();

  const bool debug = true;
  double m = properties.scalar("POWER_LAW_INDEX");
  assert(m > 0);
  double reg = properties.scalar("REG_COEFF");
  assert(reg >= 0);

  if (debug) {
    cout << "indice loi puissance      : " << m << endl;
    cout << "facteur de regularisation : " << reg << endl;
  }



}

//_____________________________________________________________
 // Sensitivity to strain rate : mechanical stiffness

  ///Useful for xpFluidMechanics::setJacobian (strain_rate)

void  xpPowerLawFluid::sensitivityTo(const std::string& phys_token, 
			       xTensor4& sensitivity) {


   if (phys_token == "strain_rate") {

    double Deq = curr->scalar("equ_strain_rate");
    xTensor2 D = curr->tensor2("strain_rate");
    double m = properties.scalar("POWER_LAW_INDEX"); 

    double fac1 = 2.* viscosity * pow( Deq , m - 1.0 );// 2*eta0(Deq)^(m-1)
    double fac2 = 2.* (m-1.0) * pow( Deq , - 2.0 ); // 2(m-1)/Deq²
    for (int i = 0; i < 3; ++i) { 
      for (int j = 0; j < 3; ++j) {
	for (int k = 0; k < 3; ++k) {
	  for (int l = 0; l < 3; ++l) {
	    /* NB : cette expression s'applique que l'on prenne D complet ou 
	       D deviatoire dans la loi de comportement.
	       Si on veut prendre la loi de comportement sous sa forme deviatoire,
	       il suffit de modifier T4const en prenant :
	       - Identite4 -1/3(Identite2(x)Identite2 
	       - au lieu de Identite4 */
	    sensitivity(i,j,k,l) = fac1 * (T4const(i,j,k,l) + fac2*D(i,j)*D(k,l) );
	  }
	}
      }
    } 
  }
  else { 
    fprintf(stderr, "No sensitivity coded\n"); assert(0); 
  }
  return; 
}

//_____________________________________________

void xpPowerLawFluid::computeCurrentState()
{
  const bool mydbg=false;
 
      // constantes et variables utiles
      double m       = properties.scalar("POWER_LAW_INDEX"); 
      double eps_reg = properties.scalar("REG_COEFF");
      xTensor2 D     = curr->tensor2("strain_rate");

      // equivalent strain rate
      double Deq2 = 2.* D.contract(D);
      double newton_viscosity = viscosity;

      // regularization : useful especially when D=0
      if (Deq2<eps_reg) { Deq2 =  eps_reg; }
      double Deq = sqrt(Deq2);
      curr->scalar("equ_strain_rate") = Deq;

      //Viscosity
      viscosity *= pow( Deq , m - 1.0 );
      xpNewtonianIncompFluid::computeCurrentState();
      viscosity = newton_viscosity;
      return;   
}

}
