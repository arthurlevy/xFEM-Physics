/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include <cstdio>
#include "xpCarreau.h"
#include "xTensor2.h"
#include "xDebug.h"

//using namespace std;

namespace xfem
{
xpCarreau::xpCarreau() 
{
  variables_signature.register_scalar("equ_strain_rate");
  properties_signature.register_scalar("CARREAU_TIME");
  properties_signature.register_scalar("CARREAU_POWER");
  properties_signature.register_scalar("INFINITE_VISCOSITY");
  properties.setSignature(&properties_signature);
  properties.astring("MATERIAL_CLASS") = "MATERIAL_CARREAU_FLUID"; 
}

void xpCarreau::checkProperties() 
{

  xpNewtonianIncompFluid::checkProperties();
  const bool debug = true;
  double m = properties.scalar("CARREAU_POWER");
  assert(m > 0);
  double reg = properties.scalar("CARREAU_TIME");
  assert(reg >= 0);
  double infvis = properties.scalar("INFINITE_VISCOSITY");

  if (debug) 
    {
      cout << "carreau power             : " << m << endl;
      cout << "carreau time              : " << reg << endl;
      cout << "carreau infinite viscosity: " << infvis << endl;
    }
}
//_____________________________________________________________
 // Sensitivity to strain rate : mechanical stiffness

  ///Useful for xmElastic::setJacobian (strain) or xmFluidMechanics::setJacobian (strain_rate)

void  xpCarreau::sensitivityTo(const std::string& phys_token, 
			       xTensor4& sensitivity) {

   if (phys_token == "strain_rate") {

    double Deq = curr->scalar("equ_strain_rate");
    xTensor2 D = curr->tensor2("strain_rate");
    double m = properties.scalar("CARREAU_POWER"); 
    double lambda = properties.scalar("CARREAU_TIME");
    double eta_inf = properties.scalar("INFINITE_VISCOSITY");

    double fac1 = viscosity * (m-1)/2 * pow(lambda,2) * pow(1 + pow(lambda*Deq,2) , (m-3)/2 );
    double fac2 = viscosity * ( pow(  1 +  pow(Deq*lambda,2.0) , (m-1)/2) + eta_inf ); 
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
	    sensitivity(i,j,k,l) = 8*fac1*D(i,j)*D(k,l)  +  2*(fac2)*T4const(i,j,k,l) ;
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


  void xpCarreau::computeCurrentState()
{
  const bool mydbg=false;
      // constantes et variables utiles
      double m       = properties.scalar("CARREAU_POWER"); 
      double lambda = properties.scalar("CARREAU_TIME");
      xTensor2 D     = curr->tensor2("strain_rate");
      double eta_inf = properties.scalar("INFINITE_VISCOSITY");

      // equivalent strain rate
      double Deq = sqrt (2.* D.contract(D));
      //   if (Deq>1e4*lambda) Deq=1e4*lambda; //seuillage haut. pour maillage fin
      curr->scalar("equ_strain_rate") = Deq;
      double newton_viscosity = viscosity;

      //curr->scalar("viscosity") = properties.scalar("VISCOSITY") *   pow( Deq , m - 1.0 );
      viscosity *= pow( 1 + pow(lambda*Deq, 2) , (m - 1.0)/2 ) + eta_inf;

      xpNewtonianIncompFluid::computeCurrentState();
      viscosity = newton_viscosity;
      return;   
}
}
