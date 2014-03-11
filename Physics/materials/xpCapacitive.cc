/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/


#include <cstdio>
#include "xpCapacitive.h"

using namespace std;

namespace xfem
{
  xpCapacitive::xpCapacitive() : xpConductive()
  {
    properties_signature.register_string("MATERIAL_CLASS");
    properties_signature.register_string("NAME");
    properties_signature.register_scalar("THERMIC_CAPACITY_SLOPE");
    properties_signature.register_scalar("REF_THERMIC_CAPACITY");
    properties.setSignature(&properties_signature);
    properties.astring("MATERIAL_CLASS") = "MATERIAL_THERMO_TRANS";

    // conduction variables
    variables_signature.register_scalar("temperature");
    variables_signature.register_scalar("temperature_timederiv");
    variables_signature.register_scalar("thermic_capacity");


    // Propagation
    variables_signature.register_scalar("characteristic_time");
    variables_signature.register_scalar("element_size");
    variables_signature.register_vector("artificial_temperature_gradient");
  }

  //_______________________________________________________________________

  /// Checks properties, and defines the identity matrix and the diffusion tensorial direction (here it is identity which means that the material is isotropic)
   void xpCapacitive::checkProperties() 
  {

    xpConductive::checkProperties();
    const bool debug = true;
    if (debug) cout << "checking" << endl;

  }
 //_______________________________________________________________________
  // Sensitivities

  ///Useful for xpThermicTransient::SetJacobian.
    void  xpCapacitive::sensitivityTo(const std::string& phys_token, 
				       double& sensitivity)
   {
     if (phys_token == "temperature_time_derivate")
       {	 sensitivity = curr->scalar("thermic_capacity");}

     else if (phys_token == "drho_c/dT")
       {
	 sensitivity =  properties.scalar("THERMIC_CAPACITY_SLOPE");// is actualy d(rhoc)/dT
       }
     else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }

     return;
   }

 void  xpCapacitive::sensitivityTo(const std::string& phys_token, 
				       xVector& sensitivity)
   {
     if (phys_token ==  "dk/dT")
       {
	 sensitivity =  curr->vector("temperature_gradient")    *  0;//0 is actually d(conductivity)/dT
       }
     else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
     return;
   }

  void  xpCapacitive::sensitivityTo(const std::string& phys_token, 
				    xTensor2& sensitivity)
  {
    if (phys_token == "supg_term")
      {
	xVector V = curr->vector("velocity");
	sensitivity = tensor_product(V,V) * curr->scalar("characteristic_time");
      }
    else { xpConductive::sensitivityTo(phys_token, sensitivity); }
  }
  

  //_______________________________________________________________________
  
  /// Compute current state is divided in 2 steps depending on which problem needs the update.
  void xpCapacitive::computeCurrentState(string phys_to_update)
  {
    const bool mydbg=false;

    if (phys_to_update=="diffusion")
      {
	if(mydbg) cout << "computing current state - thermic" << endl;
	xpConductive::computeCurrentState();

	curr->scalar("thermic_capacity")=properties.scalar("REF_THERMIC_CAPACITY")+curr->scalar("temperature") * properties.scalar("THERMIC_CAPACITY_SLOPE");
	if(mydbg) cout << "current state computed for thermic-diffusion" << endl;
      }
    ////////////////////////////
    else if (phys_to_update=="convection")
      {
  	if(mydbg) cout << "computing current state - convection" << endl;
	// constantes et variables utiles
      	xVector V = curr->vector("velocity");
	double normV = V.mag();
	if (normV<1e-10) normV=1e-10;

	// in the case where there is no diffusion, ts = h/(2*V.norm) because peclet = inf.
	curr->scalar("characteristic_time")= curr->scalar("element_size")/(2*normV);//*(1/tanh(peclet)-1/peclet);
	curr->vector("artificial_temperature_gradient") =   curr->vector("temperature_gradient")*V*V  * curr->scalar("characteristic_time");
	if(mydbg) cout << "current state computed for convection" << endl;
      }
    
    else
      {
	cout << "there is a problem with the ComputeCurrentState phys_to_update tag" << endl;
	assert(0);
      }
    return;
  }

 
  void xpCapacitive::computeCurrentState()
  {
    const bool mydbg=false;    
    xpCapacitive::computeCurrentState("diffusion");    
    if ( curr->getSignature()->exist_vector("velocity") )
    {    xpCapacitive::computeCurrentState("convection"); }
  }
  



}
