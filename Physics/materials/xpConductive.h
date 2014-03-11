/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#ifndef __XPCONDUCTIVE_H
#define __XPCONDUCTIVE_H
#include <iostream>
#include <cassert>
#include <string>

#include "xMaterial.h"

namespace xfem
{

  ///This class defines a  material for use with an xmThermic formulation.
  /*! The  xpConductive material follows the Fourier law q=-lambda*grad(T) and allows a thermal source.
This class needs material file that contains : 
THERMAL_SOURCE:s
THERMIC_CONDUCTIVITY : lambda
  */
class xpConductive : virtual public xMaterial {

public:
  xpConductive();
    
  virtual void checkProperties();
  void  sensitivityTo(const std::string& phys_token, 
		      xTensor2& sensitivity);
  
  virtual void computeCurrentState();
  void SetThermicConductivityIsotropic(double k, xTensor2& conductivity);
  
  inline static double delta(int i, int j) {return ((i == j)?1.0:0.0);}
  
protected:
  xTensor2 thermic_conductivity;
};

} 

#endif


