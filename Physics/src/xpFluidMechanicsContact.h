/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPFLUIDMECHANICCONTACT_H
#define _XPFLUIDMCHANICCONTACT_H

#include "xpFluidMechanics.h"

using namespace xfem;

/// This class enables to treat contact in a fluid mechanic formulation. The contact is treated using a penalization method which consists in applying a force depending  on the eventual penetration on the boundary. Finally only to new methods are rewritten : setFunction and setJacobian.
class xpFluidMechanicsContact : public xpFluidMechanics
{
public :
  /// lstool represents a tool on which the contact will be possible.
  xpFluidMechanicsContact(xData *data,xLevelSet& lsmat,xLevelSet& lstool,  int interpolation_order=1);

  
  /// only adds a residue on the iso_zero of the levelset of the material to take into account the eventual penalization force. This force is -xi*penetration*normal_tool . penetration is the distance between the point X+V*dt and the tool interface.
  void setFunction(lalg::xCSRVector& F);
  /// The new force is non linear and imposSe an additional term on the sensitivity which is -xi*normal_x_normal
  void setJacobian(lalg::xCSRMatrix& J);


  //  void updateContactCoeff();
 
  /// this method optimizes the value of the penalization coefficient in regard to the general Jacobian matrix of the fluid mechanic volumic problem. In an second order framework, the double allowed_penetration is goal penetration searched.
  void setCoeffContact(double allowed_penetration);


protected:
  /// The tool is represented by this levelset
  xLevelSet LSTool;
  /// this is the proportionality coefficient XI between the penetration distance and the force applyed
  double coeff_contact;

  /// equals 1 for a penalization force f = -xi*d*n  and equals 2 for f = -xi * d*d * n
  int order;
  
  xRegion toolsurface;

}
;
#endif

		       
