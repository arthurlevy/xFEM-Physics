/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPELASTICCONTACT_H
#define _XPELASTICCONTACT_H

#include "xpElastic.h"

using namespace xfem;

/// This class enables to treat contact in a linear elastic formulation. The contact is treated using a penalization method which consists in applying a force proportional to the eventual penetration on the boundary. Finally only to new methods are rewritten : setFunction and setJacobian.
class xpElasticContact : public xpElastic
{
public :
  /// lstool represents a tool on which the contact will be done.
  xpElasticContact (xData *data,xLevelSet& lsmat,xLevelSet& lstool,  int interpolation_order=1, int penalization_order=2);
  

  /// only adds a residue on the iso_zero of the levelset of the material to take into account the eventual penalization force. This force is -xi*penetration*normal_tool
  void setFunction(lalg::xCSRVector& F);
  /// The new force is non linear and impose an additional term on the sensitivity which is -xi*normal_x_normal
  void setJacobian(lalg::xCSRMatrix& J);
 
  void setCoeffContact(double c);

protected:
  /// The tool is represented by this levelset
  xLevelSet LSTool;
  /// this is the proportionality coefficient xi between the penetration distance and the force applyed
  double coeff_contact;

  /// equals 1 for a penalization force f = -xi*d*n  and equals 2 for f = -xi * d*d * n
  int order; 

}
;
#endif

		       
