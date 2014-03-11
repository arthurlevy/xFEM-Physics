/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPLSCONVECTOR_H
#define _XPLSCONVECTOR_H

#include "xpThermic.h"
//#include "MaterialCommand.h"

using namespace xfem;

/*! xpLevelSetConvector.h
An instance of this class defines a convection problem.
It is solved using a theta method for the time integration and therefore need the xpEval.h "theta" evaluators xEvalThetaField xEvalGradThetaField.
This class is useful for solving a pure convection problem eventually using SUPG.
works with an xpFluidWithLevelSet material.
*/
/// formulation that allows to convect a levelset using an SUPG method (similar to the 'cousin' class formulation of thermic convection)
class xpLevelSetConvector : public xpThermic 
{
public :
  xpLevelSetConvector (xData *data, double theta_ = 1, int interpolation_degree =1) :
    xpThermic(data, interpolation_degree), Theta(theta_), first_update(true)//,   n(&double_manager)

    { }
  
  /// passage du champs courant au champs old necessaire au changement de pas de temps:
  void ShiftFieldCurrentToOld();
  void updateInternalVariables();
  
  template <class VECTTYPE>
  void setFunction(VECTTYPE& F);
  void setFunction(lalg::xCSRVector& F) {setFunction<lalg::xCSRVector>(F);}

  template <class MATTYPE>
  void setJacobian(MATTYPE&  J);
  void setJacobian(lalg::xCSRMatrix&  J) {setJacobian<lalg::xCSRMatrix>(J);}
  
  /// solution export
  void exportFields(int details, const std::string& extension, bool binary=false);

  /// load an external levelset in the class xField (which is called T because the class inheritates of Thermics)
  void loadLevelSet(xLevelSet external_LS);
  /// exports the xField of the class as a levelset.
  xLevelSet getAsLevelSet();

private:

  /// once the element size is computed (in updateInternalVariables), this bool is set to false so that we don't recompute it at each update.
  bool first_update;
  
 ///time integration is done using the theta method.  Theta = 1 for a fully implicite method, and Theta = 0 for a fully explicit method
  double Theta; 
};

 
#endif
