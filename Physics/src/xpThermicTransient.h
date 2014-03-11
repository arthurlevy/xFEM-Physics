/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPTHERMICTRANSIENT_H
#define _XPTHERMICTRANSIENT_H

#include "xpThermic.h"

using namespace xfem;

/*! xpThermicTransient.h
An instance of this class defines a transient thermal problem.
It is solved using a theta method for the time integration and therefore need the xmEval.h "theta" evaluators xEvalThetaField xEvalGradThetaField.
This class can either solve a diffusion problem aor a pure convection problem eventually using SUPG.
This class heritates from xmMultPhysSingle
can simply work with an xThermoTrans material.
*/
/// formulation for a thermal transient problem (diffusion or convection)
class xpThermicTransient : public xpThermic 
{
public :
  xpThermicTransient (xData *data, double theta_, int interpolation_degree =1) :
    xpThermic(data, interpolation_degree), Theta(theta_) , TypeOfProblem("diffusion"), SUPG(true), first_update(true) {}
  xpThermicTransient (xData *data,  xLevelSet& lsmat, double theta_, int interpolation_degree =1):
    xpThermic(data, lsmat, interpolation_degree), Theta(theta_), TypeOfProblem ("diffusion"), SUPG(true), first_update(true) {}

  // ~xmThermic ();
  
  /// passage du champs courant au champs old necessaire au changement de pas de temps:
  void ShiftFieldCurrentToOld();
  void updateInternalVariables();
  template <class VECTTYPE>
  void setFunction(VECTTYPE& F);
  void setFunction(lalg::xCSRVector& F) {setFunction<lalg::xCSRVector>(F);}

  template <class MATTYPE>
  void setJacobian(MATTYPE&  J);
  void setJacobian(lalg::xCSRMatrix&  J) {setJacobian<lalg::xCSRMatrix>(J);}
  
  /// Two types of problems determined by the TypeOfProblem string
  //! useful only in setFunction and setJacobian
  void setToConvection() { TypeOfProblem = "convection"; name = "Thermique - Convection"; }
  //! useful only in setFunction and setJacobian
  void setToDiffusion() {TypeOfProblem = "diffusion" ; name = "Thermique - Diffusion"; }  
  void setSUPG() {SUPG = true;}
  void unsetSUPG() {SUPG = false;}

  /// solution export
  void exportFields(int details, const std::string& extension, bool binary=false, bool sorted = true);


private:
  ///time integration is done using the theta method.  Theta = 1 for a fully implicite method, and Theta = 0 for a fully explicit method
  double Theta; 
  
  /// this tag will help switching between different formulations. useful only in setFunction and setJacobian
  string TypeOfProblem;

  /// will apply SUPG method instead of galerkin. Useful only in the "convection" problem
  bool SUPG;

  /// once the element size is computed (in updateInternalVariables), this bool is set to false so that we don't recompute it at each update.
  bool first_update;
};


 
#endif
