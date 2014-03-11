/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPTHERMIC_H
#define _XPTHERMIC_H

#include "xpPhysicalFormulation.h"
#include "xRenumbering.h"


using namespace xfem;

class xpThermic : public xpPhysicalFormulation
/// this class  defines a thermal steady-state conductive problem with a source term. It works with an xConductiveWithSource material (or subclass material).
{
public :
  /// In the case of one single material
  xpThermic  (xData * d,int order=1);
  /// In the case of 2 materials (matter and air) separated by a levelset
  xpThermic (xData * d, xLevelSet& lsmat, int order=1);

  ///boundary conditions
  void TreatmentOfEssEnv(xEntityFilter  filter=xAcceptAll(),bool firstDeclaration=false);

  void TreatmentOfNatEnv(xAssembler& ass);

  void declareApproximation();
  void declareUsefulDofs();

  ///un_allocates the blocked keys of the essential boundary condition
  void resetDirichletDofs(xEntityFilter filter=xAcceptAll());

  /// gestion variables internes
  void declareInternalVariables();

  void initializeFields(xEval<double> &eval_initial);

  virtual void updateInternalVariables();
  virtual void setFunction(lalg::xCSRVector& F);
  virtual void setJacobian(lalg::xCSRMatrix&  J);

  virtual pair<int, int> renumber()
  {
    pair< int, int > bws;
    bws =  xReverseCutHillMcKeeNumberingBoost(T, all.begin(), all.end());
    std::cout << name << " - Renumbering ( bandwidth " << bws.first << " -> "
	      << bws.second << " )" << endl;
	renumbered=true;
    return bws;
  }

  virtual int get_nb_nonzero()
  {
    int toto=ComputeNnz( T, all.begin(), all.end() ); cout << name << " number of non zeros :" <<toto<< endl; return toto;
  }  

  virtual  void exportFields(int details, const std::string& extension, bool binary=false, bool sorted =false);

  xDoubleManager getDoubleManager() { return double_manager; } 

  xField* get_field() { return &T;}
  
protected:
  xField T;
  /// usefull for boundary condition with a P2 interpolation
  xField Tlin , Tbub; 
  xIntegrationRulePartition intrule_TT;

};
 
#endif
