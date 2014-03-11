/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPELASTIC_H
#define _XPELASTIC_H
#include "xpPhysicalFormulation.h"

using namespace xfem;

/*! xpElastic.h
This file declares the class xpElastic
An instance of this class defines a linear elastic problem.
The diplacement fields and stress fields are seeked
This class heritates from xmMultPhysSingle
can simply work with an xElastic material.
*/
/// formulation for a  compressible elastic problem expressed in displacement (classical).
class xpElastic : public xpPhysicalFormulation
{
public :
  ///For a domain containing a single material
  xpElastic (xData *data, int order=1);
  /// For two domains separetd by an interface given by lsmat 
  xpElastic (xData *data, xLevelSet& lsmat, int order=1);
  // Destructor ~xmFluidMechanics ();

  ///Treatment of boundary conditions
  void TreatmentOfEssEnv(xEntityFilter  filter=xAcceptAll(),bool firstDeclaration=false);
  ///Treatment of boundary conditions
  virtual void TreatmentOfNatEnv(xAssembler& ass);

  void declareApproximation();
  void declareUsefulDofs();
  ///un_allocates the blocked keys of the essential boundary condition
  void resetDirichletDofs(xEntityFilter  filter=xAcceptAll());
  /// gestion variables internes
  void declareInternalVariables();
  /// gestion variables internes
  void updateInternalVariables();

  //initial condition imposition
  void initializeFields(xEval<xVector> &eval_initial);
  void initializeFields(double u_initial, double v_initial, double w_initial=0);

  void setFunction(lalg::xCSRVector& F);
  void setJacobian(lalg::xCSRMatrix&  J);

  virtual pair<int, int> renumber()
  {
    pair< int, int > bws;
    bws =  xReverseCutHillMcKeeNumberingBoost(U, all.begin(), all.end());
    std::cout << name << " - Renumbering ( bandwidth " << bws.first << " -> "
	      << bws.second << " )" << endl;
	renumbered = true;
    return bws;
  }

  virtual int get_nb_nonzero()
  {
    return ComputeNnz( U, all.begin(), all.end() );
  }  
  
  /// solution export
  void exportFields(int details,const std::string& extension, bool binary=false, bool sorted =false);

  xField* get_field() {    return &U;   }

protected:
  xField U ;
  /// useful for non zero boundary conditions
  xField Ulin , Ubub;
  xIntegrationRulePartition intrule_uu;

};

#endif
