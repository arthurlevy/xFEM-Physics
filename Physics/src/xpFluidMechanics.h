/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPFLUIDMECHANICS_H
#define _XPFLUIDMECHANICS_H
#include "xpPhysicalFormulation.h"

using namespace xfem;
/*!This file declares the class xpFluidMechanics
An instance of this class defines a fluid mechanic problem.
The problem is solved using a mixed velocity / pressure formulation.
This class heritates from xmMultPhysSingle.
Can easily work with an xNewtonianIncomp material
*/
/// formulation for a  incompressible fluid mechanic problem with a mixed V/P formulation.
class xpFluidMechanics : public xpPhysicalFormulation {
public :
  ///For a domain containing a single material
  xpFluidMechanics (xData *data, int order=1);
  /// For two domains separetd by an interface given by lsmat 
  xpFluidMechanics (xData *data, xLevelSet& lsmat, int order=1);
  // ~xpFluidMechanics ();

  ///Treatment of boundary conditions
  void TreatmentOfEssEnv(xEntityFilter  filter=xAcceptAll(),bool firstDeclaration=false);
  ///Treatment of boundary conditions
  void TreatmentOfNatEnv(xAssembler& ass);
  void declareApproximation();
  void declareUsefulDofs();

  void resetDirichletDofs(xEntityFilter  filter=xAcceptAll());
  
  /// gestion variables internes
  void declareInternalVariables();
  /// gestion variables internes
  void updateInternalVariables();

  //initial condition imposition
  ///imposes an initial condition on the velocity field current
  void initializeFields(xEval<xVector> &eval_initial);
  /// imposes an initial homogeneous condition on the velocity field current
  void initializeFields(double u_initial, double v_initial, double w_initial=0);

  void setFunction(lalg::xCSRVector& F);
  void setJacobian(lalg::xCSRMatrix&  J);

  ///calls the parent method and export the different levelsets
 // virtual void  updateDomain(xLevelSet lsmat);
  //  virtual void  updateDomain(xLevelSet lsmat, xLevelSet lstable);
  
  virtual int get_nb_nonzero()
  {
    int toto=ComputeNnz( V, all.begin(), all.end() )+ComputeNnz( P, all.begin(), all.end() ); cout << name << " number of non zeros :" <<toto<< endl; return toto;
  }

  virtual pair<int, int> renumber()
  {
    pair< int, int > bws(0,0);
    // bws =  xReverseCutHillMcKeeNumberingBoost(V, all.begin(), all.end());
    //std::cout << name << " - Velocity Renumbering ( bandwidth " << bws.first << " -> "
    //	      << bws.second << " )" << endl;
	cout << "renumbering not coded" << endl;
    return bws;
    
  }

  /// solution export
  virtual   void exportFields(int details, const std::string& extension, bool binary=false, bool sorted = false);

  /// compute the normal velocity to the level-set and stores it in a new levelset
  xLevelSet compute_vnorm(xLevelSet lsmat );

  ///propagates the levelset (whose iso-zero defines a free surface) using the velocity field and the timestep. The Hamilton Jacobi method is used.
  xLevelSet propagateLevelSet(xLevelSet lsmat);
  ///Same as above but in the case where the norm velocity is already known (vit_norm), tolr is the tolerance for the reinitilisation step (back to grad(phi)=1)
  xLevelSet propagateLevelSet(xLevelSet lsmat, xLevelSet veloc_norm, double tolr = 1e-6);


  xField* get_field() { return &V;}
  xDoubleManager getDoubleManager() { return double_manager; } 
 

protected:

  xField V , P;
  ///for non zero boundary conditions
  xField Vlin , Vbub; 
  xIntegrationRulePartition intrule_uu, intrule_up;
};



#endif
