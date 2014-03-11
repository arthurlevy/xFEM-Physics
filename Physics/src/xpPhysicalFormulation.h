/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
/*! \mainpage Physics : environment to solve physical and multiphysical problems with the Xfem library.
* This module Physics besides the classical Xfem dependencies needs Xext.
*
* A physical resolution needs a formulation which inheritates of xpPhysicalFormulation and a material which inheritates of xMaterial
*
* The formulation contains the fields and mainly allow to assemble a residual thanks to setFunction method and the stiffness matrix thanks to setJacobian method.
*
* Therefore it can be used as an argument in the solve method of a nonLinearSystemSolver (today only exists nonLinearSystemSolverNR).
*
* The material defines the behaviour. One may use a material among the inheritent of xMaterial define in ./material sub-directory which is agreement with the xpPhysicalFormulation used.
*
* In a multiphysical framework one has to define several formulations and code one and only material that may inheritate of two different ones (please refer to ./test/thermoelasticity/ThermoElasticMaterial.h for an example). The different formulations must then be solved using an xpResolutionScheme. Today, the only "multiphysical solver" available is xpLoopSolver which successively solve each formulation untill all the residual are converged.
*
* \see For more info, please refer to my PhD thesis http://tel.archives-ouvertes.fr/tel-00464071
*\author Arthur LEVY

*\date 7th of July 2010
*/


/*! \class xMaterial
 *  \brief defined in Xfem library
*/
/*! \class xElastic
 *  \brief defined in Xfem library
 */
 /*! \include ThermoElasticMaterial
 */
 

#ifndef _MULTPHYSSINGLE_H
#define _MULTPHYSSINGLE_H

#include <string>
#include <cassert>
#include "xAlgorithm.h"
#include "xDoubleManager.h"
#include "xVariabManager.h"
#include "xVector.h"
#include "xEval.h"
#include "xData.h"
#include "xField.h"
#include "xRegion.h"
#include "xPhysSurf.h"
#include "xAssembler.h"
#include "xRenumbering.h"
#include "xSolVisitor.h"

#include "ValueOldAndCurrent.h"

#include "xpTimePilot.h"

namespace lalg
{
    class xCSRMatrix;
    class xCSRVector;
}

using namespace xfem;

/// This class is general framework for defining a physical problem. Its main purpose is to store a field (in the inherited class), and to be able to return a residue (setFunction) and a Jacobian (setjacobian) for applying a Newton Raphson method.
class xpPhysicalFormulation {
public :

  xpPhysicalFormulation(xData * d);
  
  virtual ~xpPhysicalFormulation() {;}
  
  void initialize(xEntityFilter  filter=xAcceptAll());
  
  template <class TYPE>
  void initialize(xEval<TYPE> &eval_initial, xEntityFilter  filter=xAcceptAll())
  { 
    declareApproximation();
    initializeFields(eval_initial);
    TreatmentOfEssEnv(filter, true);
    declareInternalVariables();     
    if(renumbered) pair<int, int> res_renum = renumber();
    updateInternalVariables();
    return;
  }
  
  void initialize(xpTimePilot *external_tp, xEntityFilter  filter=xAcceptAll())
  {setTimePilot(external_tp); initialize(filter);}
  void initialize( xVariabManager *external_VM, xEntityFilter  filter=xAcceptAll())
  {setVarManager(external_VM); initialize(filter);}
  void initialize(xpTimePilot *external_tp, xVariabManager *external_VM, xEntityFilter  filter=xAcceptAll())
  {setTimePilot(external_tp); setVarManager(external_VM); initialize(filter);}
  
  template <class TYPE>
  void initialize(xpTimePilot *external_tp, xEval<TYPE> &eval_initial, xEntityFilter  filter=xAcceptAll())
  {setTimePilot(external_tp); initialize(eval_initial, filter);}
 
  
  /// Dirichlet BC eventually on a restricted region (be sure to clear state with  resetDirichletDofs if needed)
  virtual void TreatmentOfEssEnv(xEntityFilter  filter=xAcceptAll(), bool firstDeclaration=false)
  {     error_in("TreatmentOfEssEnv");     return;   }
  /// just calls the previous one with xAcceptAll
  void TreatmentOfEssEnv(bool firstDeclaration) { TreatmentOfEssEnv(xAcceptAll(), firstDeclaration);}

  ///Neumann BC
  virtual void TreatmentOfNatEnv(xAssembler& ass)
  {     error_in("TreatmentOfNatEnv");     return;   }

  virtual void declareApproximation()
  {     error_in("declareApproximation");     return;   }

  virtual void declareUsefulDofs()
  {     error_in("declareUsefulDofs");     return;   }

  ///unallocates the fixed keys of the essential boundary condition.
  virtual void resetDirichletDofs(xEntityFilter filter=xAcceptAll())
  {     error_in("resetDirichletDofs");     return;   }
  
  /// gestion variables internes
  virtual void declareInternalVariables()
  {     error_in("declareInternalVariables");     return;   }

  /// gestion variables internes
  virtual void updateInternalVariables()
  {     error_in("updateInternalVariables");     return;   }
 
  /// updates the field.
  /// Handles eventual use of renumbering !!
  template <class VECTOR_TYPE >
  void updateFields(VECTOR_TYPE Solu, const double coeff)
{  
  if (renumbered)
    {
      Visit(xAddIndexedScaledSolutionVisitor<typename VECTOR_TYPE::iterator>(Solu.begin(),coeff), 
	    double_manager.begin("dofs"),
	    double_manager.end("dofs")   );
    }
  else
    {
      Visit(xAddScaledSolutionVisitor<VECTOR_TYPE>(Solu.begin(),coeff), 
	    double_manager.begin("dofs"),
	    double_manager.end("dofs")   );
    }
    updateInternalVariables();
  return;
}

  /// initialize the field (for initial conditions, for instance) with the evaluator eval_initial. Note that this should be done before the imposition of Essential boundary condition and the reduction of dof, as in initialize(xEval<TYPE> &eval_initial, xEntityFilter  filter=xAcceptAll())
  virtual void initializeFields(xEval<double> &eval_initial)
  {         error_in("InitializeField<DOUBLE>");     return;   }
  
  void  initializeFields(double T)
  { xEvalConstant<double> my_initializer(T); initializeFields(my_initializer);}
  
  virtual void initializeFields(xEval<xVector> &eval_initial)
  {         error_in("InitializeField<VECTOR>");     return;   }

  void  initializeFields(double x, double y, double z)
  { xEvalConstant<xVector> my_initializer(xVector(x,y,z));  initializeFields(my_initializer);}
  
  // building of the discretized system
  ///sets the residue in the vector F
  virtual void setFunction(lalg::xCSRVector& F)
  {     error_in("setFunction");     return;    }

  /// set the Jacobian (stifness matrix) in the matrix J, usefull for the Newton Raphson method
  virtual void setJacobian(lalg::xCSRMatrix&  J) 
  {     error_in("setJacobian");     return;    }

  /// in a transient framework copies current value to old
  virtual void ShiftFieldCurrentToOld()
  {     error_in("ShiftFieldsCurrentToOld");     return;   }

  ///solution export. The integer details defines the degree of data that will be exported. For instance details = 0 will export temperature only. binary=true for a binary export in Gmsh instead of ASCII. sorted will sort the output in the .pos file. usefull for comparing results compued on differents architectures.
  virtual void exportFields(int details ,const std::string& extension, bool binary=false, bool sorted = false)
  {     error_in("exportFields");     return;   }

  /// updates the xPhysSurf domain, using the levelset that defines the interface between 2 materials.
  void  updateDomain(xLevelSet lsmat);
  /// same with three different materials and two interfaces.
  void  updateDomain(xLevelSet lsmat, xLevelSet lstable);
  /// in a contact framework, we have 2 LS, one for the material, one for the support. They have to be treated separately. example in Appli/soudage_us
  void  updateDomain(xLevelSet lsmat, xLevelSet lstable, xLevelSet lscompo);
  
  /// if the levelset(s) have(has) moved, updates the classification accordingly.
  void updateDomain();

//this still need additional work  virtual void  updateDomain(xLevelSet lsmat, xLevelSet lstable, xLevelSet lscompo);

  ///renumbers dof in order to get a matrix with smaller bandwidth.
  virtual pair<int, int> renumber()
  { error_in("renumber");     return pair<int,int>(0,0);   }

  /// returns the number of non zero elements in the tangent matrix. Usefull to handle better memory with xCSRMatrix::ReAllocate()
  virtual int get_nb_nonzero()
  {  error_in("get_nb_nonzero(");     return 0;   }

  /// In a multiphysical framework one and only variable manager handles every formulations. Therefore, it is external to the formulation.
  void setVarManager(xVariabManager *external_vm) 
  {  variab_manager = external_vm;   return;  }

  xVariabManager* getVarManager()
  { return variab_manager;}

  /// In a multiphysical framework one and only time pilot handles every formulations. Therefore, it is external to the formulation.
  void setTimePilot(xpTimePilot *external_tp)
  {    Multphys_time_pilot = external_tp;   return; }

  xpTimePilot getTimePilot()
  {    return *Multphys_time_pilot;}

  void setOrderOfMagnitude(double M)
  {   order_of_magnitude = M; return; }
  
  double getOrderOfMagnitude()
  {   return order_of_magnitude;}

  string get_name()
  {     return name;   }
  
  int get_ndofs() 
  { return double_manager.size("dofs"); }

  xRegion get_iso_zero()
  { return iso_zero; }    
  
  void set_iso_zero(xRegion iso_)
  { iso_zero = iso_;}

  virtual xField* get_field()
  // to avoid warning message transforme 
  // {     error_in("get_field");     return &xField(&double_manager);   }
  // into where return value is set to null
  {     error_in("get_field");     return (xField *)NULL;   }

  xDoubleManager getDoubleManager() { return double_manager; } 
  
protected:
  /// defined as a pointer since different physics will use the same variable manager.
  xVariabManager* variab_manager;
  ///  defined as a pointer since different physics will use the same time pilot.
  xpTimePilot *Multphys_time_pilot;
  /// The residue will be multiplied by this double to get it dimensionless. Different physics will then have comparable residuals (usefull so far in NewtonRaphson and xpLoopSolver).
  double order_of_magnitude;
  /// In a multi-material framework, the constructor of this physSurf will allocate the type of material to each node of the mesh.
  xPhysSurf* domain;
  /// in a multimaterial framework the boundary region iso_zero may be usefull for applying boundary forces by integrating on this region only.
  xRegion all, iso_zero;
  string name;
  xDoubleManager double_manager;
  xData * thedata;
  xEntityFilter filter_integration;
  xIter domain_beg, domain_end;

  ///dimension of the problem (initialized in the constructor using the dimension of the mesh).
  int spaceDim ;

  ///Degree of interpolation
  int interpolation_degree;

  ///If the useful dofs have been declared, we need to reset them before redeclaring. usefull in TreatmentOfEssDof
  bool usefulDofsDeclared;
  
  /// renumbering of the dofs for optimizing the bandwidth has been operated ?
  bool renumbered;


private:
  /// Bug return function in the case where the method of the inherited class is note coded
  void error_in(string  fct) const
  {   std::cerr << fct << " not coded " << std::endl;   assert(0);   return; }}
  ;
 
#endif
