/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
/*xpElastic.cc
  Implementation d'un modele de type xpPhysicalFormulation
  Formulation elements finis d'un probleme elastique linéaire
  Résolution autonome (hors contexte multiphysiques.
  SLC : 07/2007
*/
#include <fstream>
#include "xAlgorithm.h"
#include "xLinearSystemSolverSuperLU.h"
#include "xEnv.h"
#include "xEval.h"
#include "xExportGmsh.h"
#include "xCSRMatrix.h"
#include "xCSRVector.h"

#include "xForm.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xMaterialSensitivity.h"
#include "xOperators.h"
#include "xPhysSurf.h"
#include "xSimpleGeometry.h"
#include "xSpace.h"
#include "xValue.h"
#include "xApproxFunction.h"

#include "mPoint.h"
#include "MaterialCommand.h"
#include "NonUniformMaterialSensitivity.h"

#include "xpElastic.h"
#include "xpEval.h"


using namespace xfem;
using namespace lalg;
using namespace AOMD;


// ===============================================================================
// Constructeur de la formulation
xpElastic::xpElastic (xData * d,int order) :
  xpPhysicalFormulation(d),
  U(&double_manager),
  Ulin(&double_manager),
  Ubub(&double_manager),
  intrule_uu(2*order)
{
  name = "Elasticity" ;
  interpolation_degree = order;
}

xpElastic::xpElastic (xData * d, xLevelSet& lsmat, int order) :
  xpPhysicalFormulation(d),
  U(&double_manager),
  Ulin(&double_manager),
  Ubub(&double_manager),
  intrule_uu(2*order)
//  intrule_uu(3)
{ 
  name = "Elasticity" ;
  domain = new xPhysSurf(lsmat, xClassifyOn("matter"), xClassifyOn("air"));
  xRegion toto(domain->getMesh_bnd());
  iso_zero= toto;
  interpolation_degree = order;
}


// ===============================================================================
void xpElastic:: declareApproximation() {
  const bool debug= 0;
  if (debug) cout<<"----- DECLARING APPROXIMATION"<<endl;


  switch (interpolation_degree)
    {
    case 1:
      {
	if (debug) cout << "interpolation d'ordre 1" << endl;
	xSpaceLagrange lagux("DISPLACEMENT_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE);
	xSpaceLagrange laguy("DISPLACEMENT_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
	xSpaceComposite  lag_u(lagux, laguy);
	if (spaceDim ==3)
	  { 
	    xSpaceLagrange laguz("DISPLACEMENT_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_ONE);
	    lag_u.insert(laguz);
	  }
	//inserting the approximation in U
	U.insert(lag_u);// filter_domain_element));
       

	//Ulin.insert(xSpaceFiltered(lag_u,filter_integration )) ;// filter_domain_element));
	//	Ubub.insert(xSpaceFiltered(, filter_integration));//filter_domain_element));
      }
      break;
    case 2:
      {
	if (debug) cout << "interpolation d'ordre 2" << endl;
	xSpaceLagrange lagux("DISPLACEMENT_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_TWO);
	xSpaceLagrange laguy("DISPLACEMENT_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_TWO);
	xSpaceComposite  lag_u(lagux, laguy);
	if (spaceDim ==3)
	  { 
	    xSpaceLagrange laguz("DISPLACEMENT_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_TWO);
	    lag_u.insert(laguz);
	  }
	//inserting the approximation in U
	U.insert(lag_u);// filter_domain_element));
	
	// auxiliary fields for non zero boundary conditions imposition
	xSpaceLagrange laglinux("DISPLACEMENT_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE);
	xSpaceLagrange laglinuy("DISPLACEMENT_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
	xSpaceComposite  laglin_U(laglinux, laglinuy);
	if (spaceDim ==3)
	  {
	    xSpaceLagrange laglinuz("DISPLACEMENT_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_ONE);
	    laglin_U.insert(laglinuz);
	  }
	xSpaceDifference lagbub_U(lag_u, laglin_U);

	Ulin.insert(laglin_U) ;// filter_domain_element));
	Ubub.insert(lagbub_U);//filter_domain_element));
      } 
      break;
    }


  // declaration des valeurs
  xValueCreator<ValueOldAndCurrentDouble_c>  creator;
  // xValueCreator<xValueDouble>  creator;
  DeclareInterpolation(U,  creator, all.begin(), all.end());
  if (debug) double_manager.PrintForDebug("dcl.dbg");
  return;
}

// ===============================================================================
void xpElastic :: declareInternalVariables() {
  const bool debug=false;
  if (debug) cout <<"----- DECLARING INTERNAL VARIABLES"<<endl;
  xTensorsValueCreator create_tensors;

  DeclareMaterialVariablesCommand_c decl_command(create_tensors, *variab_manager);
  ApplyCommandOnIntegrationRule(decl_command, intrule_uu, all.begin(), all.end()); 
  if (debug) variab_manager->PrintForDebug("Var_initiales.dbg" );
  return;
}

// ===============================================================================
void xpElastic :: updateInternalVariables() {
  const bool debug=false;
  if (debug) cout <<"----- UPDATING INTERNAL VARIABLES\n";
  //when declared, Oldandcurrent is switched to old, switch it to current:
  OldAndCurrent_c::current();

  // strain : computed from current value of U
  xEvalGradField< xSymmetrize > gradsym_U(U);
  SetMaterialVariablesVisitor_c<xEvalGradField< xSymmetrize> > visitor_E("strain",gradsym_U); 
  VisitMaterialVariablesCommand_c visit_E(visitor_E,* variab_manager);  

  // Other material variables:
  UpdateMaterialVariablesVisitor_c update_visitor("elastic");
  VisitMaterialVariablesCommand_c command_update_visitor(update_visitor,* variab_manager);   
  
  if (debug) cout << "----- APPLYING UPDATE RULES" << endl;
  ApplyCommandOnIntegrationRule(visit_E, intrule_uu, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(command_update_visitor, intrule_uu, domain_beg, domain_end);
  if (debug) cout << "----- Internal variables updated" << endl;
 return;
}

// ===============================================================================

void xpElastic::initializeFields(xEval<xVector> &eval_initial)
{
  xLinearSystemSolverSuperLU<> solver_l2;
  xAssemblerBasic<> assembler_basic_l2;
  OldAndCurrent_c::current();
  L2Projection(U, eval_initial, assembler_basic_l2, intrule_uu, solver_l2, all.begin(), all.end()); 
  return;
}

void xpElastic::initializeFields(double u_initial, double v_initial, double w_initial)
{
  xEvalConstant<xVector> eval_ini(xVector(u_initial, v_initial, w_initial));
  initializeFields(eval_ini);
  return;
}


// ===============================================================================
void xpElastic :: setFunction( xCSRVector& b ) {

  const bool debug = false;
  if (debug) cout << "----- ASSEMBLING THE EQUATIONS SYSTEM" << endl;
  xAssemblerBasic<> assembler_F(b);
  double COEFF_ASSMB = 1.0; 

  OldAndCurrent_c::current();

  if (debug) cout << "\t SYSTEM FUNCTION : reset equation system\n";
  b.ZeroArray();

  if (debug) cout << "\t SYSTEM FUNCTION : imposed tensions\n";
  assembler_F.setCoeff(-COEFF_ASSMB);
  TreatmentOfNatEnv(assembler_F);

  /*if (debug) cout << "\t SYSTEM FUNCTION : volume forces\n";
  assembler_F.setCoeff(COEFF_ASSMB);
  GetMaterialVariable_c < xVector > body_force("acceleration", *variab_manager);
  xFormLinearWithLoad< xValOperator<xIdentity < xVector > >, 
	               GetMaterialVariable_c < xVector >   >     body_force_contrib(body_force);
  Assemble(body_force_contrib, assembler_F, intrule_uu, U, all.begin(), all.end());
  */ 

  if (debug) cout << "\t SYSTEM FUNCTION :";
  assembler_F.setCoeff(COEFF_ASSMB);
  GetMaterialVariable_c<xTensor2> stress("stress",* variab_manager);
  xFormLinearWithLoad< xGradOperator< xIdentity<xTensor2> >, 
                       GetMaterialVariable_c<xTensor2> > stress_contrib(stress);
  Assemble(stress_contrib, assembler_F, intrule_uu, U, domain_beg, domain_end);

  if(debug)
    {
      cout << " ----  IN ELASTIC -----" <<endl;
      b.PrintArray(cout);
    }

  return;
}

// ===============================================================================
void xpElastic :: setJacobian(xCSRMatrix& M ) {
  const bool debug = false;
  M.ReAllocate(get_ndofs() ,  get_nb_nonzero());
  
  if (debug) cout << "----- ASSEMBLING THE JACOBIAN MATRIX" << endl;
  xAssemblerBasic<> assmb_M(M);
  xAssemblerBasicAndTranspose<> assmb_Mtr(M);

  if (debug) cout << "\t SYSTEM MATRIX   : reset jacobian\n";
  M.SoftZeroMatrix ( );

  if (debug) cout << "\t SYSTEM MATRIX   : stiffness \n";
  NonUniformMaterialSensitivity_c<xTensor4> tangent_module("strain",*variab_manager);
  xFormBilinearWithLaw<xGradOperator<xSymmetrize >, 
	                 NonUniformMaterialSensitivity_c<xTensor4>,
	                 xGradOperator<xSymmetrize > >     Kuu(tangent_module);
  Assemble(Kuu, assmb_M, intrule_uu, U, U, domain_beg, domain_end);  

  return;
}

// ===============================================================================
void xpElastic :: declareUsefulDofs() {
  const bool debug = false;
  if (debug) cout << "----- DECLARE DOFS STATE" << endl;
  // creation des degres de liberte utiles
  xStateDofCreator<> snh(double_manager, "dofs");
  DeclareState(U, snh, all.begin(), all.end());
  //  DeclareState(P, snh, all.begin(), all.end());
  if (debug) double_manager.PrintForDebug("sym.dbg");
}


// ===============================================================================
void xpElastic :: resetDirichletDofs(xEntityFilter filter){
  const bool debug = false;
  if (debug) cout << "----- RESET DOFS" << endl;
  
    for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin(); 
                                   it != thedata->PhysicalEnv->end(); ++it) 
    {
      const xEnv& env = *it;
      string phys = env.Phys; 
      int type = env.Type; 
      int entity = env.Entity;
      if ((phys == "DISPLACEMENT_X" || phys == "DISPLACEMENT_Y" || phys == "DISPLACEMENT_Z" ) && type == FIX ) 
	{
	  xClassRegion bc(thedata->mesh, entity, env.getDimension());
	  xFilteredRegion<xClassIter, xEntityFilter>  whereToClear(bc.begin(), bc.end(), filter);
	  DeleteState(U,  whereToClear.begin(), whereToClear.end());
	  for (xDoubleManager::vIter itdofs = double_manager.begin("dofs") ; itdofs!= double_manager.end("dofs"); ++itdofs){
	    (*itdofs)->delState();
	  }
	  double_manager.clear_subset("dofs");
	}
    }
  if (debug) double_manager.PrintForDebug("resU.dbg");
}


// ========================================================================================================
void xpElastic :: exportFields(int details, const std::string& extension, bool binary, bool sorted) {
  const bool debug = false;

  xExportGmshAscii  pexportascii;
  xExportGmshAsciiSort  pexportasciisort;
  xExportGmshBinary pexportbin;
  xExport *ppexport = &pexportascii;
  if (sorted) ppexport =&pexportasciisort;
  xExport & pexport = *ppexport;
  OldAndCurrent_c::current();


  /////ASCII EXPORT///////
  if (!binary)
    {      
      pexport.openFile("Elastic_"  +  extension);
      //to get a better precision with high order interpolation
      //pexport.setNbSplitDefault(interpolation_degree);
      //pexport.setNbSplit(interpolation_degree+1);
      
      // Displacement
      xEvalField<xIdentity<xVector> > eval_disp(U);
      Export(eval_disp, pexport, "DISPLACEMENT", intrule_uu, domain_beg, domain_end);
 
      
      if (details > 0)
	{
	  //Norme
	  xEvalField<xExtractNormVector > eval_unorm(U);
	  Export(eval_unorm, pexport, "Disp_norm", intrule_uu,domain_beg, domain_end );
	  
	  // strain	  
	  xEvalGradField<xSymmetrize > eval_eps(U);
	  Export(eval_eps, pexport, "epsilon", intrule_uu,domain_beg, domain_end);
 
	  //stress 
	  SmoothMaterialVariable_c < xTensor2 > sigma("stress", *variab_manager);
	  Export(sigma, pexport, "sigma", intrule_uu, domain_beg, domain_end );
	}
      pexport.closeFile();
    }
  
  ///////BINARY EXPORT

  else
    {
      pexportbin.openFile("Elastic_"  +  extension);
      //to get a better precision with high order interpolation
      pexportbin.setNbSplitDefault(interpolation_degree);
      pexportbin.setNbSplit(interpolation_degree);
      
      // Displacement
      xEvalField<xIdentity<xVector> > eval_disp(U);
      Export(eval_disp, pexportbin, "DISPLACEMENT", intrule_uu, domain_beg, domain_end);    
      
      
      if (details > 0)
	{

	  
	  xEvalField<xExtractNormVector > eval_unorm(U);
	  Export(eval_unorm, pexportbin, "Disp_norm", intrule_uu,domain_beg, domain_end );
	  
	  // strain	  
	  xEvalGradField<xSymmetrize > eval_eps(U);
	  Export(eval_eps, pexportbin, "epsilon", intrule_uu,domain_beg, domain_end);

	  
	  //stress 
	  SmoothMaterialVariable_c < xTensor2 > sigma("stress", *variab_manager);
	  Export(sigma, pexportbin, "sigma", intrule_uu, domain_beg, domain_end );
	}
      pexportbin.closeFile();
    }
  return;
}

// ========================================================================================================
void xpElastic :: TreatmentOfEssEnv(xEntityFilter  filter,bool firstDeclaration){

  if(!firstDeclaration) resetDirichletDofs();
  for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin(); 
                                   it != thedata->PhysicalEnv->end(); ++it) 
    {
      const xEnv& env = *it;
      string phys = env.Phys; 
      int type = env.Type; 
      int entity = env.Entity;
      double VAL;
      if (phys == "DISPLACEMENT_X" || phys == "DISPLACEMENT_Y" || phys == "DISPLACEMENT_Z" )
	{       
	  if (type == FIX) {
	    xClassRegion bc(thedata->mesh, entity, env.getDimension());
	    xFilteredRegion<xClassIter, xEntityFilter> whereToApply(bc.begin(), bc.end(), filter);
	    if (env.hasAnEvolution())
	      {
		cout << "Imposed velocity evolution : ";
		double time = (*Multphys_time_pilot)();
		double VAL = env.getEvolution()(time);
		cout << "at time " << time << ", Vd=" << VAL << endl;
		if(interpolation_degree>1)
		  {
		    DirichletBoundaryCondition (Ulin, phys, whereToApply.begin(), whereToApply.end(), VAL);
		    DirichletBoundaryCondition (Ubub, phys, whereToApply.begin(), whereToApply.end(), 0.0);
		  }
		else  DirichletBoundaryCondition (U, phys, whereToApply.begin(), whereToApply.end(), VAL);
	      }
	    else
	      {
		VAL = env.getValue();
		if(interpolation_degree>1)
		  {
		    DirichletBoundaryCondition (Ulin, phys, whereToApply.begin(), whereToApply.end(), VAL);
		    DirichletBoundaryCondition (Ubub, phys, whereToApply.begin(), whereToApply.end(), 0.0); 
		  }
		else DirichletBoundaryCondition (U, phys, whereToApply.begin(), whereToApply.end(), VAL);
	      }
	}
	else assert(1 == 0);      
      }
    } // End loop over the environment info
  declareUsefulDofs();
  return;
}

// ========================================================================================================
void xpElastic :: TreatmentOfNatEnv(xAssembler& assembler)
{
  for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin(); it != thedata->PhysicalEnv->end(); ++it) {
    const xEnv& env = *it;
    
    if (env.Phys == "ELAST_TRACTION_X" || env.Phys == "ELAST_TRACTION_Y"  || env.Phys == "ELAST_TRACTION_Z" ) {
      assert(env.Type == FIX);
      xVector val;
       if (env.hasAnEvolution())
	  {
	  	double time = (*Multphys_time_pilot)();
	  	if (env.Phys == "ELAST_TRACTION_X") val(0) =  env.getEvolution()(time);
	  	if (env.Phys == "ELAST_TRACTION_Y") val(1) =  env.getEvolution()(time);
	  	if (env.Phys == "ELAST_TRACTION_Z") val(2) =  env.getEvolution()(time);
	  }
      else
	  {
	  	if (env.Phys == "ELAST_TRACTION_X") val(0) =  env.getValue();
	  	if (env.Phys == "ELAST_TRACTION_Y") val(1) =  env.getValue();
	  	if (env.Phys == "ELAST_TRACTION_Z") val(2) =  env.getValue();
	  }
	        
    //  xIntegrationRulePartitionBoundary intrule_filtered(NatEnvFilter,interpolation_degree);
      xEvalConstant<xVector>  flux(val);
      xFormLinearWithLoad<xValOperator<xIdentity<xVector> >, xEvalConstant<xVector> > lin(flux); 
      xClassRegion bc(thedata->mesh, env.Entity, env.getDimension());
      Assemble(lin, assembler, intrule_uu, U, bc.begin(), bc.end(), xUpperAdjacency()); 
    }
    
  } 
  return;
}
