/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
/*
  Implementation d'un modele de type xpPhysicalFormulation
  Formulation elements finis d'un probleme de mécanique des fluides
  Résolution autonome (hors contexte multiphysiques.
  SLC : 07/2007
*/
#include <fstream>
#include "xAlgorithm.h"
#include "xLinearSystemSolverSuperLU.h"
#include "xEnv.h"
#include "xEval.h"
#include "xExportGmsh.h"
#include "xForm.h"
#include "xFemMatrix.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xMaterialSensitivity.h"
#include "xOperators.h"
#include "xPhysSurf.h"
#include "xSimpleGeometry.h"
#include "xSpace.h"
#include "xValue.h"

#include "mPoint.h"
#include "MaterialCommand.h"
#include "NonUniformMaterialSensitivity.h"

#include "xpEval.h"
#include "xpThermic.h"

#include "xApproxFunction.h"
#include "xDomainStringManager.h"
#include "mAttachableDataContainer.h"
#include "xDomainAttachableData.h"

#include "xCSRMatrix.h"
#include "xCSRVector.h"

using namespace xfem;
using namespace lalg;
using namespace AOMD;


// ===============================================================================
/// Constructeurs de la formulation

/// sur un domaine entièrement maillé
xpThermic::xpThermic (xData * d, int order) :
  xpPhysicalFormulation(d),
  T(&double_manager),
  Tlin(&double_manager),
  Tbub(&double_manager),
  intrule_TT(2*order)
 

{ 
  name = "Thermics";
  interpolation_degree = order;
}

/// sur un domaine limité par une level set
xpThermic::xpThermic (xData * d, xLevelSet& lsmat, int order) :
  xpPhysicalFormulation(d),
  T(&double_manager),
  Tlin(&double_manager),
  Tbub(&double_manager),
  intrule_TT(2*order)
//  intrule_TT (3)
{ 
  name = "Thermics";
  domain =  new xPhysSurf(lsmat,  xClassifyOn("matter"), xClassifyOn("air"));
  if (order==10)  // cas d'un enrichissement à l'interface.
    {
      intrule_TT = xIntegrationRulePartition(2);
    }
  interpolation_degree = order;
}

// ===============================================================================
void xpThermic :: declareApproximation() {
  const bool debug=false;
  cout<<"----- DECLARING APPROXIMATION"<<endl;


  switch (interpolation_degree)
    {
    case 0:
       {
	xSpaceLagrange lag_T("TEMPERATURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_ZERO);
	T.insert(xSpaceFiltered(lag_T, filter_integration));
	// auxiliary fields for non zero boundary conditions imposition
	xSpaceLagrange laglin_T("TEMPERATURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_ZERO);
	xSpaceDifference lagbub_T(lag_T, laglin_T);
	Tlin.insert(xSpaceFiltered(laglin_T, filter_integration));
	Tbub.insert(xSpaceFiltered(lagbub_T, filter_integration));
       }
       break;
    case 10:
      {
	xSpaceLagrange lag_T1("TEMPERATURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);

	//enrichment
	xValKeyExtend key_modifier("_INTERFACE");  
	xScalarFunctionDerivDiscXFEM enrichment_function(*domain);
       	xSpaceFiltered::filter_t filter(bind1st(mem_fun(&xPhysSurf::boundary_strict), domain));
	xSpaceXFEM space_full(lag_T1, enrichment_function, key_modifier);
       	xSpaceFiltered enriched(space_full, filter);
	xSpaceComposite lag_T(lag_T1,enriched);
	T.insert(xSpaceFiltered(lag_T, filter_integration));
	
	// auxiliary fields for non zero boundary conditions imposition
	xSpaceLagrange laglin_T("TEMPERATURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);
	xSpaceDifference lagbub_T(lag_T, laglin_T);
	Tlin.insert(xSpaceFiltered(laglin_T, filter_integration));
	Tbub.insert(xSpaceFiltered(lagbub_T, filter_integration));
      }
      break;
    case 2 :
      {
	xSpaceLagrange lag_T("TEMPERATURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_TWO);
	T.insert(xSpaceFiltered(lag_T, filter_integration));
	// auxiliary fields for non zero boundary conditions imposition
	xSpaceLagrange laglin_T("TEMPERATURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);
	xSpaceDifference lagbub_T(lag_T, laglin_T);
	Tlin.insert(xSpaceFiltered(laglin_T, filter_integration));
	Tbub.insert(xSpaceFiltered(lagbub_T, filter_integration));
      }
      break;
    case 1 :
      {
	xSpaceLagrange lag_T("TEMPERATURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);
	T.insert(lag_T);
	
	// auxiliary fields for non zero boundary conditions imposition
	xSpaceLagrange laglin_T("TEMPERATURE", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);
	xSpaceDifference lagbub_T(lag_T, laglin_T);
	Tlin.insert(laglin_T);
	Tbub.insert(lagbub_T);
      }
      break;
    }
  // declaration des valeurs
  xValueCreator<ValueOldAndCurrentDouble_c>  creator;
  DeclareInterpolation(T,  creator, all.begin(), all.end());
  if (debug) double_manager.PrintForDebug("dclT.dbg");
  return;
}


// ===============================================================================
void xpThermic :: declareInternalVariables() {
  const bool debug=false;
  if (debug) cout <<"----- DECLARING INTERNAL VARIABLES"<<endl;
  //CreateTensorsOldAndCurrentValue_c create_tensors;
  xTensorsValueCreator create_tensors;
  //dans le matériaux ne sont stockees que les variables dites utilisateur.
  //C'est le champ de température utile, c'est à dire qu'il dépend du schema d'integration temporel (il dépend de theta).
  // on stock par exemple T=theta* T(n+1) + (1-theta) T(n) comme variable de température
  DeclareMaterialVariablesCommand_c decl_command(create_tensors, *variab_manager);
  ApplyCommandOnIntegrationRule(decl_command, intrule_TT, all.begin(), all.end()); 
  if (debug) variab_manager->PrintForDebug("Var_initialesThermic.dbg" );
  return;
}

// ===============================================================================
void xpThermic :: declareUsefulDofs() {
  const bool debug = false;
  if (debug) cout << "----- DECLARE DOFS STATE" << endl;
  // creation des degres de liberte utiles
  xStateDofCreator<> snh(double_manager, "dofs");
  DeclareState(T, snh, all.begin(), all.end());
  if (debug) double_manager.PrintForDebug("symT.dbg");
}

// ===============================================================================
void xpThermic :: resetDirichletDofs(xEntityFilter filter){
  const bool debug = false;
  
  if (debug) cout << "----- RESET DIRICHLET DOFS" << endl;
  
  for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin(); 
                                   it != thedata->PhysicalEnv->end(); ++it) 
    {

      const xEnv& env = *it;
      string phys = env.Phys; 
      int type = env.Type; 
      int entity = env.Entity;
      if (phys == "TEMPERATURE" && type == FIX ) 
	{
	  xClassRegion bc(thedata->mesh, entity, env.getDimension());
	  xFilteredRegion<xClassIter, xEntityFilter>  whereToClear(bc.begin(), bc.end(), filter);
	  DeleteState(T,  whereToClear.begin(), whereToClear.end()); // delete key asociated to value
	  for (xDoubleManager::vIter itdofs = double_manager.begin("dofs") ; itdofs!= double_manager.end("dofs"); ++itdofs){
	    (*itdofs)->delState();
	  }
	  double_manager.clear_subset("dofs");
	}
    }
  if (debug) double_manager.PrintForDebug("resT.dbg");
}

// ===============================================================================



void xpThermic::initializeFields(xEval<double> &eval_initial){
  xLinearSystemSolverSuperLU<> solver_l2;
  xAssemblerBasic<> assembler_basic_l2;
  OldAndCurrent_c::current();
  L2Projection(T, eval_initial, assembler_basic_l2, intrule_TT, solver_l2, all.begin(), all.end()); 
  ShiftFieldCurrentToOld();
  return;
}

// =========================================================================================
void xpThermic :: TreatmentOfEssEnv(xEntityFilter  filter,bool firstDeclaration ){
  OldAndCurrent_c::current();
  if(!firstDeclaration) resetDirichletDofs();

  for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin(); 
                                   it != thedata->PhysicalEnv->end(); ++it) 
    {
      const xEnv& env = *it;
      string phys = env.Phys; 
      int type = env.Type; 
      int entity = env.Entity;
      double VAL;
      if (phys == "TEMPERATURE" ) 
	{       
	  if (type == FIX) 
	    {
	      xClassRegion bc(thedata->mesh, entity, env.getDimension());
	      xFilteredRegion<xClassIter, xEntityFilter> whereToApply(bc.begin(), bc.end(), filter);
	      if (env.hasAnEvolution())
		{
		  cout << "Imposed temperature evolution : ";
		  double time = (*Multphys_time_pilot)();
		  double VAL = env.getEvolution()(time);
		  cout << "at time " << time << ", Td=" << VAL << endl;
		  DirichletBoundaryCondition (Tlin, phys, whereToApply.begin(), whereToApply.end(), VAL);
		  DirichletBoundaryCondition (Tbub, phys, whereToApply.begin(), whereToApply.end(), 0.0);  
		}
	      else
		{
		  VAL = env.getValue();
		  DirichletBoundaryCondition (Tlin, phys,whereToApply.begin(), whereToApply.end(), VAL);
		  DirichletBoundaryCondition (Tbub, phys,whereToApply.begin(), whereToApply.end(), 0.0); 
		}
	    }	  
	  else assert(1 == 0);      
      }
    } // End loop over the environment info
  declareUsefulDofs();
  return;
}




// ========================================================================================================
void xpThermic :: TreatmentOfNatEnv(xAssembler& assembler)
{
  bool debug =false;
  
  for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin(); it != thedata->PhysicalEnv->end(); ++it) {
    const xEnv& env = *it;
    if (env.Phys == "SURFACIC_HEAT_FLUX") {
      assert(env.Type == FIX);
      xClassRegion bc(thedata->mesh, env.Entity, env.getDimension());
      double VAL;
      if (env.hasAnEvolution())
		{
		  double time = (*Multphys_time_pilot)();
		  VAL = env.getEvolution()(time);
		}
	  else VAL = env.getValue();
      xEvalConstant<double>  flux(VAL);
      if (debug) cout << "boundary condition imposed : FLUX = " << env.getValue() << endl;
      xFormLinearWithLoad<xValOperator<xIdentity<double> >, xEvalConstant<double> > lin(flux); 
      Assemble(lin, assembler, intrule_TT, T, bc.begin(), bc.end(), xUpperAdjacency()); 
    }
  } 
  return;
}

// ========================================================================================================

void xpThermic:: updateInternalVariables() {
  const bool debug=false;
  if (debug) cout <<"----- UPDATING INTERNAL VARIABLES\n";

  // temperature : computed from current value of T
  //   xEvalField< std::_Identity<double> > val_T(T);
  //SetMaterialVariablesVisitor_c<xEvalField<std::_Identity<double> > > visitor_T("temperature",val_T); 
  //VisitMaterialVariablesCommand_c visit_T(visitor_T, * variab_manager);  

  // temparature gradient : computed from current value of T
  xEvalGradField< std::_Identity<xVector> > grad_temp(T);
  SetMaterialVariablesVisitor_c<xEvalGradField< std::_Identity<xVector>  > > 
                               visitor_gT("temperature_gradient",grad_temp); 
  VisitMaterialVariablesCommand_c visit_gT(visitor_gT, * variab_manager);  

  // Other material variables:
  UpdateMaterialVariablesVisitor_c update_visitor;//("thermic");
  VisitMaterialVariablesCommand_c command_update_visitor(update_visitor, * variab_manager);   

  if (debug) cout << "----- APPLYING UPDATE RULES" << endl;
  //  ApplyCommandOnIntegrationRule(visit_T, intrule_TT, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(visit_gT, intrule_TT, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(command_update_visitor, intrule_TT, domain_beg, domain_end);
  return;
}

// ========================================================================================================

void xpThermic :: setFunction( xCSRVector& b ) {

  const bool debug = false;
  if (debug) cout << "----- ASSEMBLING THE EQUATIONS SYSTEM" << endl;
  xAssemblerBasic<> assembler_F(b);
  double COEFF_ASSMB = 1.0; 

  OldAndCurrent_c::current();

  if (debug) cout << "\t SYSTEM FUNCTION : reset equation system\n";
  b.ZeroArray();

  if (debug) cout << "\t SYSTEM FUNCTION : imposed fluxes\n";
  assembler_F.setCoeff(-COEFF_ASSMB);
  TreatmentOfNatEnv(assembler_F);

  if (debug) cout << "\t SECOND MEMBRE : terme diffusif: -Grad(T*).q\n";
  assembler_F.setCoeff(COEFF_ASSMB);
  GetMaterialVariable_c<xVector> flux("heat_flux", * variab_manager);
  xFormLinearWithLoad< xGradOperator< xIdentity<xVector> >, GetMaterialVariable_c<xVector> > 
                     diffusive_lin(flux);
  Assemble(diffusive_lin, assembler_F, intrule_TT, T, domain_beg, domain_end);

  if (debug) cout << "\t SECOND MEMBRE : terme source: \n";
  assembler_F.setCoeff(COEFF_ASSMB);
  GetMaterialVariable_c<double> source("thermal_source", * variab_manager);
  xFormLinearWithLoad< xValOperator <  xIdentity < double> >,
    GetMaterialVariable_c<double> > source_contrib(source);

  Assemble(source_contrib, assembler_F, intrule_TT, T, all.begin(), all.end());
  return;
}

// ===============================================================================
void xpThermic:: setJacobian(xCSRMatrix& M ) {
  const bool debug = false;
  M.ReAllocate(get_ndofs() ,  get_nb_nonzero());
  if (debug) cout << "----- ASSEMBLING THE JACOBIAN MATRIX" << endl;
  xAssemblerBasic<> assmb_M(M);
  double COEFF_ASSMB = 1.0; 
  assmb_M.setCoeff(COEFF_ASSMB);
  
  if (debug) cout << "\t SYSTEM MATRIX   : reset jacobian\n";
  M.SoftZeroMatrix ( );
  
  if (debug) cout << "\t SYSTEM MATRIX   : Diffusion part\n";
  NonUniformMaterialSensitivity_c<xTensor2> fourier("temperature_gradient",* variab_manager);
  xFormBilinearWithLaw<xGradOperator<xIdentity<xVector> >, 
    NonUniformMaterialSensitivity_c<xTensor2>,
    xGradOperator<xIdentity<xVector> > >     bilin_diffusive(fourier);
  
  Assemble(bilin_diffusive, assmb_M, intrule_TT, T, T, domain_beg, domain_end);      
}
// ========================================================================================================

void xpThermic :: exportFields(int details, const std::string& extension, bool binary, bool sorted) {
  const bool debug = false;

  xExportGmshAscii  pexportascii;  
  xExportGmshAsciiSort  pexportasciisorted;  
  xExportGmshBinary pexportbin;
  xExport *ppexport = &pexportascii;
  if (sorted) ppexport = &pexportasciisorted;
  xExport &pexport = *ppexport;
  OldAndCurrent_c::current();


  //ASCI export
  if (!binary)
    {
      pexport.openFile("Thermic_"+ extension);
      
      //to get a better precision with high order interpolation
      //pexport.setNbSplitDefault(interpolation_degree);
      //pexport.setNbSplit(interpolation_degree+1);

      if (debug) cout << "T & grad(T)" << endl;
      // Temperature
      xEvalField<xIdentity<double> > eval_temp(T);
      Export(eval_temp, pexport, "TEMPERATURE", intrule_TT, domain_beg, domain_end);
 
      SmoothMaterialVariable_c < xVector > eval_flux("heat_flux", *variab_manager);   
      Export(eval_flux, pexport, "HEAT FLUX", intrule_TT, domain_beg, domain_end);
  
      
      if (details > 0)
	{
	  //  Temperature gradient
	  xEvalGradField<xIdentity<xVector> > eval_gradtemp(T);
	  Export(eval_gradtemp, pexport, "TEMPERATURE_GRADIENT", intrule_TT, domain_beg, domain_end);
	 
	  // Internal source
	  SmoothMaterialVariable_c<double> sousource("thermal_source", *variab_manager);
	  Export(sousource, pexport, "SOURCE", intrule_TT, domain_beg, domain_end);

	  if (details > 1)
	    {
	      if (debug) cout << "components" << endl;
	      xExtractCompVector extr_x(0);
	      xExtractCompVector extr_y(1);
	      xEvalGradField<xExtractCompVector > eval_gTx(T, extr_x);
	      xEvalGradField<xExtractCompVector > eval_gTy(T, extr_y);
	      
	      Export(eval_gTx, pexport, "TEMPERATURE_GRADIENT_X", intrule_TT, domain_beg, domain_end);	      
	      Export(eval_gTy, pexport, "TEMPERATURE_GRADIENT_Y", intrule_TT, domain_beg, domain_end);
	      
	      // for 3D problem
	      if (spaceDim == 3)
		{
		  xExtractCompVector extr_z(2);
		  xEvalGradField<xExtractCompVector > eval_gTz(T, extr_z);
		  Export(eval_gTz, pexport, "TMEPERATURE_GRADIENT_Z", intrule_TT, domain_beg, domain_end);
		}
	    }
	}
      pexport.closeFile();
    }
  
  //Binary export SAME AS ABOVE
  else 
    {
      pexportbin.openFile("Thermic_"+ extension);

      pexportbin.setNbSplitDefault(interpolation_degree);
      pexportbin.setNbSplit(interpolation_degree);

      xEvalField<xIdentity<double> > eval_temp(T);
      Export(eval_temp, pexportbin, "TEMPERATURE", intrule_TT, domain_beg, domain_end);
      SmoothMaterialVariable_c < xVector > eval_flux("heat_flux", *variab_manager);
      //Export(eval_flux, pexportbin, "HEAT FLUX", intrule_TT, domain_beg, domain_end);

      
      if (details > 0)
	{
	  xEvalGradField<xIdentity<xVector> > eval_gradtemp(T);
	  Export(eval_gradtemp, pexportbin, "TEMPERATURE_GRADIENT", intrule_TT, domain_beg, domain_end);
	  
	  if (details > 1)
	    {
	      if (debug) cout << "components" << endl;
	      xExtractCompVector extr_x(0);
	      xExtractCompVector extr_y(1);
	      xEvalGradField<xExtractCompVector > eval_gTx(T, extr_x);
	      xEvalGradField<xExtractCompVector > eval_gTy(T, extr_y);
	      Export(eval_gTx, pexportbin, "TEMPERATURE_GRADIENT_X", intrule_TT, domain_beg, domain_end);
	      Export(eval_gTy, pexportbin, "TEMPERATURE_GRADIENT_Y", intrule_TT, domain_beg, domain_end);
	  
	      // for 3D problem
	      if (spaceDim == 3)
		{
		  xExtractCompVector extr_z(2);
		  xEvalGradField<xExtractCompVector > eval_gTz(T, extr_z);
		  Export(eval_gTz, pexportbin, "TMEPERATURE_GRADIENT_Z", intrule_TT, domain_beg, domain_end);
		}
	      
	    }
	
	}
      pexportbin.closeFile();
    }

  if (debug) variab_manager->PrintForDebug("Var_Thermic.dbg"+extension );
  
  return;
}
