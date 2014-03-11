/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
/*
  Implementation d'un modele de type xmMultPhysSingle
  Formulation elements finis d'un probleme de mécanique des fluides
  Résolution autonome (hors contexte multiphysiques.
  SLC : 07/2007
*/
#include <fstream>
#include "xAlgorithm.h"

#include "xLinearSystemSolverSuperLU.h"
#include "xEnv.h"
#include "xEval.h"
#include "xpEval.h"
#include "xExportGmsh.h"
#include "xForm.h"

#include "xCSRMatrix.h"
#include "xCSRVector.h"

#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xMaterialSensitivity.h"
#include "xOperators.h"
#include "xPhysSurf.h"
#include "xSimpleGeometry.h"
#include "xSpace.h"
#include "xValue.h"
#include "xLevelSetOperators.h"

#include "mPoint.h"
#include "MaterialCommand.h"
#include "NonUniformMaterialSensitivity.h"

#include "xpFluidMechanics.h"


using namespace xfem;
using namespace lalg;
using namespace AOMD;


// ===============================================================================
// Constructeurs de la formulation
xpFluidMechanics::xpFluidMechanics (xData * d, int order) :
  xpPhysicalFormulation(d),
  V(&double_manager),
  Vlin(&double_manager),
  Vbub(&double_manager),
  P(&double_manager),
  intrule_uu(2*order),
  intrule_up(2*order)
{
  name = "Fluid Mechanics";
  interpolation_degree = order;
 /* if (order==1)//P1+/P1 case
    {
      //  interpolation_degree = 3;
      //    intrule_uu = xIntegrationRulePartition(order);
      intrule_uu = xIntegrationRulePartition(3);
      intrule_up = intrule_uu;
    }
    */
}

xpFluidMechanics::xpFluidMechanics (xData * d, xLevelSet& lsmat, int order) :
  xpPhysicalFormulation(d),
  V(&double_manager),
  Vlin(&double_manager),
  Vbub(&double_manager),
  P(&double_manager),
  intrule_uu(2*order),
  intrule_up(2*order)
{ 
  name = "Fluid Mechanics";
  domain =  new xPhysSurf(lsmat, xClassifyOn("matter"),  xClassifyOn("air"));
  xExportGmshAscii pexport;
  pexport.openFile("LS_00");
  Export(lsmat, pexport, "LS_00");
  pexport.closeFile();
  xRegion toto(domain->getMesh_bnd());
  iso_zero= toto;
  interpolation_degree = order;
  if ( order==10)//XFEM case
    {
      //   interpolation_degree = 3;
      intrule_uu = xIntegrationRulePartition(2);
      intrule_up = xIntegrationRulePartition(2);
    }
}


// ===============================================================================
void xpFluidMechanics :: declareApproximation() {
  const bool debug=false;
  if (debug) cout<<"----- DECLARING APPROXIMATION"<<endl;

  //choose between a P2/P1 formulation or a P1+/P1
  switch (interpolation_degree)
    {
    case 2:
      {
	xSpaceLagrange lagvx("VELOCITY_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_TWO);
	xSpaceLagrange lagvy("VELOCITY_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_TWO);
	xSpaceComposite  lag_V(lagvx, lagvy);
	if (spaceDim == 3) 
	  {
	    xSpaceLagrange lagvz("VELOCITY_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_TWO);
	    lag_V.insert(lagvz);
	  }
	
	// auxiliary fields for non zero boundary conditions imposition
	xSpaceLagrange laglinvx("VELOCITY_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE);
	xSpaceLagrange laglinvy("VELOCITY_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
	
	xSpaceComposite  laglin_V(laglinvx, laglinvy);
	if (spaceDim == 3) 
	  {
	    xSpaceLagrange laglinvz("VELOCITY_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_ONE);
	    laglin_V.insert( laglinvz);
	  }
	xSpaceDifference lagbub_V(lag_V, laglin_V);

	V.insert(xSpaceFiltered(lag_V, filter_integration));
	Vlin.insert(xSpaceFiltered(laglin_V, filter_integration));
	Vbub.insert(xSpaceFiltered(lagbub_V, filter_integration));
      }
      break;


    case 1:
      {
	xSpaceLagrange laglinx("VELOCITY_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE);
	xSpaceLagrange lagliny("VELOCITY_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
	
	xSpaceLagrange lag3x("VELOCITY_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_THREE);
	xSpaceLagrange lag3y("VELOCITY_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_THREE);
	
	xSpaceFiltered bubx(lag3x, bubbleFunction(spaceDim));
	xSpaceFiltered buby(lag3y, bubbleFunction(spaceDim));
	
	xSpaceComposite  laglin_V(laglinx, lagliny);
	
	xSpaceComposite  lag_V(laglinx,bubx,lagliny,buby);
	if (spaceDim == 3) 
	  {
	    xSpaceLagrange laglinz("VELOCITY_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_ONE);
	    xSpaceLagrange  lag3z("VELOCITY_Z", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_THREE);
	    xSpaceFiltered bubz(lag3z, bubbleFunction(spaceDim));
	    laglin_V.insert(laglinz);
	    lag_V.insert(laglinz);
	    lag_V.insert(bubz);
	  }
	

	xSpaceDifference lagbub(lag_V, laglin_V);	
	Vlin.insert(laglin_V);
	Vbub.insert(lagbub);
       	V.insert(lag_V);
      }
      break;
      
      case 10: //P2X-FEM / P1
      {
	xSpaceLagrange laglinx("VELOCITY_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE);
	xSpaceLagrange lagliny("VELOCITY_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
	
	xSpaceLagrange lag3x("VELOCITY_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_THREE);
	xSpaceLagrange lag3y("VELOCITY_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_THREE);
	
	xSpaceFiltered bubx(lag3x, bubbleFunction(spaceDim));
	xSpaceFiltered buby(lag3y, bubbleFunction(spaceDim));
	
	xSpaceComposite  laglin_V(laglinx, lagliny);
	
	xSpaceComposite  lag_V(laglinx,bubx,lagliny,buby);
	if (spaceDim == 3) 
	  {
	    xSpaceLagrange laglinz("VELOCITY_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_ONE);
	    xSpaceLagrange  lag3z("VELOCITY_Z", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_THREE);
	    xSpaceFiltered bubz(lag3z, bubbleFunction(spaceDim));
	    laglin_V.insert(laglinz);
	    lag_V.insert(laglinz);
	    lag_V.insert(bubz);
	  }

	
	xValKeyExtend key_modifier("_INCLUSIOND");  
	xScalarFunctionDerivDiscXFEM enrichment_function(*domain);
       	xSpaceFiltered::filter_t filter(bind1st(mem_fun(&xPhysSurf::boundary_strict), domain));
	xSpaceXFEM space_full(lag_V, enrichment_function, key_modifier);
       	xSpaceFiltered enriched(space_full, filter);
	lag_V.insert(enriched);
	

	xSpaceDifference lagbub(lag_V, laglin_V);	
	Vlin.insert(laglin_V);
	Vbub.insert(lagbub);
       	V.insert(lag_V);
      }
      break;

      
    case 0: //Never use this horrible P0/P1 interpolation!
      {
	xSpaceLagrange lagvx("VELOCITY_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE);
	xSpaceLagrange lagvy("VELOCITY_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
	
	xSpaceComposite  lag_V(lagvx, lagvy);
	
	if (spaceDim == 3) 
	  {
	    xSpaceLagrange lagvz("VELOCITY_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_TWO);   	
	    lag_V.insert(lagvz);
	  }
	Vlin.insert(xSpaceFiltered(lag_V, filter_integration));
	V.insert(xSpaceFiltered(lag_V, filter_integration));
      }
      break;
    }

  xSpaceLagrange lagpr("PRESSURE",   xSpace::SCALAR,   xSpaceLagrange::DEGREE_ONE);
  P.insert(xSpaceFiltered(lagpr, filter_integration));
  
  // declaration des valeurs
  xValueCreator<ValueOldAndCurrentDouble_c>  creator;
  DeclareInterpolation(V,  creator, all.begin(), all.end());
  DeclareInterpolation(P,  creator, all.begin(), all.end());

  if (debug) double_manager.PrintForDebug("dcl.dbg");
  return;
}


// ===============================================================================
void xpFluidMechanics :: declareInternalVariables() {
  const bool debug=false;
  if (debug) cout <<"----- DECLARING INTERNAL VARIABLES"<<endl;
  //CreateTensorsOldAndCurrentValue_c create_tensors;
  xTensorsValueCreator create_tensors;
 
  DeclareMaterialVariablesCommand_c decl_command(create_tensors, * variab_manager);
  ApplyCommandOnIntegrationRule(decl_command, intrule_uu, all.begin(), all.end()); 
  if (debug) variab_manager->PrintForDebug("Var_initiales_fluid.dbg" );
  return;
}



// ===============================================================================
void xpFluidMechanics :: updateInternalVariables() {
  const bool debug=false;
  if (debug) cout <<"----- UPDATING INTERNAL VARIABLES\n";

  //when declared, Oldandcurrent is switched to old, switch it to current:
  OldAndCurrent_c::current();

  // Velocity 
  xEvalField< std::_Identity<xVector> > val_V(V);
  SetMaterialVariablesVisitor_c<xEvalField< std::_Identity <xVector> > > visitor_V("velocity",val_V); 
  VisitMaterialVariablesCommand_c visit_V(visitor_V, * variab_manager);  

  // strain rate : computed from current value of V
  xEvalGradField< xSymmetrize > gradsym_V(V);
  SetMaterialVariablesVisitor_c<xEvalGradField< xSymmetrize> > visitor_D("strain_rate",gradsym_V); 
  VisitMaterialVariablesCommand_c visit_D(visitor_D, * variab_manager);  

  // pressure : computed from current value of P
  xEvalField< std::_Identity<double> > val_P(P);
  SetMaterialVariablesVisitor_c<xEvalField<std::_Identity<double> > > visitor_P("pressure",val_P); 
  VisitMaterialVariablesCommand_c visit_P(visitor_P, * variab_manager);  


  // Other material variables:
  // deviatoric strain rate, equivalent strain rate, deviatoric stress, total stress
  UpdateMaterialVariablesVisitor_c update_visitor("fluidmechanics");
  VisitMaterialVariablesCommand_c command_update_visitor(update_visitor, * variab_manager);   
  
  if (debug) cout << "----- APPLYING UPDATE RULES" << endl;
  ApplyCommandOnIntegrationRule(visit_V, intrule_uu, all.begin(), all.end());
  ApplyCommandOnIntegrationRule(visit_D, intrule_uu, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(visit_P, intrule_uu, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(command_update_visitor, intrule_uu, domain_beg, domain_end);
  if (debug) variab_manager->PrintForDebug("Var_updatees_fluid.dbg" );
  return;
}


// ===============================================================================


//Velocity InitialCondition imposition
void xpFluidMechanics::initializeFields(xEval<xVector> &eval_initial)
{
  xLinearSystemSolverSuperLU<> solver_l2;
  xAssemblerBasic<> assembler_basic_l2;
  OldAndCurrent_c::current();
  L2Projection(V, eval_initial, assembler_basic_l2, intrule_uu, solver_l2, all.begin(), all.end()); 
  return;
}

void xpFluidMechanics::initializeFields(double u_initial, double v_initial, double w_initial)
{
  xEvalConstant<xVector> eval_ini(xVector(u_initial, v_initial, w_initial));
  initializeFields(eval_ini);
  return;
}



// ===============================================================================
void xpFluidMechanics :: setFunction( xCSRVector& b ) {

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

   if (debug) cout << "\t SYSTEM FUNCTION : volume forces\n";
  assembler_F.setCoeff(COEFF_ASSMB);
  //xUniformMaterialParameter<xVector> body_force("VOLUMIC_FORCE");
  //////
  /*
  GetMaterialVariable_c<xVector> body_force("volume_force", * variab_manager);
  xFormLinearWithLoad< xValOperator<xIdentity<xVector> >, 
	              GetMaterialVariable_c<xVector>   >     body_force_contrib(body_force);
  Assemble(body_force_contrib, assembler_F, intrule_uu, V, domain_beg, domain_end);
  */

  if (debug) cout << "\t SYSTEM FUNCTION : Rheological part - total stress approach: Sig = Sig_visc - p*Id2\n";
  assembler_F.setCoeff(COEFF_ASSMB);
  GetMaterialVariable_c<xTensor2> stress_full("stress_total", * variab_manager);
  xFormLinearWithLoad< xGradOperator< xIdentity<xTensor2> >, 
                       GetMaterialVariable_c<xTensor2> > stress_contrib(stress_full);
  Assemble(stress_contrib, assembler_F, intrule_uu, V, domain_beg, domain_end);

  if (debug) cout << "\t SYSTEM FUCNTION : Incompressibility part\n";
  assembler_F.setCoeff(-COEFF_ASSMB); // minus!!
  
  GetMaterialVariable_c<double> div_V("divergence_V", * variab_manager);
  xFormLinearWithLoad< xValOperator<xIdentity<double> >, 
    GetMaterialVariable_c<double>   >     pdiv_contrib(div_V);
    
  Assemble(pdiv_contrib, assembler_F, intrule_uu, P, domain_beg, domain_end);

  return;
}

// ===============================================================================
void xpFluidMechanics :: setJacobian(xCSRMatrix& M ) {
  const bool debug = false;
  if (debug) cout << "----- ASSEMBLING THE JACOBIAN MATRIX" << endl;
  xAssemblerBasic<> assmb_M(M);
  xAssemblerBasicAndTranspose<> assmb_Mtr(M);

  if (debug) cout << "\t SYSTEM MATRIX   : reset jacobian\n";
  M.SoftZeroMatrix ( );

  if (debug) cout << "\t SYSTEM MATRIX   : Rheological part\n";
  NonUniformMaterialSensitivity_c<xTensor4> tangent_module("strain_rate",* variab_manager);
  xFormBilinearWithLaw<xGradOperator<xSymmetrize >, 
	                 NonUniformMaterialSensitivity_c<xTensor4>,
	                 xGradOperator<xSymmetrize > >     Kvv(tangent_module);
  Assemble(Kvv, assmb_M, intrule_uu, V, V, domain_beg, domain_end);  
  if (debug)  {    
    std::ofstream out1("matK1.m"); 
    out1 << "after Rheo term.\n";
    M.OutputMatrixOctaveFormat(out1);
  }

  if (debug) cout << "\t SYSTEM MATRIX   : Pressure and Incompressibility parts\n";
  xFormBilinearWithoutLaw<xGradOperator< xTrace >, 
                          xValOperator<xIdentity<double> > > Kvp; 
  assmb_M.setCoeff(-1);
  Assemble(Kvp, assmb_M, intrule_up, V, P, domain_beg, domain_end);
  xFormBilinearWithoutLaw< xValOperator<xIdentity<double> >,
                         xGradOperator< xTrace >   > Kpv; 
  Assemble(Kpv, assmb_M, intrule_up, P, V, domain_beg, domain_end);
  ///

  if (debug)     {   
    std::ofstream out2("matK2.m"); 
    out2 << "after Pdiv terms.\n"; 
    M.OutputMatrixOctaveFormat(out2);
  }
  return;
}

// ========================================================================================================
void xpFluidMechanics :: declareUsefulDofs() {
  const bool debug = false;
  if (debug) cout << "----- DECLARE DOFS STATE" << endl;
  // creation des degres de liberte utiles
  xStateDofCreator<> snh(double_manager, "dofs");
  DeclareState(V, snh, all.begin(), all.end());
  DeclareState(P, snh, all.begin(), all.end());
  if (debug) double_manager.PrintForDebug("sym.dbg");
}

// ========================================================================================================
void xpFluidMechanics :: resetDirichletDofs(xEntityFilter filter){
  const bool debug = false;
  if (debug) cout << "----- RESET DOFS" << endl;
  for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin(); 
                                   it != thedata->PhysicalEnv->end(); ++it) 
    {
      const xEnv& env = *it;
      string phys = env.Phys; 
      int type = env.Type; 
      int entity = env.Entity;
      if ((phys == "VELOCITY_X" || phys == "VELOCITY_Y" || phys == "VELOCITY_Z" )&&type == FIX ) 
	{
	  xClassRegion bc(thedata->mesh, entity, env.getDimension());
	  xFilteredRegion<xClassIter, xEntityFilter>  whereToClear(bc.begin(), bc.end(), filter);
	  DeleteState(V,  whereToClear.begin(), whereToClear.end());
	  for (xDoubleManager::vIter itdofs = double_manager.begin("dofs") ; itdofs!= double_manager.end("dofs"); ++itdofs){
	    (*itdofs)->delState();
	  }
	  double_manager.clear_subset("dofs");
	}
    }
  if (debug) double_manager.PrintForDebug("resV.dbg");
}

// ========================================================================================================
void xpFluidMechanics :: exportFields(int details, const std::string& extension, bool binary, bool sorted) {
  const bool debug = false;

  xExportGmshAscii  pexportascii;
  xExportGmshBinary  pexportbin; 
  xExportGmshAsciiSort pexportasciisort;
  xExport *ppexport = &pexportascii;
  if (sorted) ppexport = &pexportasciisort;
  xExport &pexport =  *ppexport;
  
  if (debug)      double_manager.PrintForDebug("V"+extension+".dbg");


  ///ASCII EXPORT///
    if (!binary) 
    {
      pexport.openFile("Fluid_" + extension);
      // to get a better precision with high order interpolation
      // pexport.setNbSplitDefault(5);
      //  pexport.setNbSplit(5);
      xEvalField<xIdentity<xVector> > eval_velo(V);
      Export(eval_velo, pexport, "VELOCITY", intrule_uu, domain_beg, domain_end);

      
      if (details > 0)
	{
	  xEvalField<xIdentity<double> > eval_press(P); 
	  Export(eval_press, pexport, "PRESSURE", intrule_up, domain_beg, domain_end);
	  if (debug) cout << "stress export" << endl;
	  SmoothMaterialVariable_c < xTensor2 > sigmavis("stress_viscous", *variab_manager);
	  Export(sigmavis, pexport, "Extra_Stress", intrule_uu, domain_beg, domain_end );
	   
	  xEvalGradField<xSymmetrize> eval_D(V);
	  Export(eval_D, pexport, "D", intrule_uu, domain_beg, domain_end );
	  if (debug) cout << "div v" << endl;
	  // divergence du champ de vitesses
	  xEvalGradField<xTrace > eval_divV(V);
	  Export(eval_divV, pexport, "Div_V", intrule_uu, domain_beg, domain_end);	
	  
	  if(details > 1)
	    {

	      xEvalField<xExtractNormVector > eval_vnorm(V);
	      Export(eval_vnorm, pexport, "V_norm", intrule_uu, domain_beg, domain_end);
	      
	      SmoothMaterialVariable_c < xTensor2 > sigmatot("stress_total", *variab_manager);
	      Export(sigmatot, pexport, "Total_Stress", intrule_uu, domain_beg, domain_end );     
	      
	      if (debug) cout << "components" << endl;
	      // Velocity components & velocity norm
	      xExtractCompVector extr_x(0);  xExtractCompVector extr_y(1); xExtractCompVector extr_z(2);
	      xEvalField<xExtractCompVector > eval_vlx(V, extr_x);
	      xEvalField<xExtractCompVector > eval_vly(V, extr_y);
	      xEvalField<xExtractCompVector > eval_vlz(V, extr_z);
	      Export(eval_vlx, pexport, "V_X", intrule_uu, domain_beg, domain_end);
	      Export(eval_vly, pexport, "V_Y", intrule_uu, domain_beg, domain_end);
	      if (spaceDim == 3) Export(eval_vlz, pexport, "V_Z", intrule_uu, domain_beg, domain_end);
	      
	      if (debug) cout << "strain rate" << endl;
	      // taux de deformation
	      xExtractCompTensor extr_xx(0,0);  xExtractCompTensor extr_xy(0,1);  xExtractCompTensor extr_yy(1,1);
	      xExtractCompTensor extr_zz(2,2);  xExtractCompTensor extr_xz(0,2);  xExtractCompTensor extr_yz(1,2);
	      
	      xEvalGradField<xExtractCompTensor> eval_Dxx(V, extr_xx);
	      xEvalGradField<xExtractCompTensor > eval_Dyy(V, extr_yy);
	      xEvalGradField<xExtractCompTensor > eval_Dxy(V, extr_xy);
	      xEvalGradField<xExtractCompTensor > eval_Dzz(V, extr_zz);
	      xEvalGradField<xExtractCompTensor > eval_Dxz(V, extr_xz);
	      xEvalGradField<xExtractCompTensor > eval_Dyz(V, extr_yz);
	      
	      Export(eval_Dxx, pexport, "D_XX", intrule_uu, domain_beg, domain_end);
	      Export(eval_Dyy, pexport, "D_YY", intrule_uu, domain_beg, domain_end);
	      if(spaceDim == 3) Export(eval_Dzz, pexport, "D_ZZ", intrule_uu, domain_beg, domain_end);
	      Export(eval_Dxy, pexport, "D_XY", intrule_uu, domain_beg, domain_end);
	      if(spaceDim == 3) 
		{
		  Export(eval_Dxz, pexport, "D_XZ", intrule_uu, domain_beg, domain_end);
		  Export(eval_Dyz, pexport, "D_YZ", intrule_uu, domain_beg, domain_end);
		}
   
	      if (debug) cout << "grad_p" << endl;
	      xEvalGradField<xIdentity<xVector> > eval_grad_press(P);
	      Export(eval_grad_press, pexport, "Pressure_Grad", intrule_up, domain_beg, domain_end);  
	      
	    
	      
	      SmoothMaterialVariable_c < double > div_v("divergence_V", *variab_manager);
	      Export(div_v, pexport, "div_mat", intrule_uu, domain_beg, domain_end );    
	      
	      if(details > 4)
		{
		  SmoothMaterialVariable_c < double > visco("viscosity", *variab_manager);
		  Export(visco, pexport, "viscosity", intrule_uu, domain_beg, domain_end );  
		}
	    }
 	}
      
      pexport.closeFile();
    }
  
  ///////BINARY EXPORT
  else
    {
      pexport.openFile("Fluid_" + extension);
      //to get a better precision with high order interpolation
      pexport.setNbSplitDefault(5);
      pexport.setNbSplit(5);
      
      if (debug) cout << "v & p" << endl;
      // Velocity, Pressure & Pressure gradient
      xEvalField<xIdentity<xVector> > eval_velo(V);
      xEvalField<xIdentity<double> > eval_press(P);
      if (debug) cout << "v & p" << endl;
      Export(eval_velo, pexport, "VELOCITY", intrule_uu, domain_beg, domain_end);
      if (debug) cout << "v & p" << endl;
      Export(eval_press, pexport, "PRESSURE", intrule_up, domain_beg, domain_end);
      
      if (details > 0)
	{
	  if (debug) cout << "stress export" << endl;
	  SmoothMaterialVariable_c < xTensor2 > sigmatot("stress_total", *variab_manager);
	  // Export(sigmatot, pexport, "Total_Stress", intrule_uu, domain_beg, domain_end );     
	  SmoothMaterialVariable_c < xTensor2 > sigmavis("stress_viscous", *variab_manager);
	  // Export(sigmavis, pexport, "Extra_Stress", intrule_uu, domain_beg, domain_end );
	  
	  
	  if (debug) cout << "components" << endl;
	  // Velocity components & velocity norm
	  xExtractCompVector extr_x(0);  xExtractCompVector extr_y(1); xExtractCompVector extr_z(2);
	  xEvalField<xExtractCompVector > eval_vlx(V, extr_x);
	  xEvalField<xExtractCompVector > eval_vly(V, extr_y);
	  xEvalField<xExtractCompVector > eval_vlz(V, extr_z);
	  Export(eval_vlx, pexport, "V_X", intrule_uu, domain_beg, domain_end);
	  Export(eval_vly, pexport, "V_Y", intrule_uu, domain_beg, domain_end);
	  if (spaceDim == 3) Export(eval_vlz, pexport, "V_Z", intrule_uu, domain_beg, domain_end);
	  xEvalField<xExtractNormVector > eval_vnorm(V);
	  Export(eval_vnorm, pexport, "V_norm", intrule_uu, domain_beg, domain_end);
	  
	  
	  if (debug) cout << "strain rate" << endl;
	  // taux de deformation
	  xExtractCompTensor extr_xx(0,0);  xExtractCompTensor extr_xy(0,1);  xExtractCompTensor extr_yy(1,1);
	  xExtractCompTensor extr_zz(2,2);  xExtractCompTensor extr_xz(0,2);  xExtractCompTensor extr_yz(1,2);
	  
	  xEvalGradField<xExtractCompTensor > eval_Dxx(V, extr_xx);
	  xEvalGradField<xExtractCompTensor > eval_Dyy(V, extr_yy);
	  xEvalGradField<xExtractCompTensor > eval_Dxy(V, extr_xy);
	  xEvalGradField<xExtractCompTensor > eval_Dzz(V, extr_zz);
	  xEvalGradField<xExtractCompTensor > eval_Dxz(V, extr_xz);
	  xEvalGradField<xExtractCompTensor > eval_Dyz(V, extr_yz);
	  
	  Export(eval_Dxx, pexport, "D_XX", intrule_uu, domain_beg, domain_end);
	  Export(eval_Dyy, pexport, "D_YY", intrule_uu, domain_beg, domain_end);
	  if(spaceDim == 3) Export(eval_Dzz, pexport, "D_ZZ", intrule_uu, domain_beg, domain_end);
	  Export(eval_Dxy, pexport, "D_XY", intrule_uu, domain_beg, domain_end);
	  if(spaceDim == 3) 
	    {
	      Export(eval_Dxz, pexport, "D_XZ", intrule_uu, domain_beg, domain_end);
	      Export(eval_Dyz, pexport, "D_YZ", intrule_uu, domain_beg, domain_end);
	    }
	  
	  //SmoothMaterialVariable_c<double> Deq("equ_strain_rate", *variab_manager);
	  //Export(Deq, pexport, "equ_strain_rate", intrule_uu, domain_beg, domain_end );
	  
	  
	  if(details > 1)
	    {
	      
	      if (debug) cout << "grad_p" << endl;
	      xEvalGradField<xIdentity<xVector> > eval_grad_press(P);
	      Export(eval_grad_press, pexport, "Pressure_Grad", intrule_up, domain_beg, domain_end);
	      
	      
	      
	      if (debug) cout << "div v" << endl;
	      // divergence du champ de vitesses
	      xEvalGradField<xTrace > eval_divV(V);
	      Export(eval_divV, pexport, "Div_V", intrule_uu, domain_beg, domain_end);
	      
	      // For some type of constitutive law (newtonian for instance) viscosity doesn t need to
	      // be treated as an internal variable (but a fixed property for instance). Then it can't be exported
	      //SmoothMaterialVariable_c<double> source("viscosity", *variab_manager);
	      //Export(source, pexport, "viscosity", intrule_uu, domain_beg, domain_end );
	    }
	}  
      pexport.closeFile();
    }
  return;
}

// ========================================================================================================
void xpFluidMechanics :: TreatmentOfEssEnv(xEntityFilter  filter, bool firstDeclaration){
  bool debug = false;
  if(!firstDeclaration) resetDirichletDofs();
  for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin(); 
                                   it != thedata->PhysicalEnv->end(); ++it) 
    {
      const xEnv& env = *it;
      string phys = env.Phys; 
      int type = env.Type; 
      int entity = env.Entity;
      double VAL;
      if (phys == "VELOCITY_X" || phys == "VELOCITY_Y" || phys == "VELOCITY_Z" ) {       
	if (type == FIX) 
	  {
	    xClassRegion bc(thedata->mesh, entity, env.getDimension());
	    xFilteredRegion<xClassIter, xEntityFilter> whereToApply(bc.begin(), bc.end(), filter);
	    if (env.hasAnEvolution())
	      {
		cout << "Imposed velocity evolution : ";
		double time = (*Multphys_time_pilot)();
		double VAL = env.getEvolution()(time);
		cout << "at time " << time << ", Vd=" << VAL << endl;
		DirichletBoundaryCondition (Vlin, phys, whereToApply.begin(), whereToApply.end(), VAL);
		DirichletBoundaryCondition (Vbub, phys, whereToApply.begin(), whereToApply.end(), 0.0);  
	      }
	    else
	      {
		VAL = env.getValue();
		DirichletBoundaryCondition (Vlin, phys, whereToApply.begin(), whereToApply.end(), VAL);
		DirichletBoundaryCondition (Vbub, phys, whereToApply.begin(), whereToApply.end(), 0.0); 
	      }
	  }
	else assert(1 == 0);      
      }
    } // End loop over the environment info
  declareUsefulDofs();


  if (debug)      double_manager.PrintForDebug("VafterEssEnv.dbg");
  return;
}

// ========================================================================================================
void xpFluidMechanics :: TreatmentOfNatEnv(xAssembler& assembler)
{
  for (xPhysicalEnv::const_iterator it = thedata->PhysicalEnv->begin();
  		 it != thedata->PhysicalEnv->end(); ++it)
  {
    const xEnv& env = *it;
    if (env.Phys == "TRACTION_X" || env.Phys == "TRACTION_Y"  || env.Phys == "TRACTION_Z" )
    {
      assert(env.Type == FIX);
      xVector val;
      if (env.hasAnEvolution())
	  {
	  	double time = (*Multphys_time_pilot)();
	  	if (env.Phys == "TRACTION_X") val(0) =  env.getEvolution()(time);
	  	if (env.Phys == "TRACTION_Y") val(1) =  env.getEvolution()(time);
	  	if (env.Phys == "TRACTION_Z") val(2) =  env.getEvolution()(time);
	  }
      else
	  {
	  	if (env.Phys == "TRACTION_X") val(0) =  env.getValue();
	  	if (env.Phys == "TRACTION_Y") val(1) =  env.getValue();
	  	if (env.Phys == "TRACTION_Z") val(2) =  env.getValue();
	  }
	  xEvalConstant<xVector>  flux(val);
	  xFormLinearWithLoad<xValOperator<xIdentity<xVector> >, xEvalConstant<xVector> > lin(flux); 
	  xClassRegion bc(thedata->mesh, env.Entity, env.getDimension());
	  Assemble(lin, assembler, intrule_uu, V, bc.begin(), bc.end(), xUpperAdjacency()); 
	}
  } 
  return;
}

// ========================================================================================================

xLevelSet xpFluidMechanics::compute_vnorm( xLevelSet lsmat) {
  // returns a levelset object containing the normal velocity to the isozero of lsmat. Useful for propagating a levelset.
 
  OldAndCurrent_c::current();
  const bool debug = false; 
  if (debug) cout << "chargement de surf dans la levelset lsmat." << endl;
  //xLevelSet lsmat(all);
  //lsmat.load(surf);
  if (debug) cout << "definition de la levelset de vitesses normales" << endl;
  xLevelSet lsvn(all);
  double norm_n, vn;
  xVector  n, Vitesse;
  
  if (debug) cout << "boucle sur l'ensemble des noeuds" << endl;
  for(xIter it = all.begin(0); it != all.end(0); ++it)
    {
      mVertex *e = (mVertex*) *it;
      
      int number_of_adjacent = e->size(spaceDim);
	
      if (number_of_adjacent==0 )
	{
	  cout << "The vertex number " << e->getId() <<  " has no element linked to it" << endl;
	}

      // Get the value Of the velocity for the first adjacent (same for ant adjacent),
      mAdjacencyContainer::iter adjacent_element=e->begin(spaceDim);
      xGeomElem geo_e(*adjacent_element);
      xGeomElem geo_appro(*adjacent_element);
      mPoint pg = e->point(); 
      // on cherche la valeur du champ  au point pg
      geo_appro.setUVWForXYZ(pg);
      Vitesse=V.getVal(&geo_appro,&geo_e, Vitesse);
      
      //ON MOYENNE LES NORMALES DE TOUS LES ELEMENTS ADJACENTS
      //iteration sur l'ensemble des elements adjacents a e
      for (; adjacent_element != e->end(spaceDim); ++adjacent_element )
	{
	  //on somme les normales de tous les adjacents:
	  n +=  lsmat.getGrad(*adjacent_element);//normale à l'iso de levelset
	}
      n /= number_of_adjacent;

      if (debug) cout << "veteur normal: "<<n(0)<<"   "<<n(1)<<endl;
      norm_n = sqrt(n*n);// renormalisation de la normale.
      n = n/norm_n;
      
      // produit scalaire V.n
      vn = Vitesse*n;
      if (debug) cout << "vitesse normale: " << vn << "   Vx=" << Vitesse(0) <<"  Vy=" << Vitesse(1)<<endl;
      //enregistrement de la valeur de lsvn au point e
      lsvn(e)=vn;
    }

  if (debug) 
    {
      xExportGmshAscii  pexport;
      pexport.openFile("NORMAL_VELOC");
      Export(lsvn, pexport, "vitesses_normales");
      pexport.closeFile();
    }
  return lsvn;
}


//==================================================

xLevelSet xpFluidMechanics::propagateLevelSet(xLevelSet lsmat, xLevelSet ls_vitesses_normale, double tolr)
{
   // propagate a levelset lsmat using the normalvelocity field ls_vitesses_normale and the current time step contained in the Basic_time_step Multphys_time_pilot.

  bool debug = false;
  
  //manage different entities containing the time step. xPilotOneStep useful for levelsetOperator. BasicStep used in the xmMultPhys Appli
  xPilotOneStep my_time_step(Multphys_time_pilot->getTimeStep());
  // if (debug)  cout << Multphys_time_pilot->getTimeStep() << endl;

  // All those can be found in xLevelSetOperators.h
  xPropagationOperator my_ls_propagator(ls_vitesses_normale);
  xTimeIntegration my_time_integrator(my_ls_propagator, my_time_step);

  cout << "LEVEL SET PROPAGATION - using Hamilton Jacobi method" << endl << endl;
  lsmat.accept(my_time_integrator);

  
  cout << "   Reinitialisation" << endl;
  int max_it = 10;
  double coeff = 1 ;
  xL2norm l2rei, l2reo;
  double dtvirt = xCFL::getDt(all);
  
  xPilotError statio(l2rei, tolr, max_it, dtvirt*coeff);
  xReInitOperator    reinit;
  xEvolveToStationary rei(reinit, statio);
  lsmat.accept(rei ,all);
  
    
  if (debug) 
    {
      xExportGmshAscii  pexport;
      pexport.openFile("LS" + Multphys_time_pilot->step2string());
      cout << "LEVEL SET PROPAGATED - Exporting the new one" << endl;
      Export(lsmat,pexport ,"LS_" + Multphys_time_pilot->step2string() );
      pexport.closeFile();
    }
  return lsmat;
}
  

xLevelSet xpFluidMechanics::propagateLevelSet(xLevelSet lsmat)
{
  // propagate a levelset lsmat using the velocity field and the current time step contained in the Basic_time_step Multphys_time_pilot.

  //defines a levelset that contain the normal velocity to lsmat
  xLevelSet ls_vitesses_normale = compute_vnorm(lsmat);
  // Call the previous function:
  return propagateLevelSet(lsmat, ls_vitesses_normale);
}




/*
void  xpFluidMechanics::updateDomain(xLevelSet lsmat)
{
  xpPhysicalFormulation::updateDomain(lsmat);
  xExportGmshAscii  pexport;
  pexport.openFile("LS" + Multphys_time_pilot->step2string());
  Export(lsmat,pexport ,"LS_" + Multphys_time_pilot->step2string() );
  pexport.closeFile();
}
*/
/*
void  xpFluidMechanics::updateDomain(xLevelSet lsmat, xLevelSet lstable)
{
  xpPhysicalFormulation::updateDomain(lsmat, lstable);
  xExportGmshAscii  pexport;
  pexport.openFile("LS" + Multphys_time_pilot->step2string());
  Export(lsmat,pexport ,"LS_" + Multphys_time_pilot->step2string() );
  pexport.closeFile();
}
*/
