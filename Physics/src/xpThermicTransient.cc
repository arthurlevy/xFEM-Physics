/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include <fstream>
#include "xAlgorithm.h"
//#include "xLinearSystemSolverLU.h"
//#include "xEnv.h"
#include "xEval.h"
#include "xExportGmsh.h"
#include "xForm.h"
#include "xFemMatrix.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xMaterialSensitivity.h"
#include "xOperators.h"
//#include "xPhysSurf.h"
//#include "xSimpleGeometry.h"
//#include "xSpace.h"
//#include "xValue.h"

#include "xCSRMatrix.h"
#include "xCSRVector.h"

//#include "mPoint.h"
#include "MaterialCommand.h"
//#include "UpdateMaterialVariablesVisitorPartial_c.h"
#include "NonUniformMaterialSensitivity.h"

#include "xpEval.h"
#include "xpThermicTransient.h"


using namespace xfem;
using namespace lalg;
using namespace AOMD;

// ===============================================================================

  void xpThermicTransient::ShiftFieldCurrentToOld()
  {
    bool debug =false;
    Visit(OldEqualCurrentVisitor_c() , double_manager.begin("dofs"),double_manager.end("dofs") );
    OldAndCurrent_c::current();
    if (debug) double_manager.PrintForDebug("oldeqcurrT.dbg");
  }
// ===============================================================================
void xpThermicTransient :: updateInternalVariables() {
  const bool debug=false;
  if (debug) cout <<"----- UPDATING INTERNAL VARIABLES\n";
  
  //when declared, Oldandcurrent is switched to old, switch it to current:
  OldAndCurrent_c::current();

  //All the variables in the xmaterial are computed at "time Theta" according to the time integration scheme
  // This means for example that "temperature" variable aquals : (1-Theta)*T(old) + Theta*T(current) 

  // temperature : computed from current and old value of T via evalThetaField 
  xEvalThetaField< std::_Identity < double >   > val_T(T, Theta);
  SetMaterialVariablesVisitor_c < xEvalThetaField < std::_Identity<double>  > > visitor_T("temperature",val_T);  
  VisitMaterialVariablesCommand_c visit_T(visitor_T, *variab_manager);  

  // temparature gradient : computed from current and old value of T via evalGradThetaField 
  xEvalGradThetaField< std::_Identity<xVector> > grad_temp(T,Theta);
  SetMaterialVariablesVisitor_c<xEvalGradThetaField< std::_Identity<xVector>  > > 
                               visitor_gT("temperature_gradient",grad_temp); 
  VisitMaterialVariablesCommand_c visit_gT(visitor_gT, *variab_manager);  

 // time derivative of temperature needs the time step.
  xEvalTimeDerivativeField<std::_Identity<double> > val_dT(T,Multphys_time_pilot->getTimeStep());
  SetMaterialVariablesVisitor_c<xEvalTimeDerivativeField<std::_Identity<double> > >
                              visitor_dT("temperature_timederiv",val_dT); 
  VisitMaterialVariablesCommand_c visit_dT(visitor_dT, *variab_manager);  

  if (debug) double_manager.PrintForDebug("toto.dbg");

  // element length usefull for the supg coeff determination
  // currently code in xpThermicTransient.h
  if (1) // some time we will be able to do that once only if the mesh doesn not change.
    {
      ElementLengthMaterialVariablesVisitor_c elemlength("element_size", all);
      VisitMaterialVariablesCommand_c visit_elemlength(elemlength, *variab_manager);  
      ApplyCommandOnIntegrationRule(visit_elemlength, intrule_TT, all.begin(), all.end());
      first_update=0;
    }
   

  if (debug) variab_manager->PrintForDebug("test2.dbg");

  if (debug) cout << "----- APPLYING UPDATE RULES" << endl;
  ApplyCommandOnIntegrationRule(visit_T, intrule_TT, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(visit_dT, intrule_TT, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(visit_gT, intrule_TT, domain_beg, domain_end);

  // Other material variables:
  UpdateMaterialVariablesVisitor_c update_visitor(TypeOfProblem);
  VisitMaterialVariablesCommand_c command_update_visitor(update_visitor, *variab_manager); 
  ApplyCommandOnIntegrationRule(command_update_visitor, intrule_TT, domain_beg, domain_end);

  if (debug) variab_manager->PrintForDebug("Var_updatees_thermic.dbg" );
  if (debug) cout << "----- Internal variables updated" << endl;
  return;
}

// ===============================================================================
template <class VECTTYPE>
void xpThermicTransient :: setFunction( VECTTYPE& b ) {

  const bool debug = false;

  if (TypeOfProblem == "diffusion")
    /////////////////////////
    // DIFFUSION FRAMEWORK///
    {
      if (debug) cout << "----- ASSEMBLING THE EQUATIONS SYSTEM / DIFFUSIVE THERMAL PROBLEM" << endl;
      xAssemblerBasic<xCSRMatrix,VECTTYPE, double> assembler_F(b);
      double COEFF_ASSMB = -1.0;
      
      OldAndCurrent_c::current();
      if (debug) cout << "\t SYSTEM FUNCTION : reset equation system\n";
      b.ZeroArray();
      
      if (debug) cout << "\t SYSTEM FUNCTION : imposed fluxes\n";
      assembler_F.setCoeff(COEFF_ASSMB);
      TreatmentOfNatEnv(assembler_F);

      
      if (debug) cout << "\t SYSTEM FUNCTION : capacity terms\n";
      //management of the capacity term: int(roC*dT/dt)
      assembler_F.setCoeff(-COEFF_ASSMB);
      GetMaterialVariable_c<double> dTdt("temperature_timederiv", *variab_manager);
      //xUniformMaterialParameter<double> rho_c("THERMIC_CAPACITY");
      GetMaterialVariable_c<double> rho_c("thermic_capacity", * variab_manager);
      xEvalBinary<xMult<double, double, double> > roC_Tpoint(dTdt,rho_c);
      //      GetMaterialVariable_c<double> roC_Tpoint("rho_c_diffT", * variab_manager);
      xFormLinearWithLoad< xValOperator <  xIdentity < double> >,
	xEvalBinary<xMult<double,double,double> > > capacitif_contrib(roC_Tpoint);
      Assemble( capacitif_contrib, assembler_F, intrule_TT, T, domain_beg, domain_end);
      
       if (debug) cout << "\t SYSTEM FUNCTION : diffusive terms\n";
      // management of the first diffusif term which is : int(q.grad(T*))
      //it is the same as in the static case
      assembler_F.setCoeff(COEFF_ASSMB);
      GetMaterialVariable_c<xVector> flux("heat_flux", * variab_manager);
      xFormLinearWithLoad< xGradOperator< xIdentity<xVector> >, GetMaterialVariable_c<xVector> > 
      diffusive_lin(flux);
      Assemble(diffusive_lin, assembler_F, intrule_TT, T, domain_beg, domain_end);            
      
      if (debug) cout << "\t SYSTEM FUNCTION : terme source: \n";
      // management of the source term which is : int(s.T*)
      assembler_F.setCoeff(COEFF_ASSMB);
      GetMaterialVariable_c < double > source("thermal_source", *variab_manager);
      xFormLinearWithLoad< xValOperator <  xIdentity < double> >,
      GetMaterialVariable_c < double > > source_contrib(source);
      Assemble(source_contrib, assembler_F, intrule_TT, T, domain_beg, domain_end);
    }
  else if ( TypeOfProblem == "convection")

    ///////////////////////////
    // CONVECTION FRAMEWORK ///
    {           
      if (debug) cout << "----- ASSEMBLING THE EQUATIONS SYSTEM / CONVECTIVE THERMAL PROBLEM" << endl;
      xAssemblerBasic<xCSRMatrix, VECTTYPE, double> assembler_F(b);
      double COEFF_ASSMB = 1.0;
      assembler_F.setCoeff(COEFF_ASSMB);

      OldAndCurrent_c::current();
      if (debug) cout << "\t SYSTEM FUNCTION : reset equation system\n";
      b.ZeroArray();

      //---//
      if (debug) cout << "\t SYSTEM FUNCTION : capacity terms\n";
      //management of the capacity term: int(dT/dt)
      GetMaterialVariable_c<double> Tpoint("temperature_timederiv", * variab_manager);
      xFormLinearWithLoad< xValOperator <  xIdentity < double> >,
	GetMaterialVariable_c < double > > capacitif_contrib(Tpoint);
      Assemble( capacitif_contrib, assembler_F, intrule_TT, T, all.begin(), all.end());
      
      //---//
      if (debug) cout << "\t SYSTEM FUNCTION : convective terms\n";
      GetMaterialVariable_c<xVector> GradT("temperature_gradient", * variab_manager);
      GetMaterialVariable_c<xVector> velocity("velocity", * variab_manager);
      xEvalBinary<xMult<xVector,xVector,double> > v_gradT(GradT,velocity);
      xFormLinearWithLoad< xValOperator <  xIdentity < double> >,
	  xEvalBinary<xMult<xVector,xVector,double> > >convective_contrib(v_gradT);
      Assemble( convective_contrib, assembler_F, intrule_TT, T, all.begin(), all.end());
     
      if (SUPG)
	{
	  //we need ts*v*v*GradT
	  GetMaterialVariable_c<xVector> a_GradT("artificial_temperature_gradient", * variab_manager);
	  /* xEvalBinary<xMult<xVector, xVector, double> > v_v(velocity,velocity);
      	  GetMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
	  xEvalBinary<xMult<double,double,double> > ts_v_v(v_v,supgcoeff);
	  xEvalBinary<xMult<xVector,double, xVector> > a_gradT(GradT,ts_v_v);*/
	  xFormLinearWithLoad< xGradOperator <  xIdentity < xVector> >,
	    GetMaterialVariable_c < xVector > > convective_contrib_supg(a_GradT);
	  Assemble( convective_contrib_supg, assembler_F, intrule_TT, T, all.begin(), all.end());

	  
	  //int (grad N tau v dT/dt)
      	  GetMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
	  xEvalBinary<xMult<double,double,double> > tau_dT_dt(Tpoint,supgcoeff);
	  xEvalBinary<xMult<xVector,double,xVector> > tau_v_dT_dt(velocity,tau_dT_dt);
	  xFormLinearWithLoad< xGradOperator <  xIdentity < xVector> >,
	     xEvalBinary<xMult<xVector,double,xVector> > > transient_contrib_supg(tau_v_dT_dt);
	  Assemble(transient_contrib_supg, assembler_F, intrule_TT, T, all.begin(), all.end());
	  
	}
    }
  
  else
    {
      cout << "no type of problem set in the transient thermal problem";
      assert (0);
    }
  return;
}

// ===============================================================================
template <class MATTYPE>
void xpThermicTransient :: setJacobian(MATTYPE& M ) {
  
  const bool debug = false;
  if (debug) cout << "----- ASSEMBLING THE JACOBIAN MATRIX " << endl;
  xAssemblerBasic<MATTYPE, xCSRVector, double> assmb_M(M);
  if (debug) cout << "\t SYSTEM MATRIX   : reset jacobian\n";
  M.SoftZeroMatrix ( );
  double COEFF_ASSMB = 1.0;
  assmb_M.setCoeff(COEFF_ASSMB);
  
  if (TypeOfProblem == "diffusion")
    ///////////////////////////
    /// DIFFUSION FRAMEWORK ///
    {

      

      if (debug) cout << "\t SYSTEM MATRIX   : Capacitive part\n";
      //first part is roC/dt
      NonUniformMaterialSensitivity_c<double> capac1("temperature_time_derivate",* variab_manager);
      xFormBilinearWithLaw < xValOperator < xIdentity < double > >, 
	NonUniformMaterialSensitivity_c < double > ,
	xValOperator < xIdentity < double > > >     capacity_contrib1(capac1);
      assmb_M.setCoeff(1/Multphys_time_pilot->getTimeStep());
      Assemble(capacity_contrib1, assmb_M, intrule_TT, T, T, domain_beg, domain_end);      
      assmb_M.setCoeff(1.);

      if (debug) cout << "\t SYSTEM MATRIX   :capacitive non linear induced terms \n";
      //second part is theta*droc/dT*(Tpoint)
      NonUniformMaterialSensitivity_c<double> capac2("drho_c/dT",* variab_manager);
      xFormBilinearWithLaw<xValOperator < xIdentity < double> >, 
	NonUniformMaterialSensitivity_c<double>,
	xValOperator<xIdentity<double> > >     capacity_contrib2(capac2 );
      assmb_M.setCoeff(Theta);
      Assemble(capacity_contrib2, assmb_M, intrule_TT, T, T, domain_beg, domain_end);      
      assmb_M.setCoeff(1.);

      if (debug) cout << "\t SYSTEM MATRIX   : Diffusion part\n";
      //first part is the classical one (k) 
      NonUniformMaterialSensitivity_c<xTensor2> fourier("temperature_gradient",* variab_manager);
      xFormBilinearWithLaw<xGradOperator< xIdentity < xVector>   >, 
	NonUniformMaterialSensitivity_c<xTensor2>,
	xGradOperator< xIdentity <xVector>  > >     bilin_diffusive1(fourier);
      assmb_M.setCoeff(-Theta);
      Assemble(bilin_diffusive1, assmb_M, intrule_TT, T, T, domain_beg, domain_end);      
      assmb_M.setCoeff(1.);
      
      if (debug) cout << "\t SYSTEM MATRIX   : Diffusive non linear induced part\n";
      //second part implies the conduction sensibility: theta*dk/dT*gradT
      // and is NON-Symetric: <gradN>*....*<N>
      NonUniformMaterialSensitivity_c<xVector> conduc2("dk/dT",* variab_manager);
      xFormBilinearWithLaw<xGradOperator<xIdentity<xVector> >, 
	NonUniformMaterialSensitivity_c<xVector>,
	xValOperator<xIdentity<double> > >     bilin_diffusive2(conduc2);
      assmb_M.setCoeff(Theta);
      Assemble(bilin_diffusive2, assmb_M, intrule_TT, T, T, domain_beg, domain_end);      
      assmb_M.setCoeff(1.);
      
      if (debug) cout << "\t SYSTEM MATRIX   :Jacobian set\n";
    }
  
  else if (TypeOfProblem == "convection")
    {
      ///////////////////////////
      /// CONVECTION FRAMEWORK ///
 
      if (debug) cout << "\t SYSTEM MATRIX   : Capacitive part\n";
      xEvalConstant < double >  capacity_eval(1/Multphys_time_pilot->getTimeStep());
      xFormBilinearWithLaw < xValOperator < xIdentity < double > >, 
	xEvalConstant < double > ,
	xValOperator < xIdentity < double > > >     capacity_contrib(capacity_eval);
      Assemble(capacity_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
      
      //---//
      if (debug) cout << "\t SYSTEM MATRIX   : convective part \n";
      // d(v.gradT)/dT
      GetMaterialVariable_c<xVector> veloc("velocity", * variab_manager);
      xFormBilinearWithLaw<xValOperator < xIdentity < double> >, 
	GetMaterialVariable_c<xVector>,
	xGradOperator<xIdentity<xVector> > >     convective_contrib(veloc);
      assmb_M.setCoeff(Theta);
      Assemble(convective_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
      assmb_M.setCoeff(1.);
      
      if (SUPG)
	{
	  NonUniformMaterialSensitivity_c < xTensor2 >  sensit_supg("supg_term",* variab_manager);
	  /*	  GetMaterialVariable_c<xVector> velocity("velocity", *variab_manager);
	  GetMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
	  xEvalBinary < xTensorProduct > veloc_matrix(velocity,velocity);
	  xEvalBinary<xMult<xTensor2, double, xTensor2 > > sensit_supg(veloc_matrix,supgcoeff);*/
	  xFormBilinearWithLaw < xGradOperator < xIdentity < xVector > >, 
	    NonUniformMaterialSensitivity_c < xTensor2 > ,
	    xGradOperator < xIdentity < xVector > > >     artificial_convective_contrib(sensit_supg);
	  assmb_M.setCoeff(Theta);
	  Assemble(artificial_convective_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
	  assmb_M.setCoeff(1.);

	  
	  //int gradN  1/dt v N
	  GetMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
	  xEvalBinary<xMult<xVector,double,xVector> > tau_v(veloc,supgcoeff);
	  xFormBilinearWithLaw<xGradOperator < xIdentity < xVector> >, 
	    xEvalBinary<xMult<xVector,double,xVector> >,
	    xValOperator<xIdentity<double> > >     transient_supg_contrib(tau_v);
	  assmb_M.setCoeff(1/Multphys_time_pilot->getTimeStep());
	  Assemble(transient_supg_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
	  assmb_M.setCoeff(1.);
	  
	  
	}

      if (debug) cout << "\t SYSTEM MATRIX   :Jacobian set\n";
    }
  else
    {
      cout << "no type of problem set in the transient thermal problem";
      assert (0);
    } 
  return;
}


// ========================================================================================================

void xpThermicTransient :: exportFields(int details, const std::string& extension, bool binary, bool sorted) {

  bool debug = false;
  xExportGmshAscii pexportascii;
  xExportGmshAsciiSort pexportasciisort;
  xExportGmshBinary pexportbin;
  xExport *ppexport = &pexportascii;
  if (sorted) ppexport =&pexportasciisort;
  xExport & pexport = *ppexport;
  OldAndCurrent_c::current();

  ////////ASCII EXPORT///////////
  if(!binary)
    {
      pexport.openFile("Thermic_"+extension);
      
      //to get a better precision with high order interpolation
      //    pexport.setNbSplitDefault(interpolation_degree);
      //    pexport.setNbSplit(interpolation_degree+1);
      
      // Temperature
      xEvalField<xIdentity<double> > eval_temp(T);
      Export(eval_temp, pexport, "TEMPERATURE", intrule_TT, domain_beg, domain_end);
      
      if (details > 0)
	{
	  //heat Flux
	  SmoothMaterialVariable_c < xVector > eval_flux("heat_flux", *variab_manager);
	  Export(eval_flux, pexport, "HEAT FLUX", intrule_TT, domain_beg, domain_end);
	  
	  // Temperature time derivative
	  xEvalTimeDerivativeField<std::_Identity<double> > val_dT(T,Multphys_time_pilot->getTimeStep()); 
	  Export(val_dT, pexport, "TEMP_deriv", intrule_TT, domain_beg, domain_end);     
	  
	  
	  if (details > 1)
	    {
	      //  Temperature gradient
	      xEvalGradField<xIdentity<xVector> > eval_gradtemp(T);
	      Export(eval_gradtemp, pexport, "TEMPERATURE_GRADIENT", intrule_TT, domain_beg, domain_end);
	      
	      // Internal source
	      SmoothMaterialVariable_c<double> sousource("thermal_source", *variab_manager);
	      Export(sousource, pexport, "SOURCE", intrule_TT, domain_beg, domain_end);
	      
	      if (details > 2)
		{ 
		  //  Temperature gradient
		  xEvalGradField<xIdentity<xVector> > eval_gradtemp(T);
		  Export(eval_gradtemp, pexport, "TEMPERATURE_GRADIENT", intrule_TT, domain_beg, domain_end);
		  // temp gradient components
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
		  if (SUPG && TypeOfProblem =="convection")
		    {
		      // supg coeff
		      SmoothMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
		      Export(supgcoeff, pexport, "SUPG_TIME", intrule_TT, domain_beg, domain_end);
		      if (details > 3)
			{
			  SmoothMaterialVariable_c<double> el_size("element_size", *variab_manager);
			  Export(el_size, pexport, "ELEM_SIZE", intrule_TT, domain_beg, domain_end);
			}
		    } 
		  /////////
		   SmoothMaterialVariable_c<double> roc("thermic_capacity", *variab_manager);
		   Export(roc, pexport, "RHO_C", intrule_TT, domain_beg, domain_end);
		  
		  /////////
		   if (details > 4)
		     {
		       SmoothMaterialVariable_c<double> tmat("temperature", *variab_manager);
		       Export(tmat, pexport, "TEMP_mat", intrule_TT, domain_beg, domain_end);
		     }
		}
	    }
	}
      pexport.closeFile();
    }


  /////BINARY EXPORT//////
  else
    {
      pexportbin.openFile("Thermic_"+extension);
      
      //to get a better precision with high order interpolation
      pexportbin.setNbSplitDefault(interpolation_degree);
      pexportbin.setNbSplit(interpolation_degree);
      
      // Temperature
      xEvalField<xIdentity<double> > eval_temp(T);
      Export(eval_temp, pexportbin, "TEMPERATURE", intrule_TT, domain_beg, domain_end);
      
      if (details > 0)
	{
	  //heat Flux
	  SmoothMaterialVariable_c < xVector > eval_flux("heat_flux", *variab_manager);
	  //	  Export(eval_flux, pexportbin, "HEAT FLUX", intrule_TT, domain_beg, domain_end);
	  
	  // Temperature time derivative
	  xEvalTimeDerivativeField<std::_Identity<double> > val_dT(T,Multphys_time_pilot->getTimeStep()); 
	  Export(val_dT, pexportbin, "TEMP_deriv", intrule_TT, domain_beg, domain_end);     
	  
	  
	  if (details > 1)
	    {
	      //  Temperature gradient
	      xEvalGradField<xIdentity<xVector> > eval_gradtemp(T);
	      Export(eval_gradtemp, pexportbin, "TEMPERATURE_GRADIENT", intrule_TT, domain_beg, domain_end);
	      
	      // Internal source
	      SmoothMaterialVariable_c<double> sousource("thermal_source", *variab_manager);
	      //	      Export(sousource, pexportbin, "SOURCE", intrule_TT, domain_beg, domain_end);
	      
	      if (details > 2)
		{ // temp gradient components
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
		  if (SUPG && TypeOfProblem =="convection")
		    {
		      // supg coeff
		      SmoothMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
		      //  Export(supgcoeff, pexportbin, "SUPG_TIME", intrule_TT, domain_beg, domain_end);
		    } 
		}
	    }
	}
      pexportbin.closeFile();
    }
  
  if (debug)
    {
      variab_manager->PrintForDebug("Var_Thermic.dbg"+extension );
      double_manager.PrintForDebug("doubleMngatExport.dbg"+extension);
    }

  return;
}

