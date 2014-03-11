/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xpConvection.h"
#include "xpEval.h"
#include "xExportGmsh.h"
#include "xMaterialSensitivity.h"
#include "MaterialCommand.h"
#include "xCSRMatrix.h"
#include "xCSRVector.h"

using namespace lalg;


void xpConvection::initialize()
{
	declareApproximation();
    TreatmentOfEssEnv(xAcceptAll(), true);
    declareInternalVariables();
    if(renumbered) pair<int, int> res_renum = renumber();
    computeElementLength();
    updateInternalVariables();
    return;
}

 void xpConvection::setFunction(xCSRVector& F)
  {
    //parent calculation
    xpThermic::setFunction(F);
    
    bool debug =false;
    xAssemblerBasic<> assembler_F(F);
    double COEFF_ASSMB = 1.0; 
    assembler_F.setCoeff(-COEFF_ASSMB);
    
    if (debug) cout << "\t SYSTEM FUNCTION : convective terms\n";
    GetMaterialVariable_c<xVector> GradT("temperature_gradient", * variab_manager);
    GetMaterialVariable_c<xVector> veloc("velocity", * variab_manager);
    xUniformMaterialParameter<double> rho_c("THERMIC_CAPACITY");
    xEvalBinary<xMult<xVector,xVector,double> > v_gradT(GradT, veloc);
    xEvalBinary<xMult<double,double,double> > rho_c_v_gradT(rho_c, v_gradT);
    xFormLinearWithLoad< xValOperator <  xIdentity < double> >,
      xEvalBinary<xMult<double,double,double> > >convective_contrib(rho_c_v_gradT);
    Assemble( convective_contrib, assembler_F, intrule_TT, T, all.begin(), all.end());
    
    if (SUPG)
      {
	if (debug) cout << "\t SYSTEM FUNCTION : SUPG terms\n";
	assembler_F.setCoeff(-COEFF_ASSMB);
	xEvalBinary < xTensorProduct > veloc_matrix(veloc,veloc);
	GetMaterialVariable_c<double> ts("characteristic_time", *variab_manager);
	xEvalBinary<xMult<double,double,double> > rho_c_ts(rho_c, ts);
	xEvalBinary<xMult<xVector,double,xVector> > rho_c_ts_GradT(GradT,rho_c_ts);
	xEvalBinary<xMult<xTensor2,xVector,xVector> > supgterm(veloc_matrix,rho_c_ts_GradT);
	xFormLinearWithLoad< xGradOperator <  xIdentity < xVector> >,
	  xEvalBinary<xMult<xTensor2,xVector,xVector> > > convective_contrib_supg(supgterm);
	Assemble( convective_contrib_supg, assembler_F, intrule_TT, T, all.begin(), all.end());

	if (debug) cout << "\t SYSTEM FUNCTION : SUPG induced source terms\n";
	assembler_F.setCoeff(COEFF_ASSMB);
	GetMaterialVariable_c<double> source("thermal_source", * variab_manager);
	xEvalBinary<xMult<double,double,double> > ts_source(source,ts);
	xEvalBinary<xMult<xVector,double,xVector> > supgsourceterm(veloc,ts_source);
	xFormLinearWithLoad< xGradOperator <  xIdentity < xVector> >,
	  xEvalBinary<xMult<xVector,double, xVector> > > source_supg(supgsourceterm);
	Assemble( source_supg, assembler_F, intrule_TT, T, all.begin(), all.end());
      }
    return; 
  }

 void xpConvection::setJacobian(xCSRMatrix&  J)
  {
    bool debug =false;
    // parent calculation
    xpThermic::setJacobian(J);

    xAssemblerBasic<> assmb_M(J);
    double COEFF_ASSMB = 1.0;
    assmb_M.setCoeff(COEFF_ASSMB);

    if (debug) cout << "\t SYSTEM MATRIX   : convective part \n";
    // rho_c * d(v.gradT)/dT
    GetMaterialVariable_c<xVector> veloc("velocity", * variab_manager);
    xUniformMaterialParameter<double> rho_c("THERMIC_CAPACITY");
    xEvalBinary<xMult<xVector,double,xVector> > rho_c_v(veloc,rho_c);

    xFormBilinearWithLaw<xValOperator < xIdentity < double> > , 
      xEvalBinary<xMult<xVector,double,xVector> >,
      xGradOperator<xIdentity<xVector> > >     convective_jac_contrib(rho_c_v);
    assmb_M.setCoeff(-COEFF_ASSMB);
    Assemble(convective_jac_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());       
    
    if(SUPG)
      {
	if (debug) cout << "\t SYSTEM MATRIX   : SUPG part \n";
	// rho_c_ts_V_x_V
	xEvalBinary < xTensorProduct > veloc_matrix(veloc,veloc);
       	GetMaterialVariable_c<double> ts("characteristic_time", *variab_manager);
	xEvalBinary<xMult<double,double,double> > rho_c_ts(rho_c, ts);
	xEvalBinary<xMult<xTensor2,double,xTensor2> > supgterm(veloc_matrix,rho_c_ts);
	xFormBilinearWithLaw < xGradOperator < xIdentity < xVector > >, 
	  xEvalBinary<xMult<xTensor2,double,xTensor2> >,
	  xGradOperator < xIdentity < xVector > > >     artificial_convective_contrib(supgterm);
	Assemble(artificial_convective_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
      }
    return;
  }

void  xpConvection::computeElementLength()
{
  const bool debug=false;
  if (debug) cout <<"----- COMPUTING ELEMENTS LENGTH"<<endl;
  OldAndCurrent_c::current();
  ElementLengthMaterialVariablesVisitor_c elemlength("element_size", all);
  VisitMaterialVariablesCommand_c visit_elemlength(elemlength, *variab_manager);  
  ApplyCommandOnIntegrationRule(visit_elemlength, intrule_TT, all.begin(), all.end()); 
  if (debug) variab_manager->PrintForDebug("elements_length.dbg");
  return;
}

void  xpConvection::setVelocity(xEvalConstant<xVector> my_eval)
{
  SetMaterialVariablesVisitor_c< xEvalConstant<xVector > > velocvisitor("velocity",my_eval);
  VisitMaterialVariablesCommand_c  visit_velo(velocvisitor, *variab_manager);
  ApplyCommandOnIntegrationRule(visit_velo, intrule_TT, all.begin(), all.end());
  updateInternalVariables();
}


void xpConvection :: exportFields(int details, const std::string& extension, bool binary,bool sorted) {
  const bool debug = false;
  xExportGmshAscii  pexportascii;  
  xExportGmshAsciiSort  pexportasciisorted;  
  xExport *ppexport = &pexportascii;
  if (sorted) ppexport = &pexportasciisorted;
  xExport &pexport = *ppexport;
  OldAndCurrent_c::current();
  pexport.openFile("Thermic_"+ extension);

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
	  
	  // Internal source
	  SmoothMaterialVariable_c<double> sousource("thermal_source", *variab_manager);
	  Export(sousource, pexport, "SOURCE", intrule_TT, domain_beg, domain_end);

	  SmoothMaterialVariable_c<xVector> velocity("velocity", *variab_manager);
	  Export(velocity, pexport, "VELO", intrule_TT, domain_beg, domain_end);

	  SmoothMaterialVariable_c<double> TS("characteristic_time", *variab_manager);
	  Export(TS, pexport, "Tsupg", intrule_TT, domain_beg, domain_end);
	  
	}
    }
  pexport.closeFile();

  return;
}
