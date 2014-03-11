/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include <fstream>
#include "xAlgorithm.h"
#include "xEnv.h"
#include "xEval.h"
#include "xExportGmsh.h"
#include "xForm.h"
#include "xFemMatrix.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xMaterialSensitivity.h"
#include "xOperators.h"
//#include "xPhysSurf.h"
#include "xSimpleGeometry.h"
//#include "xSpace.h"
#include "xValue.h"
//#include "xLevelSetOperators.h"


#include "mPoint.h"
#include "MaterialCommand.h"
//#include "UpdateMaterialVariablesVisitorPartial_c.h"
#include "NonUniformMaterialSensitivity.h"

#include "xpEval.h"
#include "xpLevelSetConvector.h"

#include "xCSRVector.h"
#include "xCSRMatrix.h"


using namespace xfem;
using namespace lalg;
using namespace AOMD;


// ===============================================================================

  void xpLevelSetConvector::ShiftFieldCurrentToOld()
  {
    bool debug =false;
    Visit(OldEqualCurrentVisitor_c() , double_manager.begin("dofs"),double_manager.end("dofs") );
    OldAndCurrent_c::current();
    if (debug) double_manager.PrintForDebug("oldeqcurrT.dbg");
  }
// ===============================================================================
void xpLevelSetConvector :: updateInternalVariables() {
  const bool debug=false;
  if (debug) cout <<"----- UPDATING INTERNAL VARIABLES\n";
  
  //when declared, Oldandcurrent is switched to old, switch it to current:
  OldAndCurrent_c::current(); 

  // levelset : computed from current and old value of T via evalThetaField 
  xEvalThetaField< std::_Identity < double >   > val_T(T, Theta);
  SetMaterialVariablesVisitor_c < xEvalThetaField < std::_Identity<double>  > > visitor_T("levelset",val_T);  
  VisitMaterialVariablesCommand_c visit_T(visitor_T, *variab_manager);  

  // temparature gradient : computed from current and old value of T via evalGradThetaField 
  xEvalGradThetaField< std::_Identity<xVector> > grad_temp(T, Theta);
  SetMaterialVariablesVisitor_c<xEvalGradThetaField< std::_Identity<xVector>  > > 
                               visitor_gT("levelset_gradient",grad_temp); 
  VisitMaterialVariablesCommand_c visit_gT(visitor_gT, *variab_manager);  

 // time derivative of temperature needs the time step.
  xEvalTimeDerivativeField<std::_Identity<double> > val_dT(T,Multphys_time_pilot->getTimeStep());
  SetMaterialVariablesVisitor_c<xEvalTimeDerivativeField<std::_Identity<double> > >
                              visitor_dT("levelset_timederiv",val_dT); 
  VisitMaterialVariablesCommand_c visit_dT(visitor_dT, *variab_manager);  

  if (debug) double_manager.PrintForDebug("toto.dbg");

  // element length usefull for the supg coeff determination
  // currently code in xpLevelSetConvector.h
  if (first_update) // some time we will be able to do that once only if the mesh doesn not change.
    {
    	if (debug) cout <<"-----  compute and store size of each element\n";
    	ElementLengthMaterialVariablesVisitor_c elemlength("element_size", all);
      VisitMaterialVariablesCommand_c visit_elemlength(elemlength, *variab_manager);  
      ApplyCommandOnIntegrationRule(visit_elemlength, intrule_TT, all.begin(), all.end());
      if (debug) cout <<"compute associated characteristic SUPG time\n";
      UpdateMaterialVariablesVisitor_c update_visitor("characteristic_time");
	  VisitMaterialVariablesCommand_c command_update_visitor(update_visitor, *variab_manager); 
	  ApplyCommandOnIntegrationRule(command_update_visitor, intrule_TT, domain_beg, domain_end);  
      first_update=0;
    }
   

  //if (debug) variab_manager->PrintForDebug("test2.dbg");

  if (debug) cout << "----- APPLYING UPDATE RULES" << endl;
  ApplyCommandOnIntegrationRule(visit_T, intrule_TT, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(visit_dT, intrule_TT, domain_beg, domain_end);
  ApplyCommandOnIntegrationRule(visit_gT, intrule_TT, domain_beg, domain_end);

  // update material variables
  UpdateMaterialVariablesVisitor_c update_visitor("convector");
  VisitMaterialVariablesCommand_c command_update_visitor(update_visitor, *variab_manager);
  ApplyCommandOnIntegrationRule(command_update_visitor, intrule_TT, domain_beg, domain_end);  

  //if (debug) variab_manager->PrintForDebug("Var_updatees_thermic.dbg" );
  if (debug) cout << "----- Internal variables updated" << endl;
  return;
}


// ===============================================================================
template <class VECTTYPE>
void xpLevelSetConvector :: setFunction( VECTTYPE& b ) {

  const bool debug = false;

           
      if (debug) cout << "----- ASSEMBLING THE EQUATIONS SYSTEM / CONVECTIVE THERMAL PROBLEM" << endl;
      xAssemblerBasic<xCSRMatrix, VECTTYPE, double> assembler_F(b);
      double COEFF_ASSMB = 1.0; //newton Raphson framework
      assembler_F.setCoeff(COEFF_ASSMB);

      OldAndCurrent_c::current();
      if (debug) cout << "\t SYSTEM FUNCTION : reset equation system\n";
      b.ZeroArray();
	
	  //---//
	  if (debug) cout << "\t SYSTEM FUNCTION : convective flux on boundaries\n";
	  assembler_F.setCoeff(COEFF_ASSMB);
      TreatmentOfNatEnv(assembler_F);

      //---//
      if (debug) cout << "\t SYSTEM FUNCTION : capacity terms\n";
      //management of the capacity term: int(dT/dt)
      GetMaterialVariable_c<double> Tpoint("levelset_timederiv", * variab_manager);
      xFormLinearWithLoad< xValOperator <  xIdentity < double> >,
	  GetMaterialVariable_c < double > > capacitif_contrib(Tpoint);
      Assemble( capacitif_contrib, assembler_F, intrule_TT, T, all.begin(), all.end());
      
      //---//
      if (debug) cout << "\t SYSTEM FUNCTION : convective terms\n";
      GetMaterialVariable_c<xVector> GradT("levelset_gradient", * variab_manager);
      GetMaterialVariable_c<xVector> velocity("modified_velocity", * variab_manager);
      xEvalBinary<xMult<xVector,xVector,double> > v_gradT(GradT,velocity);
      xFormLinearWithLoad< xValOperator <  xIdentity < double> >,
	  		xEvalBinary<xMult<xVector,xVector,double> > >  convective_contrib(v_gradT);
      Assemble( convective_contrib, assembler_F, intrule_TT, T, all.begin(), all.end());
     
     //////////SUPG terms/////////////
	
	  //we need ts*v*v*GradT
	  //GetMaterialVariable_c<xVector> a_GradT("artificial_levelset_gradient", * variab_manager);
	  xEvalBinary<xMult<xVector, xVector, double> > v_v(velocity,velocity);
	  GetMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
	  xEvalBinary<xMult<double,double,double> > ts_v_v(v_v,supgcoeff);
	  xEvalBinary<xMult<xVector,double, xVector> > a_GradT(GradT,ts_v_v);
	  
	  xFormLinearWithLoad< xGradOperator <  xIdentity < xVector> >,
	    xEvalBinary<xMult<xVector,double, xVector> > > convective_contrib_supg(a_GradT);
	  Assemble(convective_contrib_supg, assembler_F, intrule_TT, T, all.begin(), all.end());


	  //int (grad N ts v dT/dt)
 //    GetMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
	  xEvalBinary<xMult<double,double,double> > ts_dT_dt(Tpoint,supgcoeff);
	  xEvalBinary<xMult<xVector,double,xVector> > ts_v_dT_dt(velocity,ts_dT_dt);
	  xFormLinearWithLoad< xGradOperator <  xIdentity < xVector> >,
	     xEvalBinary<xMult<xVector,double,xVector> > > transient_contrib_supg(ts_v_dT_dt);
	  Assemble(transient_contrib_supg, assembler_F, intrule_TT, T, all.begin(), all.end());
 

  return;
}

// ===============================================================================
template <class MATTYPE>
void xpLevelSetConvector :: setJacobian(MATTYPE& M ) {
  
  const bool debug = false;
  if (debug) cout << "----- ASSEMBLING THE JACOBIAN MATRIX " << endl;
  xAssemblerBasic<MATTYPE, xCSRVector, double> assmb_M(M);
  if (debug) cout << "\t SYSTEM MATRIX   : reset jacobian\n";
  M.SoftZeroMatrix ( );
  double COEFF_ASSMB = 1.0;
  assmb_M.setCoeff(COEFF_ASSMB);
  
   
      if (debug) cout << "\t SYSTEM MATRIX   : Capacitive part\n";
      xEvalConstant < double >  capacity_eval(1/Multphys_time_pilot->getTimeStep());
      xFormBilinearWithLaw < xValOperator < xIdentity < double > >, 
	xEvalConstant < double > ,
	xValOperator < xIdentity < double > > >     capacity_contrib(capacity_eval);
      Assemble(capacity_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
      
      //---//
      if (debug) cout << "\t SYSTEM MATRIX   : convective part \n";
      // d(v.gradT)/dT
      GetMaterialVariable_c<xVector> veloc("modified_velocity", * variab_manager);
      xFormBilinearWithLaw<xValOperator < xIdentity < double> >, 
	GetMaterialVariable_c<xVector>,
	xGradOperator<xIdentity<xVector> > >     convective_contrib(veloc);
      assmb_M.setCoeff(Theta);
      Assemble(convective_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
      assmb_M.setCoeff(1.);
      
      if (1)//SUPG terms
	{
	// the next line is OK if xThermoTrans material is used ie. if a convective thermal problem exists in the multiphysical problem to solve. Indeed the sensitivity to "supg_term" is coded in xThermoTrans. In other cases, please remove the next line and uncomment the 4 others.
	  //NonUniformMaterialSensitivity_c < xTensor2 >  sensit_supg("supg_term",* variab_manager);
	  GetMaterialVariable_c<xVector> velocity("modified_velocity", *variab_manager);
	  GetMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
	  xEvalBinary < xTensorProduct > veloc_matrix(velocity,velocity);
	  xEvalBinary<xMult<xTensor2, double, xTensor2 > > sensit_supg(veloc_matrix,supgcoeff);
	  xFormBilinearWithLaw < xGradOperator < xIdentity < xVector > >, 
	     xEvalBinary<xMult<xTensor2, double, xTensor2 > > ,
	    xGradOperator < xIdentity < xVector > > >     artificial_convective_contrib(sensit_supg);
	  assmb_M.setCoeff(Theta);
	  Assemble(artificial_convective_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
	  assmb_M.setCoeff(1.);

	  
	  //int gradN  1/dt v N
//  GetMaterialVariable_c<double> supgcoeff("characteristic_time", *variab_manager);
	  xEvalBinary<xMult<xVector,double,xVector> > tau_v(veloc,supgcoeff);
	  xFormBilinearWithLaw<xGradOperator < xIdentity < xVector> >, 
	    xEvalBinary<xMult<xVector,double,xVector> >,
	    xValOperator<xIdentity<double> > >     transient_supg_contrib(tau_v);
	  assmb_M.setCoeff(1/Multphys_time_pilot->getTimeStep());
	  Assemble(transient_supg_contrib, assmb_M, intrule_TT, T, T, all.begin(), all.end());      
	  assmb_M.setCoeff(1.);
	  
	  
	}

      if (debug) cout << "\t SYSTEM MATRIX   :Jacobian set\n";
   

  return;
}


// ========================================================================================================

void xpLevelSetConvector :: exportFields(int details, const std::string& extension, bool binary) {

  bool debug = false;
  xExportGmshAscii pexport;
  xExportGmshBinary pexportbin;
  OldAndCurrent_c::current();

  ////////ASCII EXPORT///////////
  if(!binary)
    {
      pexport.openFile("LS_as_formulation"+extension);
      // Temperature
      xEvalField<xIdentity<double> > eval_temp(T);
      Export(eval_temp, pexport, "LS", intrule_TT, domain_beg, domain_end);
      
      if (details > 0)
		{
	  	//normal
	   xEvalGradField<xIdentity<xVector> > eval_gradtemp(T);
	   Export(eval_gradtemp, pexport, "LS_normal", intrule_TT, domain_beg, domain_end);
	   
	    SmoothMaterialVariable_c < xVector > gr("levelset_gradient", *variab_manager);
      	Export(gr, pexport, "gradphi", intrule_TT, domain_beg, domain_end);

	    SmoothMaterialVariable_c < xVector > m_v("modified_velocity", *variab_manager);
      	Export(m_v, pexport, "modified veloc", intrule_TT, domain_beg, domain_end);
      		   
		}
	   pexport.closeFile();
    }


  /////BINARY EXPORT//////
  else
    {
   	pexportbin.openFile("LS_as_formulation"+extension);
      // Temperature
      xEvalField<xIdentity<double> > eval_temp(T);
      Export(eval_temp, pexportbin, "LS", intrule_TT, domain_beg, domain_end);
      
      if (details > 0)
		{
	  	//normal
	   xEvalGradField<xIdentity<xVector> > eval_gradtemp(T);
	   Export(eval_gradtemp, pexportbin, "LS_normal", intrule_TT, domain_beg, domain_end);
	   
	    SmoothMaterialVariable_c < xVector > m_v("modified_velocity", *variab_manager);
      	Export(m_v, pexportbin, "modified veloc", intrule_TT, domain_beg, domain_end);
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

// ========================================================================================================
void xpLevelSetConvector :: loadLevelSet(xLevelSet external_LS)
{
	bool debug = 1;
	if (debug) cout << "boucle sur l'ensemble des noeuds" << endl;
	for(xIter it = all.begin(0); it != all.end(0); ++it)
	{
	//node:
	mVertex *node = (mVertex*) *it;
	
	//value at the node:
	double node_value;
	node_value = external_LS(node);
	
	//copy it to the field
	T.setVal(node, node_value);
	}

	// usually, we load the LS, then convct it so it has to be loaded as the "old" value. After conection, we ought to export the "current" value.
	ShiftFieldCurrentToOld();
	return;
}

// ========================================================================================================

xLevelSet xpLevelSetConvector :: getAsLevelSet()
{
	xLevelSet external_LS(all);
	for(xIter it = all.begin(0); it != all.end(0); ++it)
	{
		//node:
		mVertex *node = (mVertex*) *it;
	
		//value at the node:
		double node_value;
		
		mEntity *parent_element=node->get(spaceDim,1); // get 1st element of size SpaceDim linked to node
  		xGeomElem geo(parent_element);
  		mPoint pg = node->point(); 
      	// on cherche la valeur du champ  au point pg
     	geo.setUVWForXYZ(pg);
		node_value = T.getVal(&geo,&geo,node_value);
	
		//copy it to the LS
		external_LS(node) = node_value;
	}
	
	  cout << "   Reinitialisation" << endl;
  	int max_it = 100;
  	double coeff = 0.9 ;
  	xL2norm l2rei, l2reo;
  	double dtvirt = xCFL::getDt(all);
  	double tolr =1e-4;
  
  	xPilotError statio(l2rei, tolr, max_it, dtvirt*coeff);
  	xReInitOperator    reinit;
  	xEvolveToStationary rei(reinit, statio);
  	external_LS.accept(rei ,all);
	
	return external_LS;
}

