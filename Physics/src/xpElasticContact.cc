/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xAlgorithm.h"
#include "xpElasticContact.h"
#include "xpEval.h"
#include "xOperators.h"
#include "xExportGmsh.h"
#include "xAlgorithm.h"
#include "xCSRMatrix.h"
#include "xCSRVector.h"


using namespace xfem;
using namespace lalg;
using namespace AOMD;

xpElasticContact::xpElasticContact (xData *data,xLevelSet& lsmat,xLevelSet& lstool, int interpolation_order, int penalization_order):
  xpElastic(data, lsmat, interpolation_order), LSTool(lstool),  coeff_contact(1), order(penalization_order)
{
  xExportGmshAscii  pexport;
  pexport.openFile("LS_INIT");
  Export(LSTool, pexport, "LS_CONTACT");
  Export(lsmat, pexport, "LS_00");
  pexport.closeFile();
}

void xpElasticContact::setFunction(xCSRVector &F)
{
  const bool debug = false;
  xpElastic::setFunction(F);
  
  xAssemblerBasic<> assembler_F(F);
  OldAndCurrent_c::current();
  
  // evaluates the distance of the displaced point to the tool
  xEvalLevelSetAtDisplacedPoint distance(LSTool, U);
  // keeps only the negative part (if there is a penetration)
  xEvalUnary<xNegativePart<double> > neg_distance(distance);
  // normal to the tool
  xEvalGradLevelSet<xIdentity<xVector> >  normal(LSTool);
  assembler_F.setCoeff(coeff_contact);

  if (order==1)
    {
      // artificial force applied (penetration * normal)
      xEvalBinary<xMult<xVector, double, xVector> > distance_normal(normal , neg_distance);

      //EXPORT
      if (debug)
	{
	  xExportGmshAscii  pexport;
	  pexport.openFile("TEST_FUNCT" );
	  Export(distance, pexport, "distance_at_disp_pnt", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	  Export(neg_distance, pexport, "distance_penet", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator()); 
	  Export(distance_normal, pexport, "penetr", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	  pexport.closeFile();
	}

      xFormLinearWithLoad<xValOperator<xIdentity<xVector> >, 
	xEvalBinary<xMult<xVector, double, xVector> > > contact_contrib(distance_normal); 
      // assembling is done using the UpperCreator which are the elements containing the levelSet.
      //Indeed iso_zero are entities of dim n-1.
      Assemble(contact_contrib, assembler_F, intrule_uu, U, iso_zero.begin(), iso_zero.end(), xUpperCreator());
    }

  else if (order==2)
    {
      // parabolic (smoother at distance=0) force
      xEvalBinary<xMult<double, double, double> > neg_square_distance1(neg_distance, neg_distance);
      xEvalUnary<xChangeSign <double> > neg_square_distance(neg_square_distance1);
      // artificial force applyed (penetration * normal)
      xEvalBinary<xMult<xVector, double, xVector> > distance_normal(normal , neg_square_distance);

      //EXPORT
      if (debug)
	{
	  xExportGmshAscii  pexport;
	  pexport.openFile("CONTACT_FUNC_EL");
	  Export(LSTool, pexport, "LStool")	;
	  xEvalField<xIdentity<xVector> > eval_disp(U);
	  Export(eval_disp, pexport, "U", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	  Export(distance, pexport, "distance_at_disp_pnt", intrule_uu,domain_beg, domain_end );
	  Export(distance, pexport, "distance_at_disp_pnt_surf", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	  Export(neg_distance, pexport, "distance_penet", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator()); 	      
	  Export(distance_normal, pexport, "effort_contact", 
		 intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	  pexport.closeFile();
	}
      
      xFormLinearWithLoad<xValOperator<xIdentity<xVector> >, 
	xEvalBinary<xMult<xVector, double, xVector> > > contact_contrib(distance_normal); 
      // assembling is done using the UpperCreator which are the elements containing the levelSet.
      //Indeed iso_zero are entities of dim n-1.
      Assemble(contact_contrib, assembler_F, intrule_uu, U, iso_zero.begin(), iso_zero.end(), xUpperCreator());
    }
  
  else 
    {
      cout << "order of penalization can be 1 or 2" << endl;
      assert(0);
    }
  
  return;
}


void xpElasticContact::setJacobian(xCSRMatrix & J)
{
  const bool debug = false;
   xpElastic::setJacobian(J);

   xAssemblerBasic<> assmb_M(J);
   //normal to the tool
   xEvalGradLevelSet<xIdentity<xVector> >  normal(LSTool);
   //tensor 2 normal x normal
   xEvalBinary<xTensorProduct > normal_x_normal(normal , normal);


   if (order==1)
     {
       assmb_M.setCoeff(coeff_contact);
       // returns 0 if distance of tool at displaced point is negative, and 1 if positive (if there is penetration.
       xEvalIfPenetrationAtDisplacedPoint penetration(LSTool, U);
       xEvalUnary<xChangeSign <double> > neg_penetration(penetration);
       // our sensitivity
       xEvalBinary<xMult<xTensor2, double, xTensor2> > sensitivity(normal_x_normal, penetration);

       //EXPORT
       if (debug)
	 {
	   xExportGmshAscii  pexport;
	   pexport.openFile("TEST_JACOB" );
	   Export(penetration, pexport, "penetration_bool", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	   Export(normal, pexport, "normal", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());   
	   Export(sensitivity, pexport, "sensitivity", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	   pexport.closeFile();
	 }
						 

       xFormBilinearWithLaw<xValOperator<xIdentity<xVector> >, 
	 xEvalBinary<xMult<xTensor2, double, xTensor2> >,
	 xValOperator<xIdentity<xVector> > >     contact_contrib( sensitivity);
       // assembling is done using the UpperCreator which are the elements containing the levelSet.
       //Indeed iso_zero are entities of dim n-1.
       Assemble(contact_contrib, assmb_M, intrule_uu, U, U, iso_zero.begin(), iso_zero.end(), xUpperCreator());  
     }
   
  else if (order == 2)
     {
       assmb_M.setCoeff(-coeff_contact);
       // evaluates the distance of the displaced point to the tool
       xEvalLevelSetAtDisplacedPoint distance(LSTool, U);
       // keeps only the negative part (if there is a penetration)
       xEvalUnary<xNegativePart<double> > neg_distance(distance);
       // our sensitivity
       xEvalBinary<xMult<xTensor2, double, xTensor2> > sensitivity(normal_x_normal, neg_distance);
       
       //EXPORT
       if (debug)
	 {
	   xExportGmshAscii  pexport;
	   pexport.openFile("CONTACT_JACOB"  );
	   Export(sensitivity, pexport, "sensitivity", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	   pexport.closeFile();
	 }
       
       xFormBilinearWithLaw<xValOperator<xIdentity<xVector> >, 
	 xEvalBinary<xMult<xTensor2, double, xTensor2> >,
	 xValOperator<xIdentity<xVector> > >     contact_contrib(sensitivity);
       // assembling is done using the UpperCreator which are the elements containing the levelSet.
       //Indeed iso_zero are entities of dim n-1.
       Assemble(contact_contrib, assmb_M, intrule_uu, U, U, iso_zero.begin(), iso_zero.end(), xUpperCreator());  
     }
   
   else 
     {
      cout << "order of penalization can be 1 or 2" << endl;
      assert(0);
    }
  return;
}



void xpElasticContact::setCoeffContact(double allowed_penetration)
{
  int ndofs=get_ndofs();
  xCSRMatrix J(ndofs);
  updateInternalVariables();
  xpElastic::setJacobian(J);
  double norm;
  
  for (xCSRMatrix::iterator  it = J.begin() ; it!=J.end() ; ++it )
    {
      norm += (*it ) *( * it );
    }

  norm = sqrt(norm);
  if (order==1)   coeff_contact = abs(norm)/allowed_penetration;
  else if (order==2)   coeff_contact = abs(norm)/(allowed_penetration);
  else 
    {
      cout << "order of penalization must be 1 or 2" << endl;
      assert(0);
    }
      
  cout << "Elasticity : Coefficient of contact set to : " << coeff_contact << endl;
  return;
}
