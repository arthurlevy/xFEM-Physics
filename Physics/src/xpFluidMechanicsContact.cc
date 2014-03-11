/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xAlgorithm.h"
#include "xpFluidMechanicsContact.h"
#include "xpEval.h"
#include "xOperators.h"
#include "xExportGmsh.h"
#include "xAlgorithm.h"
#include "xCSRVector.h"
#include "xCSRMatrix.h"

using namespace lalg;
using namespace AOMD;

xpFluidMechanicsContact::xpFluidMechanicsContact(xData *data,xLevelSet& lsmat,xLevelSet& lstool, int interpolation_degree):
  xpFluidMechanics(data, lsmat), LSTool(lstool),  coeff_contact(1), order (2)
{
  xPhysSurf* domaintool = new xPhysSurf(LSTool,  xIdentity<mEntity*>(),  xIdentity<mEntity*>());
  xRegion toto(domaintool->getMesh_bnd());
  toolsurface = toto;
  xExportGmshAscii  pexport;
  pexport.openFile("LS_INIT");
  Export(LSTool, pexport, "LS_CONTACT");
  Export(lsmat, pexport, "LS_00");
  pexport.closeFile();
}

void xpFluidMechanicsContact::setFunction(xCSRVector &F)
{
  const bool debug = true;
  xpFluidMechanics::setFunction(F);
 
  double dt = Multphys_time_pilot->getTimeStep();
  xAssemblerBasic<> assembler_F(F);
  OldAndCurrent_c::current();
  
  // evaluates the distance of the displaced point to the tool
  xEvalLevelSetAtDisplacedPoint distance(LSTool, V, dt);
  // keeps only the negative part (if there is a penetration)
  xEvalUnary<xNegativePart<double> > neg_distance(distance);
  // Normal to the tool
  xEvalGradLevelSet<xIdentity<xVector> >  normal(LSTool);

  
   if (order==2)
    {
      // parabolic (smoother at distance=0) force
      xEvalBinary<xMult<double, double, double> > neg_square_distance1(neg_distance, neg_distance);
      xEvalUnary<xChangeSign <double> > neg_square_distance(neg_square_distance1);
      // artificial force applyed (penetration * normal)
      xEvalBinary<xMult<xVector, double, xVector> > distance_normal(normal , neg_square_distance);

      ///EXPORT
      if (debug)
	{
	  xExportGmshAscii  pexport;
	  pexport.openFile("CONTACT_FUNC" +  Multphys_time_pilot->step2string());
	  Export(distance_normal, pexport, "penetration", 
		 intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	  pexport.closeFile();
	}
      
      xFormLinearWithLoad<xValOperator<xIdentity<xVector> >, 
	xEvalBinary<xMult<xVector, double, xVector> > > contact_contrib(distance_normal); 
      //in order to assemble a force (ksi * penetration * normal)
      assembler_F.setCoeff(coeff_contact);
      // assembling is done using the UpperCreator which are the elements containing the levelSet.
      //Indeed iso_zero are entities of dim n-1.
       Assemble(contact_contrib, assembler_F, intrule_uu, V, iso_zero.begin(), iso_zero.end(), xUpperCreator());
    }

   else if (order == 1)
     {
       // artificial force applyed (penetration * normal)
       xEvalBinary<xMult<xVector, double, xVector> > distance_normal(normal , neg_distance);
       
       ///EXPORT
       if (debug)
	 {
	   xExportGmshAscii  pexport;
	   pexport.openFile("CONTACT_FUNC" +  Multphys_time_pilot->step2string());
	   Export(distance_normal, pexport, "penetration", 
		  intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	   pexport.closeFile();
	}
       
      xFormLinearWithLoad<xValOperator<xIdentity<xVector> >, 
	xEvalBinary<xMult<xVector, double, xVector> > > contact_contrib(distance_normal); 
      //in order to assemble a force (ksi * penetration * normal)
      assembler_F.setCoeff(coeff_contact);
      // assembling is done using the UpperCreator which are the elements containing the levelSet.
      //Indeed iso_zero are entities of dim n-1.
      Assemble(contact_contrib, assembler_F, intrule_uu, V, iso_zero.begin(), iso_zero.end(), xUpperCreator());
    }
   
  else 
    {
      cout << "order of penalization must be 1 or 2" << endl;
      assert(0);
    }
  return;
}


void xpFluidMechanicsContact::setJacobian(xCSRMatrix & J)
{ 
  const bool debug = false;
  xpFluidMechanics::setJacobian(J);  
  double dt = Multphys_time_pilot->getTimeStep();
  xAssemblerBasic<> assmb_M(J);
  //normal to the tool
  xEvalGradLevelSet<xIdentity<xVector> >  normal(LSTool);
  //tensor 2 normal x normal
  xEvalBinary<xTensorProduct > normal_x_normal(normal , normal);
  assmb_M.setCoeff(-coeff_contact*dt);

   if (order==1)
     {
       // returns 0 if distance of tool at displaced point is negative, and 1 if positive (if there is penetration).
       xEvalIfPenetrationAtDisplacedPoint penetration(LSTool, V, dt);
       xEvalUnary<xChangeSign <double> > neg_penetration(penetration);
       xEvalBinary<xMult<xTensor2, double, xTensor2> > sensitivity(normal_x_normal, neg_penetration);

       //EXPORT
       if (debug)
	 {
	   xExportGmshAscii  pexport;
	   pexport.openFile("CONTACT_JACOB"+ Multphys_time_pilot->step2string() );
	   Export(sensitivity, pexport, "sensitivity", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	   pexport.closeFile();
	 }
       
       xFormBilinearWithLaw<xValOperator<xIdentity<xVector> >, 
	 xEvalBinary<xMult<xTensor2, double, xTensor2> >,
	 xValOperator<xIdentity<xVector> > >     contact_contrib(sensitivity);
     
       // assembling is done using the UpperCreator which are the elements containing the levelSet.
       //Indeed iso_zero are entities of dim n-1.
       Assemble(contact_contrib, assmb_M, intrule_uu, V, V, iso_zero.begin(), iso_zero.end(), xUpperCreator());  
     }
   

   else if (order == 2)
     {
       // evaluates the distance of the displaced point to the tool
       xEvalLevelSetAtDisplacedPoint distance(LSTool, V, dt);
       // keeps only the negative part (if there is a penetration)
       xEvalUnary<xNegativePart<double> > neg_distance(distance);
       // our sensitivity
       xEvalBinary<xMult<xTensor2, double, xTensor2> > sensitivity(normal_x_normal, neg_distance);

       //EXPORT
       if (debug)
	 {
	   xExportGmshAscii  pexport;
	   pexport.openFile("CONTACT_JACOB" + Multphys_time_pilot->step2string() );
	   Export(sensitivity, pexport, "sensitivity", intrule_uu, iso_zero.begin(), iso_zero.end(), xUpperCreator());
	   pexport.closeFile();
	 }
       
       xFormBilinearWithLaw<xValOperator<xIdentity<xVector> >, 
	 xEvalBinary<xMult<xTensor2, double, xTensor2> >,
	 xValOperator<xIdentity<xVector> > >     contact_contrib(sensitivity);
       
       // assembling is done using the UpperCreator which are the elements containing the levelSet.
       //Indeed iso_zero are entities of dim n-1.
       Assemble(contact_contrib, assmb_M, intrule_uu, V, V, iso_zero.begin(), iso_zero.end(), xUpperCreator());  
     }
   
   else 
     {
      cout << "order of penalization can be 1 or 2" << endl;
      assert(0);
    }
   
  return;
}



void xpFluidMechanicsContact::setCoeffContact(double allowed_penetration)
{
  int ndofs=get_ndofs();
  xCSRMatrix J(ndofs);
  updateInternalVariables();
  xpFluidMechanics::setJacobian(J);
  double norm;


  for (xCSRMatrix::iterator  it = J.begin() ; it!=J.end() ; ++it )
    {
      norm += (*it ) *(  *it );
    }

  norm = sqrt(norm) ;

  if (order==1)   coeff_contact = abs(norm);
  else if (order==2)   coeff_contact =  abs(norm)/(allowed_penetration);
  else 
    {
      cout << "order of penalization must be 1 or 2" << endl;
      assert(0);
    }
      
  cout << "Fluid Mechanics : Coefficient of contact set to : " << coeff_contact << endl;
  return;
}

