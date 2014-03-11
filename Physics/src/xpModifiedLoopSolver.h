/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPMODCLOOPSOLVER_H
#define _XPSMODLOOPSOLVER_H

#include "xpLoopSolver.h"

using namespace xfem;

/// this class allows to solve successively THREE non linear problems (xpMultPhys entities) untill every residue are simultaneously smaller than the NLsolver tolerance. This is an iterative solver for multiphyscal problems that will be efficient only if different problems are weakly coupled.It is Modified compared to its parent because it loops on the two first physics util both are converged then solves the third one and starts again. Usefull if the last physic does not influence the two first ones.

template < typename SOLVER_TYPE >
class xpModifiedLoopSolver : public xpLoopSolver<SOLVER_TYPE>
{
  
public:
 typedef typename SOLVER_TYPE::vector_type  vector_type;
  typedef vector<xpPhysicalFormulation*>::iterator phys_iterator;
  
  xpModifiedLoopSolver(nonLinearSystemSolverNR<SOLVER_TYPE >& s, int maxit):
    xpLoopSolver<SOLVER_TYPE>(s, maxit)
    { ;} 

  bool solve(vector<xpPhysicalFormulation*> pbToLoopSolve);


};



////////// THE SOLVE METHOD ///////////////////////
template < class SOLVER_TYPE >
bool xpModifiedLoopSolver<SOLVER_TYPE>::solve(vector<xpPhysicalFormulation*> pbToSolve)
{
  bool two_first_conv = false;
  bool last_conv = false;
  bool convOK = false;
  
  for (int i=0 ; i <= xpLoopSolver<SOLVER_TYPE>::NRmaxit ; i++)
  	{	
  		cout << "Multiphysical loop with 3 physics. Loop number : " << i<< endl;
  		//iterative solving of the two first problems
  		two_first_conv = false;
  		vector<xpPhysicalFormulation*>  deux_premiers_pb;
		deux_premiers_pb.push_back(pbToSolve[0]);
		deux_premiers_pb.push_back(pbToSolve[1]);  
		two_first_conv = xpLoopSolver<SOLVER_TYPE>::solve(deux_premiers_pb);
		
		if (!two_first_conv) break;
		//solve the last problem
		pbToSolve[2]->updateInternalVariables();
		xpLoopSolver<SOLVER_TYPE>::NLsolver->solve(pbToSolve[2]);
	
	
	     //Let's say the convergence is Ok and check if every problem are solved.
	     convOK = true;
	     //To proceed, one calculates the residues of all the  problems:
	     for (phys_iterator iter=pbToSolve.begin(); iter!=pbToSolve.end(); ++iter)
	       {
		 vector_type current_residu((*iter)->get_ndofs());
		 //in the optimized framework, the internal variables of the Modified formulation needs to be updated before.
		 (*iter)->updateInternalVariables();
		 //then the residue can be computed.
		 (*iter)->setFunction(current_residu);
		 //if the current residue is to big,  problem is not solved.
		 if (current_residu.two_norm() / (*iter)->getOrderOfMagnitude()   >    xpLoopSolver<SOLVER_TYPE>::ResolutionSchemeTolerance)   convOK = false;  
	       }
	   
	 // if the conv is still Ok, it means that all the problems are solved: go away from the Loop
	 if (convOK)
	   {
	     break;
	   }
		
		
  	}


 //output of the conclusion
  if (convOK)
    {
      cout << "The iterative solver allowed to solve simultaneously the problems : " ; 
      for (phys_iterator iter=pbToSolve.begin(); iter!=pbToSolve.end(); ++iter)
	  {	  cout   << (*iter)->get_name() << " , " ;	}

    }
  
  else
    {
      cout << "The iterative solver didn't reach convergence, the coupling between problems ";
      for (phys_iterator iter=pbToSolve.begin(); iter!=pbToSolve.end(); ++iter)

	  {	  cout   <<  (*iter)->get_name() << " , " ;	}
      cout    << " may be too strong."<< endl;
      
    }
  return (convOK);
}
#endif

