/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPLOOPSOLVER_H
#define _XPLOOPSOLVER_H

#include <string>
#include "xpPhysicalFormulation.h"
#include "xpResolutionScheme.h"

#include "nonLinearSystemSolverNR.h"


using namespace xfem;

/// this class allows to solve successively different non linear physics (xpPhysicalFormulation entities) untill every residuals are simultaneously smaller than the NLsolver tolerance. This is an iterative solver for multiphyscal problems that will be efficient only if different problems are weakly coupled.

template < typename SOLVER_TYPE, typename M, typename V >
class xpLoopSolver : public xpResolutionScheme
{
  
public:
  typedef V  vector_type;
  typedef vector<xpPhysicalFormulation*>::iterator phys_iterator;
  
  xpLoopSolver(nonLinearSystemSolverNR<SOLVER_TYPE,M,V >& s, int maxit):
    NRmaxit(maxit),
    NLsolver(&s),
    xpResolutionScheme(s.get_tol())
  { ;} 

  ~xpLoopSolver(){;}

  bool solve(vector<xpPhysicalFormulation*> pbToLoopSolve);

public:
  bool convOK;

protected:
///maximum number of loops
  int NRmaxit;
  nonLinearSystemSolverNR<SOLVER_TYPE,M,V> *NLsolver;

};



////////// THE SOLVE METHOD ///////////////////////
template < typename SOLVER_TYPE, typename M, typename V >
bool xpLoopSolver<SOLVER_TYPE,M,V>::solve(vector<xpPhysicalFormulation*> pbToSolve)
{
  
  int NbPb = pbToSolve.size();
  int nbiterfinal = 0;
  int ndofs[NbPb];
  vector_type toto;
  vector<vector_type> residu;

  
  for (phys_iterator iter=pbToSolve.begin(); iter!=pbToSolve.end(); ++iter)
    {
      vector_type monresidu( (*iter)->get_ndofs() );
      residu.push_back(monresidu);  
    }
  bool convOK=false;
  
  for (int i=0 ; i <= NRmaxit ; i++)
    {
      for (phys_iterator PbIterator=pbToSolve.begin(); PbIterator != pbToSolve.end(); ++PbIterator)
	{
	  cout << "LOOP SOLVER - Solving the problem   :  " << (*PbIterator)->get_name();

	 printf("    iter: %3d  \n",i);
	 (*PbIterator)->updateInternalVariables();
	 NLsolver->solve(*PbIterator);

	 //All problems need to be solve at least once.
	 if (PbIterator == pbToSolve.end() || i > 0)
	   {	 
	     //Then let's say the convergence is Ok and check if the other problems are solved.
	     convOK = true;

	     //To proceed, one calculates the residues of all the other problems:
	     for (phys_iterator iter=pbToSolve.begin(); iter!=pbToSolve.end(); ++iter)
	       {
		 vector_type current_residu((*iter)->get_ndofs());
		 //in the optimized framework, the internal variables of the specific formulation needs to be updated before.
		 (*iter)->updateInternalVariables();
		 //then the residue can be computed.
		 (*iter)->setFunction(current_residu);
		 //if the current residue is to big,  problem is not solved.

		 if (current_residu.two_norm() / (*iter)->getOrderOfMagnitude()   >    ResolutionSchemeTolerance)
		   //note that the residue normalization is performed by dividing  by the order of magnitude of the problem.
		   {
		     convOK = false;
		   }
	       }
	   }
	 // if the conv is still Ok, it means that all the problems are solved: go away from the pbIterator Loop
	 if (convOK)
	   {
	     break;
	   }
	 // else solve the next problem in the same Loop i:
       }
     //  if conv is Ok go away from the i Loop as well:
     if (convOK)
       {
	 nbiterfinal = i;
	 break;
       }
     // Else start a new i iteration of all problems:
   }


 //output of the conclusion
  if (convOK)
    {
      cout << "The iterative solver allowed to solve simultaneously the problems : " ; 
      for (phys_iterator iter=pbToSolve.begin(); iter!=pbToSolve.end(); ++iter)
	{	  cout   << (*iter)->get_name() << " , " ;	}
      cout << " after " << nbiterfinal << " iterations" << endl << endl;

    }
  
  else
    {
      cout << "The iterative solver didn't reach convergence, the coupling between problems ";
      for (phys_iterator iter=pbToSolve.begin(); iter!=pbToSolve.end(); ++iter)

	{	  cout   <<  (*iter)->get_name() << " , " ;	}
      cout    << " may be too strong."<< endl;
    }
  return convOK;
}
	       


 
#endif

