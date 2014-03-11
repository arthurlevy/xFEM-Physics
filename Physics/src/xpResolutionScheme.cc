/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xpResolutionScheme.h"

bool xpResolutionScheme::solve(vector<xpPhysicalFormulation*> pbToLoopSolve)
  {
    cout << "No resolution scheme affected" << endl;
    return false;
  }

bool xpResolutionScheme::solve(xpPhysicalFormulation& problem1, xpPhysicalFormulation& problem2)
  {
  	vector<xpPhysicalFormulation*> my_problems;
  	my_problems.push_back(&problem1);
 	my_problems.push_back(&problem2);
    return solve(my_problems);
  }
  
bool xpResolutionScheme::solve(xpPhysicalFormulation& problem1,xpPhysicalFormulation& problem2, xpPhysicalFormulation& problem3)
    {  
      vector<xpPhysicalFormulation*>  maListedePb;
      maListedePb.push_back(&problem1);
      maListedePb.push_back(&problem2);
      maListedePb.push_back(&problem3);
      return solve(maListedePb);
    }
