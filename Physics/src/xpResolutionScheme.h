/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPRESOLUTIONSCHEME_H
#define _XPRESOLUTIONSCHEME_H 

#include <string>
#include "xpPhysicalFormulation.h"


using namespace xfem;

/// this class defines how to solve a multiphysic problem. It's only purpose is to set a general class from which will heritate the different schemes (iterative sequential...)
class xpResolutionScheme
{
 
public:
  xpResolutionScheme(): ResolutionSchemeTolerance(1) {;}

  xpResolutionScheme(double tol): ResolutionSchemeTolerance(tol) {;}

  /// pattern for herited class
  virtual bool solve(vector<xpPhysicalFormulation*> pbToLoopSolve);
  /// this method simply calls the previous one
  virtual bool solve(xpPhysicalFormulation& problem1, xpPhysicalFormulation& problem2);
  /// this method simply calls the previous one
  bool solve(xpPhysicalFormulation& problem1,xpPhysicalFormulation& problem2, xpPhysicalFormulation& problem3);

protected:
  double ResolutionSchemeTolerance;

};

#endif
