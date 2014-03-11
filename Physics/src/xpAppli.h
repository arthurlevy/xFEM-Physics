/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPAPPLI_H
#define _XPAPPLI_H

#include "xpPhysicalFormulation.h"
#include "xpResolutionScheme.h"
#include "xpTimePilot.h"
#include "xLevelSet.h"

using namespace xfem;

/// This class defines a multiphysc application. Mainly, it contains a vector of Physical problem (xpPhysicalFormulation) and Resolution Scheme xpResolutionScheme.
class xpAppli {
public :
  xpAppli (const string aname,  xpTimePilot t_p, xpResolutionScheme &solver  ): name(aname) , Timepilot(t_p), Applisolver(&solver) {;}
  /// if no Time pilot is given in the constructor (SteadyState Problem, for instance) it will initialize to 1,1,1.
  xpAppli (const string aname,xpResolutionScheme & solver  ): name(aname), Timepilot(1, 1, 1),Applisolver(&solver)  {;}
  ~xpAppli() {;}
  
  void info();
  /// this method calls each subproblem exportFields method. details is the level of details that should be exported (details = 0 will export temperature only for instance)
  void exportFields(int details, const std::string& extension, bool binary=false, bool sorted = false);
  void add_time_pilot(xpTimePilot t_p)
  {Timepilot = t_p ; return;}

  void add_model(xpPhysicalFormulation &amodel);

  /// this declares each models fields, approximation and internal variables and transfers time pilot and variable manager to each subproblem. It also initialize the time_history.txt file with then date and data files used in the application.
  void init_models(xEntityFilter filter=xAcceptAll());
/// this declares  internal variables and transfers time pilot and variable manager to each subproblem. It also initialize the time_history.txt file with then date and data files used in the application.It is useful if you need to impose an initial condition then you should : 1-declareApprox 2-InitializeFields 3-use this method.
  void init_models_with_approx_already_declared(xEntityFilter filter=xAcceptAll());
  /// Declares approximations of each sub problem do not use with init_models but with init_models_with_approx_already_declared only.
  void declareApproximations();

  ///the transfert methods will assure the time pilot and variable manager of the xmMultphys sub-problems to be the same as the one of the appli.
  void transfer_varmanager();
  void transfer_timepilot();
  
  /// this method  is usefull to update the domain of each problems once the geometry (lsmat) has changed.
  void updateDomain(xLevelSet lsmat);

 // same as above with more materials (and more interfaces)
 void updateDomain(xLevelSet lsmat, xLevelSet tool);
 void updateDomain(xLevelSet lsmat, xLevelSet tool, xLevelSet compo);
  
  bool solve();

public:
  vector<xpPhysicalFormulation*> ListOfModels;
  xVariabManager VarMNGR;
  xpTimePilot Timepilot;
  
private:
  string name;
  xpResolutionScheme* Applisolver;

};
 
#endif
