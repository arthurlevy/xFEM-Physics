/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _steps__pilot__h
#define _steps__pilot__h

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <string>

#include "xLevelSet.h"
#include "xLevelSetOperators.h"

using namespace std;

using namespace xfem;

/// This class defines a time step iterator usefull  for transcient problems
class xpTimePilot 
{
public:
  /// results will be exported on every k iterations
  xpTimePilot(double time, int iter_m, int k=1)
    : final_time(time),current_time(0.), iter_max(iter_m), iter_current(0), iter_save(1), i_save(k), current_convergence(false){}

  
  bool finished();
  ///increment time with time step dt
  void increment_time(const double dt);
  ///increment time with time step dt
  void operator+=(const double dt) {increment_time(dt);}
  ///increment time with the value final_time/iter_max
  void operator++();
  /// returns true if the current time step is supposed to be exported (as defined with i_save)
  bool export_step();
  void DisplayStepInfo();
  /// in the time_history.txt file
  void SaveStepInfo();
  void DisplayGeneralInfo();
  
  ///increment the time_pilot using the CFL condition applied on ls_vitesses_normale, a leveset containing normal velocity.
  /*! We can also give a boundary region to restrict the calculation of the CFL condition on the iso_zero of the levelset.*/
  void increment_with_CFL_Condition(xLevelSet ls_vitesses_normale, xRegion bnd = xRegion());

  /*! This method increments the time_pilot with the CFL condition on the levelSet ls_vitesses_normale containing a normal velocity 
     OR 
  the specified_dt time step, if it is smaller.
  It returns true if CFL condition imposes the time step.
  */
  bool increment_with_CFL_Condition(xLevelSet ls_vitesses_normale,  double specified_dt, xRegion bnd = xRegion());

  /// returns the time step as a 2 digit string
  string step2string() const ; 
  
  /// divide the time step by a factor fac without chaging anything else. Useful in a case of non-convergence for startegy change.
  void reduce_time_step(double fac =0.5)
  {
    current_time += current_time_step*(fac-1);
    current_time_step *= fac;
    current_type_of_increment += " - reduced";
  }

  

  /// -------returns current time-------
  double operator()(void) {return current_time;}
  double getTimeStep() {  return current_time_step; }
  double getFinalTime() {return final_time;}
  bool getCurrentConvergence() { return current_convergence;}
  void setCurrentConvergence(bool conv) { current_convergence = conv;}
  string getType_of_increment() {return current_type_of_increment;}
  void setType_of_increment(string typ) {current_type_of_increment = typ;}
  int getIter() {return iter_current;}
  int getIterMax() {return iter_max;}

 private:
  ///
  int iter_max;
  /// maximum final time
  double final_time;
  bool current_convergence;
  string current_type_of_increment;
  double current_time;
  ///dt
  double current_time_step;
  /// time step number
  int iter_current;
  
  /// export solution every i_save step
  int i_save;
  /// how many steps have been performed since last saving step
  int iter_save; 
};


#endif

