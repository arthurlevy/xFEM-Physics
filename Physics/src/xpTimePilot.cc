/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xpTimePilot.h"

//using std::vector;
//using namespace std;
 

bool xpTimePilot::finished()
{
  if ((current_time >= final_time)||(iter_current > iter_max )) return true;
  else return false;
}  

void xpTimePilot::increment_time(const double dt)
{
  if (iter_current > 0) SaveStepInfo();
  current_time += dt;
  iter_current += 1;
  iter_save += 1;
  current_time_step = dt;
  current_type_of_increment = "user defined";
  current_convergence=false;
  if (iter_save > i_save)   iter_save = 1; 
}

bool xpTimePilot::export_step() 
{
  return (iter_save == 1||iter_current < 5);
}

void xpTimePilot::operator++()
{
  increment_time ( final_time/iter_max ) ; 
  current_type_of_increment = "constructor defined" ;  
}

void  xpTimePilot::DisplayStepInfo()
{
  cout << "================================= Current Time Step ================================================\n";
  cout << " Increment : " << iter_current <<"  time : " << current_time << " time step : " << current_time_step << " type of increment : " << current_type_of_increment <<endl  ; 
  cout << "====================================================================================================\n";
}

void  xpTimePilot::SaveStepInfo()
{
  fstream filestr;
  ofstream toto("time_history.txt", ios::app);
  toto.width(3); toto << right  << iter_current ;
  toto <<" :  time : ";
  toto.width(8);  toto.precision(4);  toto << left << fixed <<  current_time;
  toto << " | time step : ";
  toto.width(8); toto.precision(4); toto << left << fixed << current_time_step;
  toto << " | "<< current_type_of_increment << " increment " << " |  Previous step " ;
  if (current_convergence) {toto << " did " ;}
  else {toto << " -> DID NOT <- ";}
  toto << "converge." << endl;
}


void  xpTimePilot::DisplayGeneralInfo()
{
  cout << "======================== Time Pilot General Settings ==========================\n";
  cout << "\t\t Total time                 : "<<final_time<<endl;
  cout << "\t\t Total number of increments : "<<iter_max+1<<endl;
  if (iter_current==0)
    {
      current_time_step = final_time/iter_max;
    }
  cout << "\t\t Estimated time step        : "<<current_time_step<<endl;
  cout << "============================= Time Step History =================================\n";
 
}


string xpTimePilot::step2string() const 
    {  
      char str[30];
      // we would prefer iter 2 written as 02, so that alpabetically, 09 is before 10
      std::sprintf(str, "%02d",iter_current );
      //the 0 above is the filling digit and the 2 is the imposed number of digit
      return std::string(str);
    }


//=================================================

void  xpTimePilot::increment_with_CFL_Condition( xLevelSet ls_vitesses_normale, xRegion bnd)
{
  // if no bnd was given as argument then bnd = xRegion() ;  and bnd.dim == -1
  if (!bnd.dim()==-1)  ls_vitesses_normale.restrictTo(bnd);
  double my_CFL_timestep = std::fabs(xCFL::getDt(ls_vitesses_normale));
  increment_time(my_CFL_timestep);
  current_type_of_increment = " CFL imposed" ;
  current_convergence=false;
}



bool xpTimePilot::increment_with_CFL_Condition(xLevelSet ls_vitesses_normale,  double specified_dt, xRegion bnd)
{
  // if no bnd was given as argument then bnd = xRegion() ;  and bnd.dim == -1
  if (!bnd.dim()==-1)  ls_vitesses_normale.restrictTo(bnd);
  bool CFL_impose_timestep;
  double my_CFL_timestep = std::fabs(xCFL::getDt(ls_vitesses_normale));

    if (my_CFL_timestep < specified_dt)
      {
	// Then the CFL condition leads
	increment_time(my_CFL_timestep);
	CFL_impose_timestep = true;
	current_type_of_increment = " CFL imposed" ;  
      }
    else
      {
	increment_time(specified_dt);
	CFL_impose_timestep = false;
      }
    return CFL_impose_timestep;
}

