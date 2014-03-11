/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xpAppli.h"
//#include "time.h"

void xpAppli::info()
{
  cout << "Application: " << name << endl;
  cout << "Contains Models:\n";
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); iter!=ListOfModels.end(); ++iter)
    {
      cout << "\t\t || "<< (*iter)->get_name() << endl;
    }
}

void xpAppli::exportFields(int details, const std::string& extension, bool binary, bool sorted) {
  cout << "Export de la solution de : " ;
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); iter!=ListOfModels.end(); ++iter)
    {
      cout <<  (*iter)->get_name() << " , ";
      (*iter)->exportFields(details, extension, binary, sorted);
    }
  cout << endl << endl ;
}

void xpAppli::init_models(xEntityFilter filter) {
  const bool debug=true;
  cout << "Initialisation of the list of models" << endl;
  transfer_timepilot();
  // initializing the "time_history " file
  fstream filestr;
  ofstream toto("time_history.txt", ios::app);
  toto << "=======================================================" << endl ;
  toto << "        Initialisation of the new application. " << endl ;
  toto << "                       +++++++                  " << endl;

  // append time and date.
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  //  toto <<"                " << asctime(timeinfo) << endl;
  toto << endl<< endl;

  toto << " The application contains : " << endl;
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); iter!=ListOfModels.end(); ++iter)
    {
      if (debug) cout <<  (*iter)->get_name() << endl;
      (*iter)->declareApproximation();
      (*iter)->TreatmentOfEssEnv(filter, true);
      toto << (*iter)->get_name() << endl;
    }
  
  
 //append the data file.
  toto << "=================== MATERIAL DATA USED ==================" << endl;
  ifstream tata("data/matter.mat",ifstream::in);
  toto << tata.rdbuf();
  tata.close();
  ifstream tata2("data/air.mat",ifstream::in);
  toto << tata2.rdbuf();
  tata2.close();
  toto << "=============== MESH & BOUNDARY CONDITIONS  =============" << endl;
  ifstream tutu("data/main.dat",ifstream::in);
  toto << tutu.rdbuf();
  tutu.close();
  toto << "=========================================================" << endl ;
  
  transfer_varmanager();
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); 
  			iter!=ListOfModels.end(); 
  			++iter)
  	{
   	  (*iter)->declareInternalVariables();
  	}
}

void xpAppli::init_models_with_approx_already_declared(xEntityFilter filter) {
  const bool debug=true;
  cout << "Initialisation of the list of models" << endl;
  transfer_timepilot();
  // initializing the "time_history " file
  fstream filestr;
  ofstream toto("time_history.txt", ios::app);
  toto << "=======================================================" << endl ;
  toto << "        Initialisation of the new application. " << endl ;
  toto << "                       +++++++                  " << endl;

  // append time and date.
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  toto <<"                " << asctime(timeinfo) << endl;
  toto << endl<< endl;

  toto << " The application contains : " << endl;
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); iter!=ListOfModels.end(); ++iter)
    {
      if (debug) cout <<  (*iter)->get_name() << endl;
      (*iter)->TreatmentOfEssEnv(filter, true);
      toto << (*iter)->get_name() << endl;
    }
  
  
 //append the data file.
  toto << "=================== MATERIAL DATA USED ==================" << endl;
  ifstream tata("data/matter.mat",ifstream::in);
  toto << tata.rdbuf();
  tata.close();
  ifstream tata2("data/air.mat",ifstream::in);
  toto << tata2.rdbuf();
  tata2.close();
  toto << "=============== MESH & BOUNDARY CONDITIONS  =============" << endl;
  ifstream tutu("data/main.dat",ifstream::in);
  toto << tutu.rdbuf();
  tutu.close();
  toto << "=========================================================" << endl ;
  
  transfer_varmanager();
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); 
  		iter!=ListOfModels.end(); 
  		++iter)
   {
   	  (*iter)->declareInternalVariables();
   }

}


void xpAppli::transfer_varmanager()
{
  const bool debug=true;
  cout << "models get the variable manager" << endl;
  int nbmodel = ListOfModels.size();
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); iter!=ListOfModels.end(); ++iter)
    {
      if (debug) cout <<  (*iter)->get_name() << endl;
      (*iter)->setVarManager(&VarMNGR);
    }
}

void xpAppli::transfer_timepilot()
{
  const bool debug=true;
  cout << "models get the time pilot" << endl;
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); iter!=ListOfModels.end(); ++iter)
    {
      if (debug) cout <<  (*iter)->get_name() << endl;
      (*iter)->setTimePilot(&Timepilot);
    }
}

void xpAppli::updateDomain(xLevelSet lsmat)
{
  const bool debug=false;
  ListOfModels[0]->updateDomain(lsmat);
  xRegion iso(ListOfModels[0]->get_iso_zero() );
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin();
  			 iter!=ListOfModels.end();
  			 ++iter)
    { 
      if (debug)  cout << "Appli tranfer the iso_zero to : " <<  (*iter)->get_name() << endl;
      (*iter)->set_iso_zero(iso);
      (*iter)->declareInternalVariables();
    }

}


void xpAppli::updateDomain(xLevelSet lsmat, xLevelSet tool)
{
  const bool debug=false;
  ListOfModels[0]->updateDomain(lsmat,tool);
  xRegion iso(ListOfModels[0]->get_iso_zero() );
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin();
  			 iter!=ListOfModels.end();
  			 ++iter)
    { 
      if (debug)  cout << "Appli tranfer the iso_zero to : " <<  (*iter)->get_name() << endl;
      (*iter)->set_iso_zero(iso);
      (*iter)->declareInternalVariables();
    }

}




void xpAppli::updateDomain(xLevelSet lsmat, xLevelSet tool, xLevelSet compo)
{
  const bool debug=false;
  if (debug)  cout << "Appli update the domain of model : " << endl;
  ListOfModels[0]->updateDomain(lsmat,tool, compo);
  xRegion iso(ListOfModels[0]->get_iso_zero() );
  for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); iter!=ListOfModels.end(); ++iter)
    { 
      if (debug)  cout << "Appli tranfer the iso_zero to : " <<  (*iter)->get_name() << endl;
      (*iter)->set_iso_zero(iso);
      (*iter)->declareInternalVariables();
    }
}



void xpAppli::declareApproximations()
  {
    bool debug=0;
    for (vector<xpPhysicalFormulation*>::iterator iter=ListOfModels.begin(); iter!=ListOfModels.end(); ++iter)
      {
	if (debug) cout <<  "declaring approximation of model :" << (*iter)->get_name() << endl;
	(*iter)->declareApproximation();
      }
  }
  
  
  
void xpAppli::add_model(xpPhysicalFormulation &amodel)
  {
    cout <<"Adding model: "<< amodel.get_name() << endl;
    ListOfModels.push_back(&amodel);
    return;
  }

bool xpAppli::solve()
  { 
    bool result=false;
    result =  Applisolver->solve(ListOfModels);
    Timepilot.setCurrentConvergence(result);
    return  result;
  }
