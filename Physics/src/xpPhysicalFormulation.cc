/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xpPhysicalFormulation.h"
#include "xExportGmsh.h"
#include "xPhysSurf3.h"




xpPhysicalFormulation::xpPhysicalFormulation(xData * d):
  order_of_magnitude(1), 
  all(d->mesh), 
  domain(0), 
  thedata(d), 
  filter_integration (xAcceptAll()),
  domain_beg(d->mesh->begin(d->mesh->dim())),
  domain_end(d->mesh->end(d->mesh->dim())),
  spaceDim(d->mesh->dim()),
  usefulDofsDeclared(false),
  renumbered(false)
{;}



void xpPhysicalFormulation::initialize(   xEntityFilter  filter)
{
    declareApproximation();
    TreatmentOfEssEnv(filter, true);
    declareInternalVariables();
    if(renumbered) pair<int, int> res_renum = renumber();
    updateInternalVariables();
    return;
}

void  xpPhysicalFormulation::updateDomain(xLevelSet lsmat)
{
  cout  << name << " problem - updating domain with new Level set ";
  if (domain) delete domain;
  domain = new xPhysSurf(lsmat, xClassifyOn("matter"), xClassifyOn("air"));
  declareInternalVariables();
  updateInternalVariables();
  iso_zero= domain->getMesh_bnd();
  cout << "..... domain updated" << endl;
}

void  xpPhysicalFormulation::updateDomain(xLevelSet lsmat, xLevelSet lstable)
{
  cout  << name << " problem - updating domain with new Level sets ";
  if (domain) delete domain;
    
  domain = new xPhysSurf3(lsmat,lstable,  xClassifyOn("matter"),xClassifyOn("air"),xClassifyOn("middle"));
  
  iso_zero= domain->getMesh_bnd();
  
  declareInternalVariables();
  updateInternalVariables();

  cout << "..... domain updated" << endl;
}

void  xpPhysicalFormulation::updateDomain()
{
  cout  << name << " problem - updating domain with new Level set position ";
  domain->update();
  declareInternalVariables();
  updateInternalVariables();
  iso_zero= domain->getMesh_bnd();
  cout << "..... domain updated" << endl;
}



void  xpPhysicalFormulation::updateDomain(xLevelSet lsmat, xLevelSet lstable, xLevelSet lscompo)
{
  cout  << name << " problem - updating domain with new Level set ";
  if (domain) delete domain;

  domain = new xPhysSurf4(lsmat,lstable,lscompo,  xClassifyOn("matter"),xClassifyOn("air"),xClassifyOn("middle"), true, true );

  iso_zero= domain->getMesh_bnd();
  declareInternalVariables();
  updateInternalVariables();

  cout << "..... domain updated" << endl;

  xExportGmshAscii  pexport;
  if (Multphys_time_pilot->export_step())
    {
      pexport.openFile("LS" + Multphys_time_pilot->step2string());
      Export(lsmat,pexport ,"LS_" + Multphys_time_pilot->step2string() );
      Export(lscompo,pexport ,"COMPO_" + Multphys_time_pilot->step2string() );
    }   
  pexport.closeFile();

}

