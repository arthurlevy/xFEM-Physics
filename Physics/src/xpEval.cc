/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xpEval.h"
using namespace AOMD;
using namespace Trellis_Util;


  void ElementLengthMaterialVariablesVisitor_c::Visit(xValue<boost::shared_ptr<xTensors> >& v) 
    { 
      const bool debug = false;
      boost::shared_ptr<xTensors> tensors = v.getVal();
      if (debug) std::cout << "before setting variables\n";
      if (debug) v.print(std::cout);

      double volume = geo_integ->getMeasure();
      double h = pow(volume,power);
      //      boost::mpl::identity<double> L=h;
      tensors->get(variable, boost::mpl::identity<double>())=h; 

      if (debug) std::cout << "after setting variables\n";
      if (debug) v.print(std::cout);
    }


//--------------------------------------------------------------------------------------------------------------------------------

void xEvalLevelSetAtDisplacedPoint::operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, result_type& res) const
{

  bool found = false;
  const bool debug = false;

  //get current point
  mPoint mcurrent_point(geo_integ->getXYZ());
  if(debug) cout << "Current point" <<  mcurrent_point(0) << ", " << mcurrent_point(1) << ", " << mcurrent_point(2) << endl;

  //get displacement at current point
  xVector displacement  ;
  displ_field.getVal<xVector>(geo_appro, geo_integ, displacement);//velocity vector
  displacement = displacement * dt ; //displacement vector
  mPoint mdisplacement( displacement(0), displacement(1) , displacement(2) );
  if(debug) cout << "Displacement" <<  displacement(0) << ", " << displacement(1) << ", " << displacement(2) << endl;
 
  // get displaced_point
  mPoint   mdisplaced_point = mcurrent_point + mdisplacement ;
 

  //get the current element
  mEntity* e_curr = geo_appro->getEntity();
  int spaceDim = e_curr->getLevel();
  if(debug) cout << "Dimension : " << spaceDim << endl ;  

  // declares an element that may contain the displaced point
  mEntity* element_containing_displaced_point = geo_appro->getEntity();

  if (debug) cout << e_curr->getId() << endl ;
 
  //LOOKING FOR ELEMENT CONTAINING DISPLACED POINT
  if (geo_appro->PointInElement(mdisplaced_point))
    {
      element_containing_displaced_point = e_curr; // current element contains displaced point
      found =true;
    }
  else  //try to find it in adjacent elements
    {
      
      //set<mEntity> adjacent_set;
      //getAdjacencySet(2, e_curr, adjacent_set);
      
      for ( AOMD::mAdjacencyContainer::iter it = e_curr->begin(0);  it != e_curr->end(0); ++it)//every vertex of e_curr
	{
	  mVertex *v = (mVertex*) *it;
	  //How many elements are linked to this vertex?
	  int number_of_adjacent_element = v->size(spaceDim);
	  if (number_of_adjacent_element==0 )
	    {
	      cout << "The vertex number " << v->getId() <<  " has no element linked to it" << endl;
	    }
	  //iterates on all those adjacent elements
	  for (int i = 0 ; i!=number_of_adjacent_element; ++i)
	    {
	    mEntity *adjacent_element=v->get(spaceDim,i); // get ith element of size SpaceDim linked to v
	    xGeomElem geo(adjacent_element);
	    if (geo.PointInElement(mdisplaced_point))// if this adjacent element contains the displaced point.
	      {
		element_containing_displaced_point = adjacent_element ;
		found = true;
		break;
	      }
	    }
	}  
    }
  //If still not found
  if ( !found )
    {
      /*
      //LOCATE ELEMENT USING xRegularGrid : THIS IS LONG
      xMesh *msh(LS.getSupport().getMesh());
      std::set<mEntity*> set_elem_containing_disp_point;
      msh->locateElement(mdisplaced_point, set_elem_containing_disp_point);
      if (set_elem_containing_disp_point.size()>0)
	{
	  element_containing_displaced_point = *(set_elem_containing_disp_point.begin());
	  found = true;
	}
      else
	{
      */
	  /*
	  //  get element size approximation:
	  mPoint min, max;
	  geo_integ->GetBoundingBox(min, max);
	  xVector diag(min ,max);
	  
	  if ( displacement.mag() >  diag.mag()) cout << "Displacement is too big compared to element size, displaced point not found " <<endl;
	  */
	  if (debug) 
	    {
	      cout << "No adjacent of element number " << e_curr->getId() 
		   << " contains the displaced point " <<  mdisplaced_point(0) 
		   << ", " << mdisplaced_point(1) << ", " << mdisplaced_point(2) << " ." << endl;
	      /*
	      if ( displacement.mag() >  diag.mag())//mag returns the norm of an xVector
		{
		  cout << "    Displacement may be too big compared to element size" << endl ;
		}
	      */
	      cout << "    It may be outside mesh"<< endl;
	      
	    }
	  
	  
	  //returns an extrapolation of  the value of the levelset at the current point
	  res= LS.getVal(e_curr, geo_appro->getUVW()) + displacement*LS.getGrad(e_curr, geo_appro->getUVW()) ;
	  
	  
	
    }
 
  else //if found    
    {
      xGeomElem geo_appro_displaced( element_containing_displaced_point);
      geo_appro_displaced.setUVWForXYZ(mdisplaced_point);

      //get relative coordinate of that displaced point in that element
      mPoint displaced_point_uvw = geo_appro_displaced.getUVW();
      
      //return the value of the levelset
      res = LS.getVal( element_containing_displaced_point, displaced_point_uvw);
    }
}


/*
void xEvalLevelSetAtDisplacedPoint::getAdjacencySet(int nb_layer, mEntity &e_curr, set<mEntity> &res)
{
  
  res.insert(e_curr);
  set<mEntity> set_of_next_adjacents_to_add;
  for (int layer_iterator=1;layer_iterator!=nb_layer; layer_iterator++)
    {
      set_of_next_adjacents_to_add.clear();
      for ( set<mEntity>::iterator iter = res.begin();//pour chaque element deja dans le set
		iter !=  res.end() ;
		iter ++)
	{
	  for ( AOJ.GetMatrix(i,j)MD::mAdjacencyContainer::iter adj = iter->begin(spaceDim);  it != iter->end(spaceDim); ++it)// on cherche ses voisins
	    {
	      set_of_next_adjacents_to_add.insert(iter);// que l'on ajoute au set next...
	    }
	}
      
      
       for ( set<mEntity>::iterator iter = set_of_next_adjacents_to_add.begin();
		iter !=  set_of_next_adjacents_to_add.end() ;
		iter ++)
	 {
	   res.insert(iter);// on ajoute finalement ce next... au set "res".
	 }
    }
  return;
}
	  
*/

//--------------------------------------------------------------------------------------------------------------------------------

 void UpdateMaterialVariablesVisitorPartial_c::Visit(xValue<boost::shared_ptr<xTensors> >& v) 
    { 
      const bool debug = false;
      if (debug) std::cout << "before update variables\n";
      if (debug) v.print(std::cout);

      OldAndCurrent_c::current();
      boost::shared_ptr<xTensors> curr = v.getVal();
      OldAndCurrent_c::old();
      boost::shared_ptr<xTensors> old = v.getVal();
      curr_mat->setOldVariables(old);
      curr_mat->setCurrentVariables(curr);
      //      curr_mat->computeCurrentState(phys_to_update); //this is to be added in xMaterial
      curr_mat->computeCurrentState();

      //on se remet � current
      OldAndCurrent_c::current();
      if (debug) std::cout << "after update variables\n";
      if (debug) v.print(std::cout);
      
    }
