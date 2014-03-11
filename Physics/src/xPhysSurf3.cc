/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#include "xPhysSurf3.h"
#include "xLevelSetOperators.h"
#include "xEntityToEntity.h"

using namespace AOMD;
using namespace Trellis_Util;

void xPhysSurf3 :: update(bool fit, bool keep_old_partition_flag)
{
  construct_(fit, keep_old_partition_flag);
}

void xPhysSurf4 :: update(bool fit, bool keep_old_partition_flag)
{
  construct_(fit, keep_old_partition_flag);
}


void xPhysSurf3 ::construct_(bool fit, bool keep_old_partition_flag, bool recursive)
{
  if (fit) 
    {
      xFitToVertices fit(fittol);
      ls.accept(fit);
      ls2.accept(fit);

    }

  if (mesh_bnd) delete mesh_bnd;
  if (mesh_bnd2) delete mesh_bnd2;
  
  mesh_bnd = new xMesh;
  mesh_bnd2 = new xMesh;
  
  mesh->cutMesh(ls,mesh_bnd, xAcceptAll(),classify_in, classify_middle, true, keep_old_partition_flag, recursive); //keep_old_partition_flag  indicates that we keep the old partition.
   mesh->cutMesh(ls2,mesh_bnd2, xAcceptAll(),classify_middle, classify_out, true, true, false); //last  indicates that we keep the old partition.

   // createSupportInfo();
   int dimmesh = mesh->dim();
   for (xIter it  = mesh->begin(dimmesh);it != mesh->end(dimmesh); ++it){
     std::set<mEntity*>  partition;
     xMesh::getPartition((*it), partition );
     for ( std::set<mEntity*>::iterator   itsub =partition.begin() ; 
	   itsub != partition.end(); ++itsub){
       xGeomElem geo_sub(*itsub);
       xGeomElem geo(*it);
       mPoint cdgxyz = geo_sub.getCDGxyz();
       geo.setUVWForXYZ(cdgxyz);
       double val  = ls.getVal((*it) , geo.getUVW() ); 
       double val2 = ls2.getVal((*it) , geo.getUVW() ); 

       if (val < 0.)  classify_in(*itsub);
       else if (val2 < 0.)  classify_out(*itsub);
       else  classify_middle(*itsub);
     
     }
   
   }
   return;
} 




void xPhysSurf4 ::construct_(bool fit, bool keep_old_partition_flag, bool recursive)
{
  if (fit) 
    {
      xFitToVertices fit(fittol);
      ls.accept(fit);
      ls2.accept(fit);
      ls3.accept(fit);

    }

  if (mesh_bnd) delete mesh_bnd;
  if (mesh_bnd2) delete mesh_bnd2;
  if (mesh_bnd2) delete mesh_bnd3;

  
  mesh_bnd = new xMesh;
  mesh_bnd2 = new xMesh;
  mesh_bnd3 = new xMesh;

  
  mesh->cutMesh(ls,mesh_bnd, xAcceptAll(),classify_in, classify_middle, true, keep_old_partition_flag, recursive); //keep_old_partition_flag  indicates that we keep the old partition.
  mesh->cutMesh(ls2,mesh_bnd2, xAcceptAll(),classify_in, classify_out, true, true, recursive); //last  indicates that we keep the old partition.
  mesh->cutMesh(ls3,mesh_bnd3, xAcceptAll(),classify_middle, classify_in, true, true, recursive); //last  indicates that we keep the old partition.


  createSupportInfo();
   int dimmesh = mesh->dim();
   for (xIter it  = mesh->begin(dimmesh);it != mesh->end(dimmesh); ++it){
     std::set<mEntity*>  partition;
     xMesh::getPartition((*it), partition );
     for ( std::set<mEntity*>::iterator   itsub =partition.begin() ; 
	   itsub != partition.end(); ++itsub){
       xGeomElem geo_sub(*itsub);
       xGeomElem geo(*it);
       mPoint cdgxyz = geo_sub.getCDGxyz();
       geo.setUVWForXYZ(cdgxyz);
       double val  = ls.getVal((*it) , geo.getUVW() ); 
       double val2 = ls2.getVal((*it) , geo.getUVW() );
       double val3 = ls3.getVal((*it) , geo.getUVW() );

       if (val3>0)  classify_middle(*itsub);
       else if (val2 > 0. && val>0)  classify_out(*itsub);
       else  classify_in(*itsub);
     
     }
   
   }
   return;
}  
