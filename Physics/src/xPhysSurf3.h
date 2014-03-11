/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef ___XPHYSURF3__H
#define ___XPHYSURF3__H

#include "xPhysSurf.h"

using namespace xfem;


/// This class allows to define a three materials limited by two interfaces (two levelsets)
class xPhysSurf3 : public xPhysSurf
 {
  public:
   xPhysSurf3( xLevelSet& ls1, xLevelSet& ls2_,
   				const xEntityToEntity& in, const xEntityToEntity& out, const xEntityToEntity& middle,
   				bool fit=true, double fittol_=1.e-2,
   				bool keep_old_partition = false, bool recursive=false )
     : xPhysSurf(ls1, in, out,  fit, fittol_,keep_old_partition,recursive), ls2(ls2_),classify_middle(middle), mesh_bnd2(0)
     {
     	construct_(fit, keep_old_partition, recursive);
     }
     

   
   ~xPhysSurf3() { delete mesh_bnd2; }
 
   xMesh* getMesh_bnd2(){return mesh_bnd2;}
   virtual void  update(bool fit=true, bool keep_old_partition_flag = false);
   xEntityToEntity &getClassifyerMiddle(){return classify_middle;}

 private:
    void construct_(bool fit=true, bool keep_old_partition_flag = false, bool recursive=false);
    
 protected:
   xLevelSet& ls2;
   xEntityToEntity classify_middle;
   xMesh * mesh_bnd2;

 };
 
 
 /// This class allows to define a three materials limited by two interfaces (three levelsets) where the two first levelsets define the same interface (matter/air) but allow to describe more complex geometry (such as corner)
 /*! for instance if two same identical material come in contact, we would like to keep the two levelsets but they define a similar interface (matter/air) \see simulation of ultrasonic welding*/

 class xPhysSurf4 : public xPhysSurf3
 {
  public:
   xPhysSurf4( xLevelSet& ls1, xLevelSet& ls2_, xLevelSet& ls3_,
   				const xEntityToEntity& in, const xEntityToEntity& out, const xEntityToEntity& middle,
   				bool fit=true, double fittol_=1.e-2,
   				bool keep_old_partition = false, bool recursive=false )
     : xPhysSurf3(ls1, ls2_, in, out, middle, fit, fittol_,keep_old_partition,recursive), ls3(ls3_),mesh_bnd3(0)
     {
     	construct_(fit, keep_old_partition, recursive);
     }
     

   
   ~xPhysSurf4() { delete mesh_bnd3; }
 
   xMesh * getMesh_bnd3(){return mesh_bnd3;}
   virtual void  update(bool fit=true, bool keep_old_partition_flag = false);

 private:
    void construct_(bool fit=true, bool keep_old_partition_flag = false, bool recursive=false);
    
 protected:
   xLevelSet& ls3;
   xMesh * mesh_bnd3;

 };

 
 
 
 #endif
