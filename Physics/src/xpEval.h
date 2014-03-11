/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef _XPEVAL_H
#define _XPEVAL_H
  
#include "ValueOldAndCurrent.h"
#include "xField.h"
#include "MaterialCommand.h"
#include "xLevelSet.h"

using namespace xfem;

/// evaluates the time derivative of a field at any point (for an approximation)
/*! see class xField
 Warning : requires OldAndCurrent value type
*/
template <class UnaryOperator,class Field=xField> 
class  xEvalTimeDerivativeField : public xEval<typename UnaryOperator::result_type>  
{
public:
 typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

  xEvalTimeDerivativeField(const Field& f_, const double dt_) : f(f_),dt(dt_) 
  {  invdt = 1./dt;}
  xEvalTimeDerivativeField(const Field& f_, const UnaryOperator& _funct, const double dt_) : 
    f(f_), funct(_funct), dt(dt_)  
  {   invdt = 1./dt;  }  
  void operator()(const xGeomElem*  appro, const xGeomElem* integ, result_type& result) const
  {
      typename UnaryOperator::argument_type vcur;
      typename UnaryOperator::argument_type vold;
      OldAndCurrent_c::current();
      f.getVal(appro, integ, vcur); 
      OldAndCurrent_c::old();
      f.getVal(appro, integ, vold); 
      OldAndCurrent_c::current();
      result = (funct(vcur)-funct(vold))*invdt;
      //      cout << "valeur current: "<< vcur <<"   valeur old: "<< vold << "   derivee: "<< result <<endl;
    }
  
private:
  const Field& f;
  double dt, invdt;
  UnaryOperator funct; 
};


/// evaluates (1-Theta)*T[old] + Theta*T[current]
/*! of an xField T.
 This evaluation is useful for the Weighted residual time integration Scheme.
*/
template <class UnaryOperator,class Field=xField> 
class  xEvalThetaField : public xEval<typename UnaryOperator::result_type>  
{
public:
 typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

  xEvalThetaField(const Field& f_, double theta_) : f(f_), Theta(theta_)  {}
  xEvalThetaField(const Field& f_, const UnaryOperator& _funct, double theta_) : f(f_), funct(_funct), Theta(theta_)  {}  
  //non optimized version, an optimized version exists when the UnaryOperator is Identity
  void operator()(const xGeomElem*  appro, const xGeomElem* integ, result_type& result) const
  {
    typename UnaryOperator::argument_type vold;
    OldAndCurrent_c::old();
    f.getVal(appro, integ, vold);    
    typename UnaryOperator::argument_type vcurrent;
    OldAndCurrent_c::current();
    f.getVal(appro, integ, vcurrent);
    result = funct(vold)*(1-Theta)+funct(vcurrent)*Theta;
  }
private:
  const Field& f;
  UnaryOperator funct;
  double Theta;

};



/// evaluates the gradient of  (1-Theta)*T[old] + Theta*T[current]
/*! of an xField T.
 This evaluation is useful for the Weighted residual time integration Scheme.
*/
template <class UnaryOperator,class Field=xField> 
class  xEvalGradThetaField : public xEval<typename UnaryOperator::result_type>  
{
public:
  typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

  xEvalGradThetaField(const Field& f_, double theta_) : f(f_), Theta(theta_)  {}
  xEvalGradThetaField(const Field& f_, const UnaryOperator& _funct, double theta_) : f(f_), funct(_funct), Theta(theta_)  {}  
  //non optimized version, an optimized version exists when the UnaryOperator is Identity
  void operator()(const xGeomElem*  appro, const xGeomElem* integ, result_type& result) const
  {
    typename UnaryOperator::argument_type vold;
    OldAndCurrent_c::old();
    f.getGrad(appro, integ, vold);

    typename UnaryOperator::argument_type vcurrent;
    OldAndCurrent_c::current();
    f.getGrad(appro, integ, vcurrent);
    result = funct(vold)*(1-Theta)+funct(vcurrent)*Theta;
  }
 
private:
  const Field& f;
  UnaryOperator funct;
  double Theta;
};


//------------------------------------------------------------------------------
/// Returns the element size, usefull for the CFL condition on Levelset propagation.
class ElementLengthMaterialVariablesVisitor_c : public MaterialVariablesVisitor_c 
{
public:
  ElementLengthMaterialVariablesVisitor_c(const string& var, xRegion r) :  variable(var) 
  {
    switch(r.getMesh()->dim())
    {
    case 3 :  power = 1./3.;    break;
    case 2 :   power = 1./2.;   break;      
    }
  }
  
  void Visit(xValue<boost::shared_ptr<xTensors> >& v);
  
private:
  const string variable;
  double power;

};




//----------------/----------------/----------------/----------------
///Evaluates the value of the levelset LS at the displaced point X+U or X+Vdt. The sign of the operator () can therefore  evaluates wheather the displaced point penetrates a tool modelized by LS or stays outside.
class xEvalLevelSetAtDisplacedPoint : public xEval<double>
{
public: 
  xEvalLevelSetAtDisplacedPoint(const xLevelSet & _LS, const xField &_U) : LS(_LS), displ_field(_U), dt(1){}
  xEvalLevelSetAtDisplacedPoint(const xLevelSet & _LS, const xField &_V, double time_step) :  LS(_LS), displ_field(_V), dt(time_step){}

  void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, result_type& res) const;
  //  void getAdjacencySet(int nb_layer, mEntity &e_curr, set<mEntity> &res);
private:
  const xLevelSet& LS;
  const xField& displ_field;
  double dt;

};

///Evaluates the value of the levelset LS at the displaced point X+U or X+Vdt. And returns one if ther is penetration and 0 if not..
class xEvalIfPenetrationAtDisplacedPoint  : public  xEvalLevelSetAtDisplacedPoint
{
public: 
  xEvalIfPenetrationAtDisplacedPoint(const xLevelSet & _LS, const xField &_U) : 
    xEvalLevelSetAtDisplacedPoint( _LS,  _U) {}
  xEvalIfPenetrationAtDisplacedPoint(const xLevelSet & _LS, const xField &_V, double time_step) : 
    xEvalLevelSetAtDisplacedPoint( _LS, _V,  time_step) {}
  void operator()(xGeomElem* geo_appro, xGeomElem* geo_integ, result_type& res) const
  {
    xEvalLevelSetAtDisplacedPoint::operator()(geo_appro, geo_integ, res);
    if  (res < 0) {  res = 1;}
    else {res = 0;}
    return;
  }



};




class UpdateMaterialVariablesVisitorPartial_c :
public MaterialVariablesVisitor_c
{
 public:
  UpdateMaterialVariablesVisitorPartial_c(string _phys)  : phys_to_update(_phys) {}
 
  void Visit(xValue<boost::shared_ptr<xTensors> >& v) ;
  
 private:
  string phys_to_update;
  //  xThermoTrans* curr_mat;
};


///This is an xFilteredRegion that is usefull to define a P1+ interpolation.
struct bubbleFunction {
  bubbleFunction(int d):dim(d) {}
  bool operator()(AOMD::mEntity* e) {
 //    cout<<"LEVEL="<<e->getType()<<endl;
    if(e->getLevel()==dim) return true;
    else return false;

  }
  int dim;
};





#endif
