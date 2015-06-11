//
//  treevertex.h
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//
// Comments for Treevertex.h coming soon.

#ifndef __Agile_RRT__treevertex__
#define __Agile_RRT__treevertex__

#include "global.h"
#include "interp.h"
#include "system.h"
#include "pendcart_3link.h"

class TreeVertex {
 private:
  ode_state_type x0_;
  InterpVector *xxedge_, *uuedge_, *xxzero_, *AA_, *BB_, *KK_, *WWK_, *SSK_;
  const NSxNS_type QQlqr_, P1lqr_, QQ_, P1_;
  const NIxNI_type RR_, RRinv_, RRlqr_, RRlqr_inv_;
  bool usezerotraj_, computedisttogoal_, save_edge_;
  double t_h_max_upper_, t_h_bar_;

  const constraints_struct constraints_;
  
  double t_h_cur_, JJapprox_cur_, JJproj_cur_, disttogoal_cur_;
  NSxNS_type Pth_cur_;
  NSx1_type xxzeroerror_cur_, xxsamp_cur_, eta_cur_;
  InterpVector *xxproj_cur_, *uuproj_cur_;
  NSx1_type xxgoal_;
  
  void Initialize();
  
 public:
  TreeVertex(const ode_state_type & x0, bool usezerotraj, double t_h_max_upper, const InterpVector & xxedge, const InterpVector & uuedge, const constraints_struct & constraints, const NSxNS_type & QQ, const NIxNI_type & RR, const NSxNS_type & P1, const NSxNS_type & QQlqr, const NIxNI_type & RRlqr, const NSxNS_type & P1lqr, double t_h_bar, bool save_edge);
  ~TreeVertex();
  double JJLinearSteer(const NSx1_type & xxsamp, double t_h);
  void get_xxzero_pt(double tt, ode_state_type * xxzero_pt) {xxzero_->pt(tt, xxzero_pt);}
  void set_xxgoal(const NSx1_type & xxgoal);
  double get_disttogoal() {return disttogoal_cur_;}
  bool Projection(const NSx1_type & xxsamp, double t_h);
  void PrintTraj(double dt, const string & name, ostream * stream);
  void PrintEdgeTraj(double dt, const string & name, ostream * stream);
  TreeVertex NewVertex();
};

#endif // __Agile_RRT__treevertex__
