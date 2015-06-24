//
//  treevertex.h
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//
//  Within is class TreeVertex which initializes, preprocesses and stores all information needed for a single vertex of a Tree.

#ifndef __Agile_RRT__treevertex__
#define __Agile_RRT__treevertex__

#include "global.h"
#include "interp.h"
#include "system.h"
#include "pendcart_3link.h"

class TreeVertex {
 private:
  ode_state_type x0_;                                                   // Vertex state.
  InterpVector *xxedge_, *uuedge_;                                      // The state and control trajectories to connect the parent to the current state.
  InterpVector *xxzero_;                                                // The zero-control trajectory.
  InterpVector *AA_, *BB_;                                              // The time-varying linearization of the zero control trajectory or x0.
  InterpVector *KK_;                                                    // Time-varying feedback gain for Agile RRT steering.
  InterpVector *WWK_, *SSK_;                                            // Reachability Gramians for Agile RRT steering.
  // Quadratic weights for solving the LQR for calculating KK.
  const NSxNS_type QQlqr_, P1lqr_;
  const NIxNI_type RRlqr_, RRlqr_inv_;
  // Quadratic weights of the cost function.
  const NSxNS_type QQ_, P1_;
  const NIxNI_type RR_, RRinv_;
  NSx1_type xxgoal_;                                                    // The goal state for computing disttogoal.
  bool usezerotraj_;                                                    // true: Linearization about xxzero.
                                                                        // false: Linearization about x0.
  bool computedisttogoal_;                                              // true: Compute distance projection is to goal.
  bool save_edge_;                                                      // true: Store xxedge and uuedge.
  double t_h_max_upper_, t_h_bar_;                                      // Max time horizons.
  const constraints_struct constraints_;                                // Constraints on state during projection.
  
  // Temporary values that may be replaced each time a distance to xxsamp is computed for nearest neighbor.
  double t_h_cur_;                                                      // Current time horizon.
  double JJapprox_cur_;                                                 // Approximate cost from solving the LQR
  double JJproj_cur_;                                                   // Projected cost computed during projection.
  double disttogoal_cur_;                                               // Distance the projection is from xxgoal.
  // Temporary values for the Efficient linear fixed  time horizon steering.
  NSxNS_type Pth_cur_;
  NSx1_type xxzeroerror_cur_, xxsamp_cur_, eta_cur_;
  InterpVector *xxproj_cur_, *uuproj_cur_;                              // Projected state and control trajectories.

  // Conducts all preprocessing.
  void Initialize();
  
 public:
  TreeVertex(const ode_state_type & x0, bool usezerotraj, double t_h_max_upper, const InterpVector & xxedge, const InterpVector & uuedge, const constraints_struct & constraints, const NSxNS_type & QQ, const NIxNI_type & RR, const NSxNS_type & P1, const NSxNS_type & QQlqr, const NIxNI_type & RRlqr, const NSxNS_type & P1lqr, double t_h_bar, bool save_edge);
  ~TreeVertex();
  // For sample state xxsamp and time horizon t_h, conducts the inexact linear fixed time horizon steering. Returns the approximate cost and populates some temporary values (i.e. some private variables that end in _cur) for use by Projection().
  double JJLinearSteer(const NSx1_type & xxsamp, double t_h);
  // Returns the zero control state at time tt in xxzero_pt.
  void get_xxzero_pt(double tt, ode_state_type * xxzero_pt) {xxzero_->pt(tt, xxzero_pt);}
  // Specifies the goal state xxgoal and sets computedisttogoal to true. Otherwise computedisttogoal==false.
  void set_xxgoal(const NSx1_type & xxgoal);
  double get_disttogoal() {return disttogoal_cur_;}
  // Conducts the projection for inexact steering toward state xxsamp with time horizon t_h. Returns true if successful. Fails if collides with an obstacle or boundary specified in struct constraints.
  bool Projection(const NSx1_type & xxsamp, double t_h);
  // Prints projected trajectories xxproj_cur_ and uuproj_cur_ in a form easily readable by Mathematica. Must only be called after a successful execution of Projection() (i.e. Projection() returns true).
  void PrintTraj(double dt, const string & name, ostream * stream);
  // Prints edge trajectories from parent to this vertex xxedge_ and uuedge_ in a form easily readable by Mathematica.
  void PrintEdgeTraj(double dt, const string & name, ostream * stream);
  void get_x0(ode_state_type * x0) {x0 = &x0_;}
  // Creates a new vertex with vertex state at time t_h_cur_ of trajectory xxproj_cur_. Must only be called after a successful execution of Projection() (i.e. Projection() returns true).
  TreeVertex NewVertex();
};

#endif // __Agile_RRT__treevertex__
