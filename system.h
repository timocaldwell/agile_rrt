//
//  system.h
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//
//  Within is the namespace sys which includes all integrations and differential equation solving for the dynamic system and equations relevant to Agile RRT and trajectory exploration.
//  Included:
//     Forward integration
//     Backward integration
//     Forward integration with constraint checking
//     Class for solving Ricatti Equations
//     Class for solving for reachability Gramians
//     Class for solving linear systems (time-varying and time-invariant)
//     Class for solving for state-transition matrix
//  All classes can be interfaced through method functions.

#ifndef __Agile_RRT__system__
#define __Agile_RRT__system__

#include "global.h"
#include "interp.h"
#include "pendcart_3link.h"

namespace sys
{
// Implemented here as a template.
// Uses odeint to integrate the system sys forward in time from initial state x0 over the time interval [0,tt_h] where tt_h>0. The timing and results of each step is returned through the pointers tt_vec and xx_vec.
// Returns true if integration is a success.
// Example integration of free dynamics from the zero state for 1.0 seconds:
//     vector<double> x0(NS,0.0)                                     // Initialize state to 0 vector
//     vector<double> *tt_vec, *xx_vec;                              // For passing results
//     PendCart pendcart_sys;                                        // Create a pendulum on a cart system object
//     bool is_success = IntegrateForward(pendcart_sys, x0, 1.0, tt_vec, xx_vec);
template <typename TT>
bool IntegrateForward(const TT & sys, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  ode_state_type xx = x0;
  double tt = 0;
  double dt = kINTEGRATION_DT_START;
  int step_num = 0;
  stepper_type stepper;
  
  tt_vec->clear();
  xx_vec->clear();
  // initial time and state.
  tt_vec->push_back(0.0);
  xx_vec->push_back(x0);
  odeint::controlled_step_result step_result;
  while (tt < tt_h && step_num < kMAX_NUM_INTEGRATION_STEPS) {
    // Shrink dt if the next step steps past the time horizon tt_h.
    if (tt_h - tt < dt)
      dt = tt_h - tt;

    step_result = stepper.try_step(sys , xx , tt , dt);
    if(step_result == 0) {                                        // successful step execution
      tt_vec->push_back(tt);
      xx_vec->push_back(xx);
    }
    step_num++;
  }
  if (tt >= tt_h)
    return true;
  return false;
}

// Implemented here as a template.
// Similar to IntegrateForward() except the system sys is integrated backward in time from final state x1 over the time interval [0,tt_h] where tt_h>0.
// The time-varying Ricatti equation and the slot 2 state-transition matrix differential equation are defined backward in time and solved using this function.
// Since odeint only integrates forward in time, time is reparameterized by ss where tt = tt_h-ss. Due to this reparameterization, the differential equation in sys must be negative of what it should be. To elaborate, if the equations are normally:
//       dxxdtt = ff(xx, tt),   xx(tt_h) = x1
//  then, the differential equation in sys must be transformed to gg(xx, ss) := -ff(xx, tt):
//       dxxdss = gg(xx, ss),   xx(0) = x1
template <typename TT>
bool IntegrateBackward(const TT & sys, const ode_state_type & x1, double tt_h, vector<double> * tts, vector<ode_state_type> * xxs) {
  ode_state_type xx = x1;
  stepper_type stepper;
  double ss = 0;
  double ds = kINTEGRATION_DT_START;
  int step_num = 0;
  tts->clear();
  xxs->clear();
  tts->push_back(tt_h);
  xxs->push_back(xx);
  odeint::controlled_step_result step_result;
  
  //integrates forward in time so setting ss = tt_h - tt
  while (ss < tt_h && step_num < kMAX_NUM_INTEGRATION_STEPS) {
    // Shrink dt if the next step steps past the time horizon tt_h.
    if (tt_h - ss < ds)
      ds = tt_h - ss;
    
    step_result = stepper.try_step(sys , xx , ss , ds);
    if(step_result == 0){
      tts->push_back(tt_h - ss);
      xxs->push_back(xx);
    }
    step_num++;
  }
  // If integration is successful, reverse to switch the parameterization back from ss to tt.
  if (ss >= tt_h) {
    reverse(tts->begin(), tts->end());
    reverse(xxs->begin(), xxs->end());
    return true;
  }
  return false;
}

// Solves dxxdtt = ff(xx, 0), xx(0) = x0
bool IntegrateFreeDynamics(const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
// Solves dxxdtt = ff(xx, uu_ff_interp), xx(0) = x0
bool IntegrateFeedForwardDynamics(const InterpVector & uu_ff_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
// Solves dxxdtt = ff(xx, KKs_interp*(xx_ref-xx)), xx(0) = x0
bool IntegratePointTrackingDynamics(const InterpVector & KKs_interp, const ode_state_type & xx_ref, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
// Solves dxxdtt = ff(xx, uu_ff_interp + KKs_interp*(xx_ref_interp-xx)), xx(0) = x0
bool IntegrateTrajectoryTrackingDynamics(const InterpVector & uu_ff_interp, const InterpVector & KKs_interp, const InterpVector & xx_ref_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);

// Similar to IntegrateForward() for just the PendCart system except it additionally:
//    1. checks at each step that the constraints specified in sys are satisfied and if not returns false.
//    2. passes the uu results to uu_vec. When PendCart is initialized, the location of uu computed for each step must be passed by uu_loc.
bool IntegrateDynamicSystemWithConstraints(const ode_state_type & uu_loc, const PendCart & sys, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec);
// Solves dxxdtt = ff(xx, uu_ff_interp + KKs_interp*(xx_ref_interp-xx)), xx(0) = x0
bool IntegrateTrajectoryTrackingDynamicsWithConstraints(const InterpVector & uu_ff_interp, const InterpVector & KKs_interp, const InterpVector & xx_ref_interp, const ode_state_type & x0, const constraints_struct & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec);
// XXXXXXXXXXX
bool IntegrateTrajectoryLinearSteeringProjectionWithConstraints(const InterpVector & xxzero, const InterpVector & BB, const InterpVector & KKlin, const InterpVector & KKproj, const InterpVector & WWK, const InterpVector & Phi, const NSx1_type & eta, const NIxNI_type & RRinv, const ode_state_type & x0, const constraints_struct & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec);
bool IntegrateTrajectoryLinearSteeringProjectionWithConstraintsWithCost(const NSx1_type & xx_samp, const InterpVector & xxzero, const InterpVector & BB, const InterpVector & KKlin, const InterpVector & KKproj, const InterpVector & WWK, const InterpVector & Phi, const NSx1_type & eta, const NSxNS_type & QQ, const NIxNI_type & RR, const NSxNS_type & P1, const NIxNI_type & RRinv, const ode_state_type & x0, const constraints_struct & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec, double * JJ);

class RicattiEq {
 public:
  RicattiEq(const NSxNS_type & AA, const NSxNI_type & BB, const NSxNS_type & QQ, const NIxNI_type & RRinv);
  RicattiEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NSxNS_type & QQ, const NIxNI_type & RRinv, double t_h);
  void operator()(const ode_state_type &PP_vec, ode_state_type &dPP_vec , const double ss);
 private:
  const NSxNS_type * AA_LTI_;
  const NSxNI_type * BB_LTI_;
  const NSxNS_type * QQ_;
  const NIxNI_type * RRinv_;
  int type_;
  const InterpVector * AAs_interp_;
  const InterpVector * BBs_interp_;
  double t_h_; // for backward integration
};
bool IntegrateLTIRicatti(const NSxNS_type & AA, const NSxNI_type & BB, const NSxNS_type & QQ, const NIxNI_type & RRinv, const ode_state_type & P1, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * PP_vec);
bool IntegrateLTVRicatti(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NSxNS_type & QQ, const NIxNI_type & RRinv, const ode_state_type & P1, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * PP_vec);

class ReachabilityEq {
 public:
  ReachabilityEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NIxNI_type & RRinv);
  ReachabilityEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const NIxNI_type & RRinv);
  ReachabilityEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const InterpVector & WWKs_interp, const NIxNI_type & RR, const NIxNI_type & RRinv);
  void operator()( const ode_state_type &WW_vec , ode_state_type &dWW_vec , const double tt );
 private:
  const InterpVector * AAs_interp_;
  const InterpVector * BBs_interp_;
  const InterpVector * KKs_interp_;
  const InterpVector * WWKs_interp_; // also for grammiantype == SSK_GRAMMIANTYPE
  const NIxNI_type * RR_;
  const NIxNI_type * RRinv_;
  int gramiantype_;
};
bool IntegrateWW(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * WW_vec);
bool IntegrateWWK(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * WWK_vec);
bool IntegrateSSK(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const InterpVector & WWKs_interp, const NIxNI_type & RR, const NIxNI_type & RRinv, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * SS_vec);

// dxx = AA*xx + BB*uu
class LinearEq {
 private:
  const NSxNS_type * AA_LTI_;
  const NSxNI_type * BB_LTI_;
  const NIxNI_type * RRinv_;
  vector<double> uu_vec_;
  int type_, controltype_;
  const InterpVector * AAs_interp_;
  const InterpVector * BBs_interp_;
  const InterpVector * uus_ff_interp_;
  const InterpVector * Phis_interp_;
  const InterpVector * KKs_interp_;
  const NSx1_type * eta_star_;

 public:
  LinearEq(const NSxNS_type & AA, const NSxNI_type & BB);
  LinearEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp);
  LinearEq(const NSxNS_type & AA, const NSxNI_type & BB, const InterpVector & uus_ff_interp);
  LinearEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & uus_ff_interp);
  LinearEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NIxNI_type & RRinv, const InterpVector & Phis_interp, const NSx1_type & eta_star);
  LinearEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const NIxNI_type & RRinv, const InterpVector & Phis_interp, const NSx1_type & eta_star);
  void operator()(const ode_state_type & xx_vec, ode_state_type & dxx_vec, const double tt);
};
bool IntegrateLTIFreeDynamics(const NSxNS_type & AA, const NSxNI_type & BB, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool IntegrateLTIFeedForwardDynamics(const NSxNS_type & AA, const NSxNI_type & BB, const InterpVector & uu_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool IntegrateLTVFreeDynamics(const InterpVector & AA_interp, const InterpVector & BB_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool IntegrateLTVFeedForwardDynamics(const InterpVector & AA_interp, const InterpVector & BB_interp, const InterpVector & uu_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool IntegrateLTVMinControlOpenLoopDynamics(const InterpVector & AA_interp, const InterpVector & BB_interp, const NIxNI_type & RRinv, const InterpVector & Phis_interp, const NSx1_type & eta_star, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool integrateLTVMinControlClosedLoopDynamics(const InterpVector & AA_interp, const InterpVector & BB_interp, const InterpVector & KK_interp, const NIxNI_type & RRinv, const InterpVector & Phis_interp, const NSx1_type & eta_star, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);

class STMslot2Eq { // dPhi(tt_h, t) = -Phi(tt_h, t) A(t)
 public:
  STMslot2Eq (const NSxNS_type & AA);
  STMslot2Eq (const NSxNS_type & AA, const NSxNI_type & BB, const NIxNS_type & KK);
  STMslot2Eq (const InterpVector & AAs_interp, double t_h);
  STMslot2Eq (const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, double t_h);
  void operator()( const ode_state_type &Phi_vec , ode_state_type &dPhi_vec , const double ss );
 private:
  const NSxNS_type * AA_LTI_;
  const NSxNI_type * BB_LTI_;
  const NIxNS_type * KK_LTI_;
  int type_;
  const InterpVector * AAs_interp_;
  const InterpVector * BBs_interp_;
  const InterpVector * KKs_interp_;
  double t_h_;  // For backward integration
};

bool IntegrateLTIOpenLoopSTMSlot2(const NSxNS_type & AA, const NSxNI_type & BB, const NIxNS_type & KK, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec);
bool IntegrateLTIClosedLoopSTMSlot2(const NSxNS_type & AA, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec);
bool IntegrateLTVOpenLoopSTMSlot2(const InterpVector & AAs_interp, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec);
bool IntegrateLTVClosedLoopSTMSlot2(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec);
} // namespace sys

#endif // __Agile_RRT__system__