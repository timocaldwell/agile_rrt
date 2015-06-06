//
//  system.h
//  Agile_RRT
//
//  Created by Timothy Caldwell on 9/25/14.
//  Copyright (c) 2014 Timothy Caldwell. All rights reserved.
//

#ifndef __Agile_RRT__system__
#define __Agile_RRT__system__

#include "global.h"
#include "interp.h"
#include "pendcart_3link.h"

namespace sys
{
//Implemented here because templates are stupid
template <typename TT>
bool IntegrateForward(const TT & sys, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  ode_state_type xx = x0;
  double tt = 0;
  double dt = kINTEGRATION_DT_START;
  int step_num = 0;
  stepper_type stepper;
  
  tt_vec->clear();
  xx_vec->clear();
  tt_vec->push_back(0.0);
  xx_vec->push_back(x0);
  
  odeint::controlled_step_result step_result;
  
  while (tt < tt_h && step_num < kMAX_NUM_INTEGRATION_STEPS) {
    if (tt_h - tt < dt)
      dt = tt_h - tt;
    step_result = stepper.try_step( sys , xx , tt , dt );
    if(step_result == 0){
      tt_vec->push_back(tt);
      xx_vec->push_back(xx);
    }
    step_num++;
  }
  if (tt >= tt_h)
    return true;
  return false;
}

//Implemented here because templates are stupid
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
    if (tt_h - ss < ds)
      ds = tt_h - ss;
    step_result = stepper.try_step( sys , xx , ss , ds );
    if(step_result == 0){
      tts->push_back(tt_h - ss);
      xxs->push_back(xx);
    }
    step_num++;
  }
  if (ss >= tt_h) {
    reverse(tts->begin(), tts->end()); //Reverse to fix integrating forward in time.
    reverse(xxs->begin(), xxs->end());
    return true;
  }
  return false;
}

//Implemented here because templates are stupid
template <typename TT>
bool IntegrateForwardWithConstraints(const ode_state_type & uu_loc, const TT & sys, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec) {
  ode_state_type xx = x0;
  double tt = 0.0;
  double dt = 0.00001;
  int step_num = 0;
  stepper_type stepper;
  
  tt_vec->clear();
  xx_vec->clear();
  uu_vec->clear();
  tt_vec->push_back(0.0);
  xx_vec->push_back(x0);
  stepper.try_step(sys , xx , tt , dt); // calling to get uu at time 0 which populating uu_loc
  uu_vec->push_back(uu_loc);
  
  xx = x0;
  tt = 0.0;
  dt = kINTEGRATION_DT_START;
  odeint::controlled_step_result step_result;
  
  while (tt < tt_h && step_num < kMAX_NUM_INTEGRATION_STEPS) {
    if (tt_h - tt < dt)
      dt = tt_h - tt;
    step_result = stepper.try_step(sys, xx, tt, dt); // note stepper also populates uu_loc (shared memory)
    if ( !sys.IsConstraintSatisfied(xx_vec->back(), xx, uu_loc) )
      return false;
    if (step_result == 0){
      tt_vec->push_back(tt);
      xx_vec->push_back(xx);
      uu_vec->push_back(uu_loc);
    }
    step_num++;
  }
  if ( tt >= tt_h ){
    //check final state for constraints
    dt = .00001;
    step_result = stepper.try_step( sys , xx , tt , dt );

    if ( !sys.IsConstraintSatisfied(xx_vec->back(), xx, uu_loc) )
        return false;
    if( step_result == 0 )
        return true;
  }
  return false;
}
bool IntegrateFreeDynamics(const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool IntegrateFeedForwardDynamics(const InterpVector & uus_vec, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool IntegratePointTrackingDynamics(const InterpVector & KKs_vec, const ode_state_type & xx_ref, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool IntegrateTrajectoryTrackingDynamics(const InterpVector & uu_ff_interp, const InterpVector & KKs_vec, const InterpVector & xx_ref_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec);
bool IntegrateTrajectoryTrackingDynamicsWithConstraints(const InterpVector & uu_ff_interp, const InterpVector & KKs_vec, const InterpVector & xx_ref_interp, const ode_state_type & x0, const kin_constraints & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec);
bool IntegrateTrajectoryLinearSteeringProjectionWithConstraints(const InterpVector & xxzero, const InterpVector & BB, const InterpVector & KKlin, const InterpVector & KKproj, const InterpVector & WWK, const InterpVector & Phi, const NSx1_type & eta, const NIxNI_type & RRinv, const ode_state_type & x0, const kin_constraints & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec);
bool IntegrateTrajectoryLinearSteeringProjectionWithConstraintsWithCost(const NSx1_type & xx_samp, const InterpVector & xxzero, const InterpVector & BB, const InterpVector & KKlin, const InterpVector & KKproj, const InterpVector & WWK, const InterpVector & Phi, const NSx1_type & eta, const NSxNS_type & QQ, const NIxNI_type & RR, const NSxNS_type & P1, const NIxNI_type & RRinv, const ode_state_type & x0, const kin_constraints & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec, double * JJ);

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
bool IntegrateWWK(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * WW_vec);
bool IntegrateSSK(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const InterpVector & WWKs_interp, const NIxNI_type & RR, const NIxNI_type & RRinv, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * WW_vec);

// dxx = AA*xx + BB*uu
class LinearEq {
 private:
  const NSxNS_type * AA_LTI_;
  const NSxNI_type * BB_LTI_;
  const NIxNI_type * RRinv_;
  vector<double> uu_vec_;
  int type_, controltype_;
  const InterpVector *AAs_interp_;
  const InterpVector *BBs_interp_;
  const InterpVector *uus_ff_interp_;
  const InterpVector *Phis_interp_;
  const InterpVector *KKs_interp_;
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
  const NSxNS_type *AA_LTI_;
  const NSxNI_type *BB_LTI_;
  const NIxNS_type *KK_LTI_;
  int type_;
  const InterpVector *AAs_interp_;
  const InterpVector *BBs_interp_;
  const InterpVector *KKs_interp_;
  double t_h_;  // For backward integration
};

bool IntegrateLTIOpenLoopSTMSlot2(const NSxNS_type & AA, const NSxNI_type & BB, const NIxNS_type & KK, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec);
bool IntegrateLTIClosedLoopSTMSlot2(const NSxNS_type & AA, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec);
bool IntegrateLTVOpenLoopSTMSlot2(const InterpVector & AAs_interp, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec);
bool IntegrateLTVClosedLoopSTMSlot2(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec);
} // namespace sys

#endif // __Agile_RRT__system__