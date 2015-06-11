//
//  system.cpp
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//

#include "system.h"

namespace sys
{
bool IntegrateFreeDynamics(const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  PendCart pendcart_sys;
  return IntegrateForward(pendcart_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegrateFeedForwardDynamics(const InterpVector & uu_ff_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  PendCart pendcart_sys(uu_ff_interp);
  return IntegrateForward(pendcart_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegratePointTrackingDynamics(const InterpVector & KKs_interp, const ode_state_type & xx_ref, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  PendCart pendcart_sys(KKs_interp, xx_ref);
  return IntegrateForward(pendcart_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegrateTrajectoryTrackingDynamics(const InterpVector & uu_ff_interp, const InterpVector & KKs_interp, const InterpVector & xx_ref_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  PendCart pendcart_sys(uu_ff_interp, KKs_interp, xx_ref_interp);
  return IntegrateForward(pendcart_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegrateDynamicSystemWithConstraints(const ode_state_type & uu_loc, const PendCart & sys, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec) {
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
  stepper.try_step(sys , xx , tt , dt); // calling to get uu at time 0 to get uu(0) at uu_loc
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
bool IntegrateTrajectoryTrackingDynamicsWithConstraints(const InterpVector & uu_ff_interp, const InterpVector & KKs_interp, const InterpVector & xx_ref_interp, const ode_state_type & x0, const constraints_struct & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec) {
  ode_state_type uu_loc(NI);
  PendCart pendcart_sys(uu_ff_interp, KKs_interp, xx_ref_interp, constraints, &uu_loc);
  return IntegrateDynamicSystemWithConstraints(uu_loc, pendcart_sys, x0, tt_h, tt_vec, xx_vec, uu_vec);
}
bool IntegrateTrajectoryLinearSteeringProjectionWithConstraints(const InterpVector & xxzero, const InterpVector & BB, const InterpVector & KKlin, const InterpVector & KKproj, const InterpVector & WWK, const InterpVector & Phi, const NSx1_type & eta, const NIxNI_type & RRinv, const ode_state_type & x0, const constraints_struct & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec) {
  ode_state_type uu_loc(NI);
  PendCart pendcart_sys(xxzero, BB, KKlin, KKproj, WWK, Phi, eta, RRinv, constraints, &uu_loc);
  return IntegrateDynamicSystemWithConstraints(uu_loc, pendcart_sys, x0, tt_h, tt_vec, xx_vec, uu_vec);
}
bool IntegrateTrajectoryLinearSteeringProjectionWithConstraintsWithCost(const NSx1_type & xx_samp, const InterpVector & xxzero, const InterpVector & BB, const InterpVector & KKlin, const InterpVector & KKproj, const InterpVector & WWK, const InterpVector & Phi, const NSx1_type & eta, const NSxNS_type & QQ, const NIxNI_type & RR, const NSxNS_type & P1, const NIxNI_type & RRinv, const ode_state_type & x0, const constraints_struct & constraints, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec, vector<ode_state_type> * uu_vec, double * JJ) {

  ode_state_type x0J0 = x0;
  x0J0.push_back(0.0);
  
  ode_state_type uu_loc(NI);
  PendCart pendcart_sys(xxzero, BB, KKlin, KKproj, WWK, Phi, eta, QQ, RR, RRinv, constraints, &uu_loc);
  bool integration_success = IntegrateDynamicSystemWithConstraints(uu_loc, pendcart_sys, x0J0, tt_h, tt_vec, xx_vec, uu_vec);
  Eigen::Map<const NSx1_type> xx_mat(xx_vec->back().data(), NS, 1);
  *JJ = xx_vec->back().back() /* ell */ + 0.5*(xx_mat-xx_samp).dot(P1*(xx_mat-xx_samp)) /* m */;

  // remove int_0^t ell dtau from xx_vec
  for(vector<ode_state_type>::iterator it = xx_vec->begin(); it!=xx_vec->end(); ++it)
      it->pop_back();
  return integration_success;
}

RicattiEq::RicattiEq(const NSxNS_type & AA, const NSxNI_type & BB, const NSxNS_type & QQ, const NIxNI_type & RRinv)
    : AA_LTI_(&AA), BB_LTI_(&BB), QQ_(&QQ), RRinv_(&RRinv), type_(kLTI_TYPE) {}
RicattiEq::RicattiEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NSxNS_type & QQ, const NIxNI_type & RRinv, double t_h)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), QQ_(&QQ), RRinv_(&RRinv), t_h_(t_h), type_(kLTV_TYPE) {}
void RicattiEq::operator()(const ode_state_type &PP_vec , ode_state_type &dPP_vec , const double ss) {
  Eigen::Map<NSxNS_type> dPP(dPP_vec.data(), NS, NS); // mat_dPP and dPP share memory
  Eigen::Map<const NSxNS_type> PP(PP_vec.data(), NS, NS);
  if (type_ == kLTV_TYPE) {
    double tt = t_h_-ss;
    vector<double> AA_vec(NS*NS);
    AAs_interp_->pt(tt, &AA_vec);
    Eigen::Map<NSxNS_type> AA(AA_vec.data(), NS, NS);
    
    vector<double> BB_vec(NS*NI);
    BBs_interp_->pt(tt, &BB_vec);
    Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);
    
    //Integrator only integrates forward in time. Otherwise dPP = -dPP.
    dPP = AA.transpose()*PP + PP*AA - PP*BB*(*RRinv_)*BB.transpose()*PP + (*QQ_);
  }
  else
    //Integrator only integrates forward in time. Otherwise dPP = -dPP.
    dPP = (*AA_LTI_).transpose()*PP + PP*(*AA_LTI_) - PP*(*BB_LTI_)*(*RRinv_)*(*BB_LTI_).transpose()*PP + (*QQ_);
}
bool IntegrateLTIRicatti(const NSxNS_type & AA, const NSxNI_type & BB, const NSxNS_type & QQ, const NIxNI_type & RRinv, const ode_state_type & P1, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * PP_vec) {
  RicattiEq ric_sys(AA, BB, QQ, RRinv);
  return IntegrateBackward(ric_sys, P1, tt_h, tt_vec, PP_vec);
}
bool IntegrateLTVRicatti (const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NSxNS_type & QQ, const NIxNI_type & RRinv, const ode_state_type & P1, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * PP_vec) {
  RicattiEq ric_sys(AAs_interp, BBs_interp, QQ, RRinv, tt_h);
  return IntegrateBackward(ric_sys, P1, tt_h, tt_vec, PP_vec);
}

ReachabilityEq::ReachabilityEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NIxNI_type & RRinv)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), RRinv_(&RRinv), gramiantype_(kWW_GRAMIANTYPE) {}
ReachabilityEq::ReachabilityEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const NIxNI_type & RRinv)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), RRinv_(&RRinv), KKs_interp_(&KKs_interp), gramiantype_(kWWK_GRAMIANTYPE) {}
ReachabilityEq::ReachabilityEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const InterpVector & WWKs_interp, const NIxNI_type & RR, const NIxNI_type & RRinv)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), RR_(&RR), RRinv_(&RRinv), KKs_interp_(&KKs_interp), WWKs_interp_(&WWKs_interp), gramiantype_(kSSK_GRAMIANTYPE) {}
void ReachabilityEq::operator()(const ode_state_type &WW_vec , ode_state_type &dWW_vec , const double tt) {
//  clock_t t_start2 = clock();

  Eigen::Map<NSxNS_type> dWW(dWW_vec.data(), NS, NS);
  Eigen::Map<const NSxNS_type> WW(WW_vec.data(), NS, NS);

  vector<double> AA_vec(NS*NS);
  AAs_interp_->pt(tt, &AA_vec);
  Eigen::Map<NSxNS_type> AA(AA_vec.data(), NS, NS);
  
  vector<double> BB_vec(NS*NI);
  BBs_interp_->pt(tt, &BB_vec);
  Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);

  if (gramiantype_ == kWW_GRAMIANTYPE) {
    dWW = AA*WW + WW*AA.transpose() + BB*(*RRinv_)*BB.transpose();
  } else {
    vector<double> KK_vec(NI*NS);
    KKs_interp_->pt(tt, &KK_vec);
    Eigen::Map<NIxNS_type> KK(KK_vec.data(), NI, NS);
    
    NSxNS_type AAK = AA-BB*KK;
    
    if (gramiantype_ == kWWK_GRAMIANTYPE) {
      dWW = AAK*WW + WW*(AAK.transpose()) + BB*(*RRinv_)*BB.transpose();
    } else {                                                           //gramiantype = SSK_GRAMMIANTYPE
      vector<double> WWK_vec(NS*NS);
      WWKs_interp_->pt(tt, &WWK_vec);
      Eigen::Map<NSxNS_type> WWK(WWK_vec.data(), NS, NS);

      NIxNS_type CC = KK*WWK - (*RRinv_)*BB.transpose();
      
      dWW = AAK*WW + WW*(AAK.transpose()) + CC.transpose()*(*RR_)*CC;
    }
  }
}
bool IntegrateWW(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * WW_vec) {
  ReachabilityEq reach_sys(AAs_interp, BBs_interp, RRinv);
  ode_state_type WW0(NS*NS, 0.0);                                     // initializes WW0 to the zero matrix
  return IntegrateForward(reach_sys, WW0, tt_h, tt_vec, WW_vec);
}
bool IntegrateWWK(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * WWK_vec) {
  ReachabilityEq reach_sys(AAs_interp, BBs_interp, KKs_interp, RRinv);
  ode_state_type WWK0(NS*NS, 0.0);                                    // initializes WWK0 to the zero matrix
  return IntegrateForward(reach_sys, WWK0, tt_h, tt_vec, WWK_vec);
}
bool IntegrateSSK(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const InterpVector & WWKs_interp, const NIxNI_type & RR, const NIxNI_type & RRinv, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * SSK_vec) {
  ReachabilityEq reach_sys(AAs_interp, BBs_interp, KKs_interp, WWKs_interp, RR, RRinv);
  ode_state_type SSK0(NS*NS, 0.0);                                    // initializes SSK0 to the zero matrix
  return IntegrateForward(reach_sys, SSK0, tt_h, tt_vec, SSK_vec);
}

LinearEq::LinearEq(const NSxNS_type & AA, const NSxNI_type & BB)
    : AA_LTI_(&AA), BB_LTI_(&BB), type_(kLTI_TYPE), controltype_(kFREE_CONTROLTYPE), uu_vec_(NI,0.0) {}
LinearEq::LinearEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), uu_vec_(NI,0.0), type_(kLTV_TYPE), controltype_(kFREE_CONTROLTYPE) {}
LinearEq::LinearEq(const NSxNS_type & AA, const NSxNI_type & BB, const InterpVector & uus_ff_interp)
    : AA_LTI_(&AA), BB_LTI_(&BB), uus_ff_interp_(&uus_ff_interp), uu_vec_(NI), type_(kLTI_TYPE), controltype_(kFEEDFORWARD_CONTROLTYPE) {}
LinearEq::LinearEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & uus_ff_interp)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), uus_ff_interp_(&uus_ff_interp), uu_vec_(NI), type_(kLTV_TYPE), controltype_(kFEEDFORWARD_CONTROLTYPE) {}
LinearEq::LinearEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const NIxNI_type & RRinv, const InterpVector & Phis_interp, const NSx1_type & eta_star)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), RRinv_(&RRinv), Phis_interp_(&Phis_interp), eta_star_(&eta_star), uu_vec_(NI), type_(kLTV_TYPE), controltype_(kMINENERGY_OPENLOOP_CONTROLTYPE) {}
LinearEq::LinearEq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, const NIxNI_type & RRinv, const InterpVector & Phis_interp, const NSx1_type & eta_star)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), KKs_interp_(&KKs_interp), RRinv_(&RRinv), Phis_interp_(&Phis_interp), eta_star_(&eta_star), uu_vec_(NI), type_(kLTV_TYPE), controltype_(kMINENERGY_CLOSEDLOOP_CONTROLTYPE) {}
// dxx = AA*xx + BB*uu
void LinearEq::operator()( const ode_state_type &xx_vec , ode_state_type &dxx_vec , const double tt ) {
  Eigen::Map<NSx1_type> dxx(dxx_vec.data(), NS, 1);
  Eigen::Map<const NSx1_type> xx(xx_vec.data(), NS, 1);
  Eigen::Map<NIx1_type> uu(uu_vec_.data(), NI, 1);

  if (controltype_ == kFEEDFORWARD_CONTROLTYPE) {
    // *** feedforward:  u = u_ff *** //
    uus_ff_interp_->pt(tt, &uu_vec_);
  } else if (controltype_ == kMINENERGY_OPENLOOP_CONTROLTYPE || controltype_ == kMINENERGY_CLOSEDLOOP_CONTROLTYPE) {
    // *** MIN ENERGY OPEN LOOP:  u = -RRinv * BB^T * Phi^T * eta *** //
    vector<double> Phi_vec(NS*NS);
    Phis_interp_->pt(tt, &Phi_vec);
    Eigen::Map<NSxNS_type> Phi(Phi_vec.data(), NS, NS);

    vector<double> BB_vec(NS*NI);
    BBs_interp_->pt(tt, &BB_vec);
    Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);
    
    uu = -(*RRinv_) * BB.transpose() * Phi.transpose() * (*eta_star_);
  }
  if (type_ == kLTV_TYPE) {
    vector<double> AA_vec(NS*NS);
    AAs_interp_->pt(tt, &AA_vec);
    Eigen::Map<NSxNS_type> AA(AA_vec.data(), NS, NS);
    
    vector<double> BB_vec(NS*NI);
    BBs_interp_->pt(tt, &BB_vec);
    Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);
    
    if (controltype_ == kMINENERGY_CLOSEDLOOP_CONTROLTYPE) {
      vector<double> KK_vec(NI*NS);
      KKs_interp_->pt(tt, &KK_vec);
      Eigen::Map<NIxNS_type> KK(KK_vec.data(), NI, NS);
      dxx = AA*xx - BB*KK*xx + BB*uu;
    }
    else
      dxx = AA*xx + BB*uu;
  }
  else
    dxx = (*AA_LTI_)*xx + (*BB_LTI_)*uu;
}
bool IntegrateLTIFreeDynamics(const NSxNS_type & AA, const NSxNI_type & BB, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  LinearEq linear_sys(AA, BB);
  return IntegrateForward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegrateLTIFeedForwardDynamics(const NSxNS_type & AA, const NSxNI_type & BB, const InterpVector & uu_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  LinearEq linear_sys(AA, BB, uu_interp);
  return IntegrateForward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegrateLTVFreeDynamics(const InterpVector & AA_interp, const InterpVector & BB_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  LinearEq linear_sys(AA_interp, BB_interp);
  return IntegrateForward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegrateLTVFeedForwardDynamics(const InterpVector & AA_interp, const InterpVector & BB_interp, const InterpVector & uu_interp, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  LinearEq linear_sys(AA_interp, BB_interp, uu_interp);
  return IntegrateForward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegrateLTVMinControlOpenLoopDynamics(const InterpVector & AA_interp, const InterpVector & BB_interp, const NIxNI_type & RRinv, const InterpVector & Phis_interp, const NSx1_type & eta_star, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  LinearEq linear_sys(AA_interp, BB_interp, RRinv, Phis_interp, eta_star);
  return IntegrateForward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}
bool IntegrateLTVMinControlClosedLoopDynamics(const InterpVector & AA_interp, const InterpVector & BB_interp, const InterpVector & KK_interp, const NIxNI_type & RRinv, const InterpVector & Phis_interp, const NSx1_type & eta_star, const ode_state_type & x0, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * xx_vec) {
  LinearEq linear_sys(AA_interp, BB_interp, KK_interp, RRinv, Phis_interp, eta_star);
  return IntegrateForward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}

STMslot2Eq::STMslot2Eq(const NSxNS_type & AA) : AA_LTI_(&AA), type_(kLTI_OPENLOOP_TYPE) {}
STMslot2Eq::STMslot2Eq(const NSxNS_type & AA, const NSxNI_type & BB, const NIxNS_type & KK)
    : AA_LTI_(&AA), BB_LTI_(&BB), KK_LTI_(&KK), type_(kLTI_CLOSEDLOOP_TYPE) {}
STMslot2Eq::STMslot2Eq(const InterpVector & AAs_interp, double t_h)
    : AAs_interp_(&AAs_interp), t_h_(t_h), type_(kLTV_OPENLOOP_TYPE) {}
STMslot2Eq::STMslot2Eq(const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, double t_h)
    : AAs_interp_(&AAs_interp), BBs_interp_(&BBs_interp), KKs_interp_(&KKs_interp), t_h_(t_h), type_(kLTV_CLOSEDLOOP_TYPE) {}
void STMslot2Eq::operator()(const ode_state_type &Phi_vec , ode_state_type &dPhi_vec , const double ss) {
  Eigen::Map<NSxNS_type> dPhi(dPhi_vec.data(), NS, NS);
  Eigen::Map<const NSxNS_type> Phi(Phi_vec.data(), NS, NS);
  if (type_ == kLTV_OPENLOOP_TYPE || type_ == kLTV_CLOSEDLOOP_TYPE) {
    double tt = t_h_-ss;
    vector<double> AA_vec(NS*NS);
    AAs_interp_->pt(tt, &AA_vec);
    Eigen::Map<NSxNS_type> AA(AA_vec.data(), NS, NS);

    //Integrator only integrates forward in time. Otherwise dPhi = -dPhi.
    if (type_ == kLTV_OPENLOOP_TYPE) {
      dPhi = Phi * AA;
    } else { //type == LTV_CLOSEDLOOP_TYPE
      vector<double> KK_vec(NI*NS);
      KKs_interp_->pt(tt, &KK_vec);
      Eigen::Map<NIxNS_type> KK(KK_vec.data(), NI, NS);

      vector<double> BB_vec(NS*NI);
      BBs_interp_->pt(tt, &BB_vec);
      Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);
      
      NSxNS_type AAK = AA-BB*KK;
      dPhi = Phi * AAK;
    }
  }
  else if(type_ == kLTI_OPENLOOP_TYPE)
    dPhi = Phi * (*AA_LTI_);
  else  // type == LTI_CLOSEDLOOP_TYPE
    dPhi = Phi * ((*AA_LTI_)-(*BB_LTI_)*(*KK_LTI_));
}
bool IntegrateLTIOpenLoopSTMSlot2(const NSxNS_type & AA, const NSxNI_type & BB, const NIxNS_type & KK, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec) {
  STMslot2Eq STMslot2_sys(AA, BB, KK);
  ode_state_type identity_vec(NS*NS);
  Eigen::Map<NSxNS_type> identity_mat(identity_vec.data(), NS, NS);
  identity_mat = NSxNS_type::Identity();
  return IntegrateBackward(STMslot2_sys, identity_vec, tt_h, tt_vec, Phi_vec);
}
bool IntegrateLTIClosedLoopSTMSlot2(const NSxNS_type & AA, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec) {
  STMslot2Eq STMslot2_sys(AA);
  ode_state_type identity_vec(NS*NS);
  Eigen::Map<NSxNS_type> identity_mat(identity_vec.data(), NS, NS);
  identity_mat = NSxNS_type::Identity();
  return IntegrateBackward(STMslot2_sys, identity_vec, tt_h, tt_vec, Phi_vec);
}
bool IntegrateLTVOpenLoopSTMSlot2(const InterpVector & AAs_interp, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec) {
  STMslot2Eq STMslot2_sys(AAs_interp, tt_h);
  ode_state_type identity_vec(NS*NS);
  Eigen::Map<NSxNS_type> identity_mat(identity_vec.data(), NS, NS);
  identity_mat = NSxNS_type::Identity();
  return IntegrateBackward(STMslot2_sys, identity_vec, tt_h, tt_vec, Phi_vec);
}
bool IntegrateLTVClosedLoopSTMSlot2 (const InterpVector & AAs_interp, const InterpVector & BBs_interp, const InterpVector & KKs_interp, double tt_h, vector<double> * tt_vec, vector<ode_state_type> * Phi_vec) {
  STMslot2Eq STMslot2_sys(AAs_interp,BBs_interp,KKs_interp, tt_h);
  ode_state_type identity_vec(NS*NS);
  Eigen::Map<NSxNS_type> identity_mat(identity_vec.data(), NS, NS);
  identity_mat = NSxNS_type::Identity();
  return IntegrateBackward(STMslot2_sys, identity_vec, tt_h, tt_vec, Phi_vec);
}
} // namespace sys