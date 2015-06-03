//
//  system.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 5/10/15.
//  Copyright (c) 2015 Timothy Caldwell. All rights reserved.
//

#include "system.h"

namespace sys
{
bool integrate_free_dynamics ( const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec ) {
    PendCart pendcart_sys;
    return integrate_forward(pendcart_sys, x0, tt_h, tt_vec, xx_vec);
}

bool integrate_feed_forward_dynamics ( InterpVector * uus_vec_in, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec ) {
    PendCart pendcart_sys(uus_vec_in);
    return integrate_forward(pendcart_sys, x0, tt_h, tt_vec, xx_vec);
}

bool integrate_point_tracking_dynamics ( InterpVector * KKs_vec_in, const ode_state_type & xx_ref_in, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec ) {
    PendCart pendcart_sys(KKs_vec_in, xx_ref_in);
    return integrate_forward(pendcart_sys, x0, tt_h, tt_vec, xx_vec);
}

bool integrate_trajectory_tracking_dynamics ( InterpVector * uu_ff_interp_in, InterpVector * KKs_vec_in, InterpVector * xx_ref_interp_in, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec) {
    PendCart pendcart_sys(uu_ff_interp_in, KKs_vec_in, xx_ref_interp_in);
    return integrate_forward(pendcart_sys, x0, tt_h, tt_vec, xx_vec);
}

bool integrate_trajectory_tracking_dynamics_w_constraints ( InterpVector * uu_ff_interp_in, InterpVector * KKs_vec_in, InterpVector * xx_ref_interp_in, const ode_state_type & x0, const kin_constraints & constraints_in, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec, vector<ode_state_type> & uu_vec) {

    ode_state_type uu_loc(NI);
    PendCart pendcart_sys(&uu_loc, uu_ff_interp_in, KKs_vec_in, xx_ref_interp_in, constraints_in);
    return integrate_forward_with_constraints(uu_loc, pendcart_sys, x0, tt_h, tt_vec, xx_vec, uu_vec);
}

bool integrate_trajectory_linear_steering_projection_w_constraints (InterpVector *xxzero, InterpVector *BB, InterpVector *KKlin, InterpVector *KKproj, InterpVector *WWK, InterpVector *Phi, const NSx1_type &eta, const NIxNI_type &RRinv, const ode_state_type & x0, const kin_constraints & constraints_in, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec, vector<ode_state_type> & uu_vec) {

    ode_state_type uu_loc(NI);
    PendCart pendcart_sys(&uu_loc, xxzero, BB, KKlin, KKproj, WWK, Phi, eta, RRinv, constraints_in);
    return integrate_forward_with_constraints(uu_loc, pendcart_sys, x0, tt_h, tt_vec, xx_vec, uu_vec);
}

bool integrate_trajectory_linear_steering_projection_w_constraints_w_cost (NSx1_type xx_samp, InterpVector *xxzero, InterpVector *BB, InterpVector *KKlin, InterpVector *KKproj, InterpVector *WWK, InterpVector *Phi, const NSx1_type &eta, const NSxNS_type &QQ, const NIxNI_type &RR, const NSxNS_type &P1, const NIxNI_type &RRinv, const ode_state_type & x0, const kin_constraints & constraints_in, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec, vector<ode_state_type> & uu_vec, double & JJ) {

    ode_state_type x0J0 = x0;
    x0J0.push_back(0.0);
    
    ode_state_type uu_loc(NI);
    PendCart pendcart_sys(&uu_loc, xxzero, BB, KKlin, KKproj, WWK, Phi, eta, QQ, RR, RRinv, constraints_in);
    bool integration_success = integrate_forward_with_constraints(uu_loc, pendcart_sys, x0J0, tt_h, tt_vec, xx_vec, uu_vec);

    Eigen::Map<const NSx1_type> xx_mat(xx_vec.back().data(), NS, 1);
    JJ = xx_vec.back().back() /* ell */ + 0.5*(xx_mat-xx_samp).dot(P1*(xx_mat-xx_samp)) /* m */;
    
    //remove int_0^t ell dtau from xx_vec
    for(vector<ode_state_type>::iterator it = xx_vec.begin(); it!=xx_vec.end(); ++it)
        it->pop_back();
    
    return integration_success;

}

void RicattiEq::operator()( const ode_state_type &PP_vec , ode_state_type &dPP_vec , const double ss ) {
    Eigen::Map<const NSxNS_type> PP(PP_vec.data(), NS, NS);
    Eigen::Map<NSxNS_type> dPP(dPP_vec.data(), NS, NS); // mat_dPP and dPP share memory

    if (type == LTV_TYPE){
        double tt = t_h-ss;
        vector<double> AA_vec(NS*NS);
        AAs_interp->pt(tt, AA_vec);
        Eigen::Map<NSxNS_type> AA(AA_vec.data(), NS, NS);
        
        vector<double> BB_vec(NS*NI);
        BBs_interp->pt(tt, BB_vec);
        Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);
        
        //Integrator only integrates forward in time. Otherwise dPP = -dPP.
        dPP = AA.transpose()*PP + PP*AA - PP*BB*RRinv*BB.transpose()*PP + QQ;
    }
    else
        //Integrator only integrates forward in time. Otherwise dPP = -dPP.
        dPP = AA_LTI.transpose()*PP + PP*AA_LTI - PP*BB_LTI*RRinv*BB_LTI.transpose()*PP + QQ;

}

bool integrate_LTI_Ricatti (const NSxNS_type & AA, const NSxNI_type & BB, const NSxNS_type & QQ, const NIxNI_type & RRinv, const ode_state_type & P1, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & PP_vec) {
    RicattiEq ric_sys(AA, BB, QQ, RRinv);
    return integrate_backward(ric_sys, P1, tt_h, tt_vec, PP_vec);
}

bool integrate_LTV_Ricatti (InterpVector * AAs_interp, InterpVector * BBs_interp, const NSxNS_type & QQ, const NIxNI_type & RRinv, const ode_state_type & P1, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & PP_vec) {
    RicattiEq ric_sys(AAs_interp, BBs_interp, QQ, RRinv, tt_h);
    return integrate_backward(ric_sys, P1, tt_h, tt_vec, PP_vec);
}

void ReachabilityEq::operator()( const ode_state_type &WW_vec , ode_state_type &dWW_vec , const double tt ) {
    Eigen::Map<const NSxNS_type> WW(WW_vec.data(), NS, NS);
    Eigen::Map<NSxNS_type> dWW(dWW_vec.data(), NS, NS);

    vector<double> AA_vec(NS*NS);
    AAs_interp->pt(tt, AA_vec);
    Eigen::Map<NSxNS_type> AA(AA_vec.data(), NS, NS);
    
    vector<double> BB_vec(NS*NI);
    BBs_interp->pt(tt, BB_vec);
    Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);

    if (gramiantype == WW_GRAMIANTYPE)
        dWW = AA*WW + WW*AA.transpose() + BB*RRinv*BB.transpose();
    else {
        vector<double> KK_vec(NI*NS);
        KKs_interp->pt(tt, KK_vec);
        Eigen::Map<NIxNS_type> KK(KK_vec.data(), NI, NS);
        
        NSxNS_type AAK = AA-BB*KK;
        
        if (gramiantype == WWK_GRAMIANTYPE) {
            dWW = AAK*WW + WW*(AAK.transpose()) + BB*RRinv*BB.transpose();
        }
        else { //gramiantype = SSK_GRAMMIANTYPE
            vector<double> WWK_vec(NS*NS);
            WWKs_interp->pt(tt, WWK_vec);
            Eigen::Map<NSxNS_type> WWK(WWK_vec.data(), NS, NS);

            NIxNS_type CC = KK*WWK - RRinv*BB.transpose();
            
            dWW = AAK*WW + WW*(AAK.transpose()) + CC.transpose()*RR*CC;
        }
    }
}

bool integrate_WW (InterpVector * AAs_interp, InterpVector * BBs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> & tts, vector<ode_state_type> & WWs) {
    ReachabilityEq reach_sys(AAs_interp, BBs_interp, RRinv);
    ode_state_type WW0(NS*NS,0.0); //WW(0.0) = 0 matrix
    return integrate_forward(reach_sys, WW0, tt_h, tts, WWs);
}
bool integrate_WWK (InterpVector * AAs_interp, InterpVector * BBs_interp, InterpVector * KKs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> & tts, vector<ode_state_type> & WWs) {
    ReachabilityEq reach_sys(AAs_interp, BBs_interp, KKs_interp, RRinv);
    ode_state_type WW0(NS*NS,0.0); //WW(0.0) = 0 matrix
    return integrate_forward(reach_sys, WW0, tt_h, tts, WWs);
}
bool integrate_SSK (InterpVector * AAs_interp, InterpVector * BBs_interp, InterpVector * KKs_interp, InterpVector * WWKs_interp, const NIxNI_type & RR, const NIxNI_type & RRinv, double tt_h, vector<double> & tts, vector<ode_state_type> & WWs) {
    ReachabilityEq reach_sys(AAs_interp, BBs_interp, KKs_interp, WWKs_interp, RR, RRinv);
    ode_state_type WW0(NS*NS,0.0); //WW(0.0) = 0 matrix
    return integrate_forward(reach_sys, WW0, tt_h, tts, WWs);
}

// dxx = AA*xx + BB*uu
void LinearEq::operator()( const ode_state_type &xx_vec , ode_state_type &dxx_vec , const double tt ) {
    Eigen::Map<const NSx1_type> xx(xx_vec.data(), NS, 1);
    Eigen::Map<NSx1_type> dxx(dxx_vec.data(), NS, 1);
    Eigen::Map<NIx1_type> uu(uu_vec.data(), NI, 1);

    if (controltype == FEEDFORWARD_CONTROLTYPE)
        // *** feedforward:  u = u_ff *** //
        uus_ff_interp->pt(tt, uu_vec);
    else if (controltype == MINENERGY_OPENLOOP_CONTROLTYPE || controltype == MINENERGY_CLOSEDLOOP_CONTROLTYPE){
        // *** MIN ENERGY OPEN LOOP:  u = -RRinv * BB^T * Phi^T * eta *** //
        vector<double> Phi_vec(NS*NS);
        Phis_interp->pt(tt, Phi_vec);
        Eigen::Map<NSxNS_type> Phi(Phi_vec.data(), NS, NS);

        vector<double> BB_vec(NS*NI);
        BBs_interp->pt(tt, BB_vec);
        Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);
        
//                Eigen::Map<const NSxNI_type> eta_star_mat(eta_star.data(), NS, NI);
        uu = - RRinv * BB.transpose() * Phi.transpose() * eta_star;
    }
    if (type == LTV_TYPE){
        vector<double> AA_vec(NS*NS);
        AAs_interp->pt(tt, AA_vec);
        Eigen::Map<NSxNS_type> AA(AA_vec.data(), NS, NS);
        
        vector<double> BB_vec(NS*NI);
        BBs_interp->pt(tt, BB_vec);
        Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);
        
        if (controltype == MINENERGY_CLOSEDLOOP_CONTROLTYPE) {
            vector<double> KK_vec(NI*NS);
            KKs_interp->pt(tt, KK_vec);
            Eigen::Map<NIxNS_type> KK(KK_vec.data(), NI, NS);
            dxx = AA*xx - BB*KK*xx + BB*uu;
        }
        else
            dxx = AA*xx + BB*uu;
    }
    else
        dxx = AA_LTI*xx + BB_LTI*uu;

}

bool integrate_LTI_free_dynamics (const NSxNS_type & AA, const NSxNI_type & BB, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec) {
    LinearEq linear_sys(AA, BB);
    return integrate_forward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}

bool integrate_LTI_feed_forward_dynamics (const NSxNS_type & AA, const NSxNI_type & BB, InterpVector * uu_interp, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec) {
    LinearEq linear_sys(AA, BB, uu_interp);
    return integrate_forward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}

bool integrate_LTV_free_dynamics (InterpVector * AA_interp, InterpVector * BB_interp, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec) {
    LinearEq linear_sys(AA_interp, BB_interp);
    return integrate_forward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}

bool integrate_LTV_feed_forward_dynamics (InterpVector * AA_interp, InterpVector * BB_interp, InterpVector * uu_interp, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec) {
    LinearEq linear_sys(AA_interp, BB_interp, uu_interp);
    return integrate_forward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}

bool integrate_LTV_min_control_open_loop_dynamics (InterpVector * AA_interp, InterpVector * BB_interp, const NIxNI_type & RRinv, InterpVector * Phis_interp, const NSx1_type & eta_star, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec) {
    LinearEq linear_sys(AA_interp, BB_interp, RRinv, Phis_interp, eta_star);
    return integrate_forward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}
bool integrate_LTV_min_control_closed_loop_dynamics (InterpVector * AA_interp, InterpVector * BB_interp, InterpVector * KK_interp, const NIxNI_type & RRinv, InterpVector * Phis_interp, const NSx1_type & eta_star, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec) {
    LinearEq linear_sys(AA_interp, BB_interp, KK_interp, RRinv, Phis_interp, eta_star);
    return integrate_forward(linear_sys, x0, tt_h, tt_vec, xx_vec);
}

void STMslot2Eq::operator()( const ode_state_type &Phi_vec , ode_state_type &dPhi_vec , const double ss ) {
    Eigen::Map<const NSxNS_type> Phi(Phi_vec.data(), NS, NS);
    Eigen::Map<NSxNS_type> dPhi(dPhi_vec.data(), NS, NS);

    if (type == LTV_OPENLOOP_TYPE || type == LTV_CLOSEDLOOP_TYPE){
        double tt = t_h-ss;

        vector<double> AA_vec(NS*NS);
        AAs_interp->pt(tt, AA_vec);
        Eigen::Map<NSxNS_type> AA(AA_vec.data(), NS, NS);

        //Integrator only integrates forward in time. Otherwise dPhi = -dPhi.
        if (type == LTV_OPENLOOP_TYPE)
            dPhi = Phi * AA;
        else { //type == LTV_CLOSEDLOOP_TYPE
            vector<double> KK_vec(NI*NS);
            KKs_interp->pt(tt, KK_vec);
            Eigen::Map<NIxNS_type> KK(KK_vec.data(), NI, NS);

            vector<double> BB_vec(NS*NI);
            BBs_interp->pt(tt, BB_vec);
            Eigen::Map<NSxNI_type> BB(BB_vec.data(), NS, NI);
            
            NSxNS_type AAK = AA-BB*KK;
            dPhi = Phi * AAK;
        }
    }
    else if(type == LTI_OPENLOOP_TYPE)
        dPhi = Phi * AA_LTI;
    else  // type == LTI_CLOSEDLOOP_TYPE
        dPhi = Phi * (AA_LTI-BB_LTI*KK_LTI);
}

bool integrate_LTI_openloop_STMslot2 (const NSxNS_type & AA, const NSxNI_type & BB, const NIxNS_type & KK, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & Phi_vec) {
    STMslot2Eq STMslot2_sys(AA, BB, KK);
    ode_state_type identity_vec(NS*NS);
    Eigen::Map< NSxNS_type > identity_mat(identity_vec.data(), NS, NS);
    identity_mat = NSxNS_type::Identity();
    return integrate_backward(STMslot2_sys, identity_vec, tt_h, tt_vec, Phi_vec);
}
bool integrate_LTI_closedloop_STMslot2 (const NSxNS_type & AA, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & Phi_vec) {
    STMslot2Eq STMslot2_sys(AA);
    ode_state_type identity_vec(NS*NS);
    Eigen::Map< NSxNS_type > identity_mat(identity_vec.data(), NS, NS);
    identity_mat = NSxNS_type::Identity();
    return integrate_backward(STMslot2_sys, identity_vec, tt_h, tt_vec, Phi_vec);
}
bool integrate_LTV_openloop_STMslot2 (InterpVector * AAs_interp, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & Phi_vec) {
    STMslot2Eq STMslot2_sys(AAs_interp, tt_h);
    ode_state_type identity_vec(NS*NS);
    Eigen::Map< NSxNS_type > identity_mat(identity_vec.data(), NS, NS);
    identity_mat = NSxNS_type::Identity();
    return integrate_backward(STMslot2_sys, identity_vec, tt_h, tt_vec, Phi_vec);
}
bool integrate_LTV_closedloop_STMslot2 (InterpVector * AAs_interp, InterpVector * BBs_interp, InterpVector * KKs_interp, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & Phi_vec) {
    STMslot2Eq STMslot2_sys(AAs_interp,BBs_interp,KKs_interp, tt_h);
    ode_state_type identity_vec(NS*NS);
    Eigen::Map< NSxNS_type > identity_mat(identity_vec.data(), NS, NS);
    identity_mat = NSxNS_type::Identity();
    return integrate_backward(STMslot2_sys, identity_vec, tt_h, tt_vec, Phi_vec);
}
}