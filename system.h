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
    bool integrate_forward (const TT & sys, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec) {
        ode_state_type xx = x0;
        
        double tt = 0;
        double dt = INTEGRATION_DT_START;
        size_t step_num = 0;
        stepper_type stepper;
        
        tt_vec.clear();
        xx_vec.clear();
        tt_vec.push_back(0.0);
        xx_vec.push_back(x0);
        
        odeint::controlled_step_result step_result;
        
        while (tt < tt_h && step_num < MAX_NUM_INTEGRATION_STEPS){
            if (tt_h - tt < dt)
                dt = tt_h - tt;
            step_result = stepper.try_step( sys , xx , tt , dt );
            if(step_result == 0){
                tt_vec.push_back(tt);
                xx_vec.push_back(xx);
            }
            step_num++;
        }
        if (tt >= tt_h)
            return true;
        return false;
    }
    
    //Implemented here because templates are stupid
    template <typename TT>
    bool integrate_backward(const TT & sys, const ode_state_type & x1, double tt_h, vector<double> & tts, vector<ode_state_type> & xxs) {
    
        ode_state_type xx = x1;
        stepper_type stepper;
        
        double ss = 0;
        double ds = INTEGRATION_DT_START;
        size_t step_num = 0;
        
        tts.clear();
        xxs.clear();
        tts.push_back(tt_h);
        xxs.push_back(xx);
        
        odeint::controlled_step_result step_result;
        
        //integrates forward in time so setting ss = tt_h - tt
        while (ss < tt_h && step_num < MAX_NUM_INTEGRATION_STEPS){
            if (tt_h - ss < ds)
                ds = tt_h - ss;
            step_result = stepper.try_step( sys , xx , ss , ds );
            if(step_result == 0){
                tts.push_back(tt_h - ss);
                xxs.push_back(xx);
            }
            step_num++;
        }
        if (ss >= tt_h){
            reverse(tts.begin(), tts.end()); //Reverse to fix integrating forward in time.
            reverse(xxs.begin(), xxs.end());
            return true;
        }
        return false;
    }

    //Implemented here because templates are stupid
    template <typename TT>
    bool integrate_forward_with_constraints (ode_state_type & uu_loc, const TT & sys, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec, vector<ode_state_type> & uu_vec) {
        
        ode_state_type xx = x0;
        double tt = 0.0;
        double dt = 0.00001;
        size_t step_num = 0;
        stepper_type stepper;
        
        tt_vec.clear();
        xx_vec.clear();
        uu_vec.clear();
        tt_vec.push_back(0.0);
        xx_vec.push_back(x0);
        stepper.try_step( sys , xx , tt , dt ); // calling to get uu at time 0 which populating uu_loc
        uu_vec.push_back(uu_loc);
        
        xx = x0;
        tt = 0.0;
        dt = INTEGRATION_DT_START;
        odeint::controlled_step_result step_result;
        
        while ( tt < tt_h && step_num < MAX_NUM_INTEGRATION_STEPS ){
            if ( tt_h - tt < dt )
                dt = tt_h - tt;
            step_result = stepper.try_step( sys , xx , tt , dt ); // note stepper also populates uu_loc (shared memory)
//            cout << "{" << tt << " , " << uu_loc.back() << "}, " << endl;

            if ( !sys.is_constraintsatisfied(xx_vec.back(), xx, uu_loc) )
                return false;
            if ( step_result == 0 ){
                tt_vec.push_back(tt);
                xx_vec.push_back(xx);
                uu_vec.push_back(uu_loc);
            }
            
            step_num++;
        }
//        ode_state_type zero_vec(NI,0.0);
//        uu_vec.push_back(zero_vec); // the control input at the final time (to make all vectors have equal length)
        if ( tt >= tt_h ){
            //check final state for constraints
            dt = .00001;
            step_result = stepper.try_step( sys , xx , tt , dt );
            if ( !sys.is_constraintsatisfied(xx_vec.back(), xx, uu_loc) )
                return false;
            if( step_result == 0 ) {
//                uu_vec.push_back(uu_loc); // the control input at the final time (to make all vectors have equal length)
                return true;
            }
        }
        return false;
    }
    
    bool integrate_free_dynamics ( const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec );

    bool integrate_feed_forward_dynamics ( InterpVector * uus_vec_in, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec );

    bool integrate_point_tracking_dynamics ( InterpVector * KKs_vec_in, const ode_state_type & xx_ref_in, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec );
    
    bool integrate_trajectory_tracking_dynamics ( InterpVector * uu_ff_interp_in, InterpVector * KKs_vec_in, InterpVector * xx_ref_interp_in, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec);
    
    bool integrate_trajectory_tracking_dynamics_w_constraints ( InterpVector * uu_ff_interp_in, InterpVector * KKs_vec_in, InterpVector * xx_ref_interp_in, const ode_state_type & x0, const kin_constraints & constraints_in, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec, vector<ode_state_type> & uu_vec);
    
    bool integrate_trajectory_linear_steering_projection_w_constraints (InterpVector *xxzero, InterpVector *BB, InterpVector *KKlin, InterpVector *KKproj, InterpVector *WWK, InterpVector *Phi, const NSx1_type &eta, const NIxNI_type &RRinv, const ode_state_type & x0, const kin_constraints & constraints_in, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec, vector<ode_state_type> & uu_vec);

    bool integrate_trajectory_linear_steering_projection_w_constraints_w_cost (NSx1_type xx_samp, InterpVector *xxzero, InterpVector *BB, InterpVector *KKlin, InterpVector *KKproj, InterpVector *WWK, InterpVector *Phi, const NSx1_type &eta, const NSxNS_type &QQ, const NIxNI_type &RR, const NSxNS_type &P1, const NIxNI_type &RRinv, const ode_state_type & x0, const kin_constraints & constraints_in, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec, vector<ode_state_type> & uu_vec, double & JJ);

    class RicattiEq {
        NSxNS_type AA_LTI;
        NSxNI_type BB_LTI;
        NSxNS_type QQ;
        NIxNI_type RRinv;
        int type;
        InterpVector * AAs_interp;
        InterpVector * BBs_interp;
        double t_h; // for backward integration

    public:
        RicattiEq (const NSxNS_type & AA_in, const NSxNI_type & BB_in, const NSxNS_type & QQ_in, const NIxNI_type & RRinv_in) : AA_LTI(AA_in), BB_LTI(BB_in), QQ(QQ_in), RRinv(RRinv_in), type(LTI_TYPE) {}
        RicattiEq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in, const NSxNS_type & QQ_in, const NIxNI_type & RRinv_in, double t_h_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), QQ(QQ_in), RRinv(RRinv_in), t_h(t_h_in), type(LTV_TYPE) {}
        void operator()( const ode_state_type &PP_vec , ode_state_type &dPP_vec , const double ss );
    };

    bool integrate_LTI_Ricatti (const NSxNS_type & AA, const NSxNI_type & BB, const NSxNS_type & QQ, const NIxNI_type & RRinv, const ode_state_type & P1, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & PP_vec);

    bool integrate_LTV_Ricatti (InterpVector * AAs_interp, InterpVector * BBs_interp, const NSxNS_type & QQ, const NIxNI_type & RRinv, const ode_state_type & P1, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & PP_vec);

    class ReachabilityEq {
        InterpVector * AAs_interp;
        InterpVector * BBs_interp;
//        InterpVector * RRs_interp;
        InterpVector * KKs_interp;
        InterpVector * WWKs_interp; // FOR grammiantype = SSK_GRAMMIANTYPE
        NIxNI_type RR;
        NIxNI_type RRinv;
        int gramiantype;

    public:
        ReachabilityEq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in, const NIxNI_type & RRinv_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), RRinv(RRinv_in), gramiantype(WW_GRAMIANTYPE) {}
        ReachabilityEq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in, InterpVector * KKs_interp_in, const NIxNI_type & RRinv_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), RRinv(RRinv_in), KKs_interp(KKs_interp_in), gramiantype(WWK_GRAMIANTYPE) {}
        ReachabilityEq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in, InterpVector * KKs_interp_in, InterpVector * WWKs_interp_in, const NIxNI_type & RR_in, const NIxNI_type & RRinv_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), RR(RR_in), RRinv(RRinv_in), KKs_interp(KKs_interp_in), WWKs_interp(WWKs_interp_in), gramiantype(SSK_GRAMIANTYPE) {}

        void operator()( const ode_state_type &WW_vec , ode_state_type &dWW_vec , const double tt );
    };
    
    bool integrate_WW (InterpVector * AAs_interp, InterpVector * BBs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> & tts, vector<ode_state_type> & WWs);
    bool integrate_WWK (InterpVector * AAs_interp, InterpVector * BBs_interp, InterpVector * KKs_interp, const NIxNI_type & RRinv, double tt_h, vector<double> & tts, vector<ode_state_type> & WWs);
    bool integrate_SSK (InterpVector * AAs_interp, InterpVector * BBs_interp, InterpVector * KKs_interp, InterpVector * WWKs_interp, const NIxNI_type & RR, const NIxNI_type & RRinv, double tt_h, vector<double> & tts, vector<ode_state_type> & WWs);
    
    // dxx = AA*xx + BB*uu
    class LinearEq {
        NSxNS_type AA_LTI;
        NSxNI_type BB_LTI;
        NIxNI_type RRinv;
        vector< double > uu_vec;
        int type, controltype;
        InterpVector * AAs_interp;
        InterpVector * BBs_interp;
        InterpVector * uus_ff_interp;
        InterpVector * Phis_interp;
        InterpVector * KKs_interp;
        NSx1_type eta_star;

    public:
        LinearEq (const NSxNS_type & AA_in, const NSxNI_type & BB_in) : AA_LTI(AA_in), BB_LTI(BB_in), type(LTI_TYPE), controltype(FREE_CONTROLTYPE), uu_vec(NI,0.0) {}
        LinearEq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), uu_vec(NI,0.0), type(LTV_TYPE), controltype(FREE_CONTROLTYPE) {}
        LinearEq (const NSxNS_type & AA_in, const NSxNI_type & BB_in, InterpVector * uus_ff_interp_in) : AA_LTI(AA_in), BB_LTI(BB_in), uus_ff_interp(uus_ff_interp_in), uu_vec(NI), type(LTI_TYPE), controltype(FEEDFORWARD_CONTROLTYPE) {}
        LinearEq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in, InterpVector * uus_ff_interp_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), uus_ff_interp(uus_ff_interp_in), uu_vec(NI), type(LTV_TYPE), controltype(FEEDFORWARD_CONTROLTYPE) {}
        LinearEq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in, const NIxNI_type & RRinv_in, InterpVector * Phis_interp_in, const NSx1_type & eta_star_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), RRinv(RRinv_in), Phis_interp(Phis_interp_in), eta_star(eta_star_in), uu_vec(NI), type(LTV_TYPE), controltype(MINENERGY_OPENLOOP_CONTROLTYPE) {}
        LinearEq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in, InterpVector * KKs_interp_in, const NIxNI_type & RRinv_in, InterpVector * Phis_interp_in, const NSx1_type & eta_star_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), KKs_interp(KKs_interp_in), RRinv(RRinv_in), Phis_interp(Phis_interp_in), eta_star(eta_star_in), uu_vec(NI), type(LTV_TYPE), controltype(MINENERGY_CLOSEDLOOP_CONTROLTYPE) {}
        void operator()( const ode_state_type &xx_vec , ode_state_type &dxx_vec , const double tt );
    };
    
    bool integrate_LTI_free_dynamics (const NSxNS_type & AA, const NSxNI_type & BB, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec);

    bool integrate_LTI_feed_forward_dynamics (const NSxNS_type & AA, const NSxNI_type & BB, InterpVector * uu_interp, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec);

    bool integrate_LTV_free_dynamics (InterpVector * AA_interp, InterpVector * BB_interp, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec);

    bool integrate_LTV_feed_forward_dynamics (InterpVector * AA_interp, InterpVector * BB_interp, InterpVector * uu_interp, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec);

    bool integrate_LTV_min_control_open_loop_dynamics (InterpVector * AA_interp, InterpVector * BB_interp, const NIxNI_type & RRinv, InterpVector * Phis_interp, const NSx1_type & eta_star, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec);
    
    bool integrate_LTV_min_control_closed_loop_dynamics (InterpVector * AA_interp, InterpVector * BB_interp, InterpVector * KK_interp, const NIxNI_type & RRinv, InterpVector * Phis_interp, const NSx1_type & eta_star, const ode_state_type & x0, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & xx_vec);

    class STMslot2Eq { // dPhi(tt_h, t) = -Phi(tt_h, t) A(t)
        NSxNS_type AA_LTI;
        NSxNI_type BB_LTI;
        NIxNS_type KK_LTI;
        int type;
        InterpVector * AAs_interp;
        InterpVector * BBs_interp;
        InterpVector * KKs_interp;
        double t_h;  // For backward integration

    public:
        STMslot2Eq (const NSxNS_type & AA_in) : AA_LTI(AA_in), type(LTI_OPENLOOP_TYPE) {}
        STMslot2Eq (const NSxNS_type & AA_in, const NSxNI_type & BB_in, const NIxNS_type & KK_in) : AA_LTI(AA_in), BB_LTI(BB_in), KK_LTI(KK_in), type(LTI_CLOSEDLOOP_TYPE) {}
        STMslot2Eq (InterpVector * AAs_interp_in, double t_h_in) : AAs_interp(AAs_interp_in), t_h(t_h_in), type(LTV_OPENLOOP_TYPE) {}
        STMslot2Eq (InterpVector * AAs_interp_in, InterpVector * BBs_interp_in, InterpVector * KKs_interp_in, double t_h_in) : AAs_interp(AAs_interp_in), BBs_interp(BBs_interp_in), KKs_interp(KKs_interp_in), t_h(t_h_in), type(LTV_CLOSEDLOOP_TYPE) {}
        void operator()( const ode_state_type &Phi_vec , ode_state_type &dPhi_vec , const double ss );
    };

    bool integrate_LTI_openloop_STMslot2 (const NSxNS_type & AA, const NSxNI_type & BB, const NIxNS_type & KK, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & Phi_vec);

    bool integrate_LTI_closedloop_STMslot2 (const NSxNS_type & AA, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & Phi_vec);

    bool integrate_LTV_openloop_STMslot2 (InterpVector * AAs_interp, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & Phi_vec);

    bool integrate_LTV_closedloop_STMslot2 (InterpVector * AAs_interp, InterpVector * BBs_interp, InterpVector * KKs_interp, double tt_h, vector<double> & tt_vec, vector<ode_state_type> & Phi_vec);
}






#endif /* defined(__Agile_RRT__system__) */