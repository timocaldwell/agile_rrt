//
//  treevertex.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 5/10/15.
//  Copyright (c) 2015 Timothy Caldwell. All rights reserved.
//

#include "treevertex.h"

TreeVertex::TreeVertex (ode_state_type x0_in, bool usezerotraj_in, double t_h_max_upper_in, InterpVector * xxedge_in, InterpVector * uuedge_in, const kin_constraints & constraints_in, const NSxNS_type & QQ_in, const NIxNI_type & RR_in, const NSxNS_type & P1_in, const NSxNS_type & QQlqr_in, const NIxNI_type & RRlqr_in, const NSxNS_type & P1lqr_in, double t_h_bar_in, bool save_edge_in)
    : usezerotraj(usezerotraj_in), t_h_max_upper(t_h_max_upper_in), QQlqr(QQlqr_in), QQ(QQ_in), RR(RR_in), RRinv(RR.inverse()), RRlqr(RRlqr_in), RRlqr_inv(RRlqr_in.inverse()), P1(P1_in), P1lqr(P1lqr_in), constraints(constraints_in), t_h_bar(t_h_bar_in), computedisttogoal(false), xxgoal(NSx1_type::Zero()), disttogoal_cur(LARGENUM), save_edge(save_edge_in)
{
    t_h_cur = LARGENUM;
    JJapprox_cur = LARGENUM;
    JJproj_cur = LARGENUM;
    x0 = x0_in;
    
    xxproj_cur = nullptr;
    uuproj_cur = nullptr;

    if ( save_edge ) {
        xxedge = new InterpVector(*xxedge_in); // deep copy
        uuedge = new InterpVector(*uuedge_in);
    }
    else {
        xxedge = nullptr;
        uuedge = nullptr;
    }
        
    initialize();
}

TreeVertex::~TreeVertex()
{
    if ( xxproj_cur )
        delete xxproj_cur;
    if ( uuproj_cur)
        delete uuproj_cur;
    if ( xxedge )
        delete xxedge;
    if ( uuedge )
       delete uuedge;
    delete xxzero;
    delete AA;
    delete BB;
    delete KK;
    delete WWK;
    delete SSK;
}

void TreeVertex::initialize()
{
    vector<double> tts;
    vector< ode_state_type > xxs;

    if ( usezerotraj ){
        //------------- FOR Linearizing about zero trajectory!--------------
        if( !sys::integrate_free_dynamics (x0, t_h_max_upper, tts, xxs) )
            cout << "ERROR in treevertex: Zero dynamics integration failed!";
        xxzero = new InterpVector(tts, xxs);
    }
    else {
        //------------- FOR Linearizing about point!--------------
        tts.push_back(0.0);
        tts.push_back(t_h_max_upper/2.0);
        tts.push_back(t_h_max_upper);
        xxs.push_back(x0);
        xxs.push_back(x0);
        xxs.push_back(x0);                
        xxzero = new InterpVector(tts,xxs);
    }
    
    PendCart pendcart;
    vector< double > uuzero (NI, 0.0);
    vector< vector<double> > AAs(tts.size(),vector<double>(NS*NS));
    vector< vector<double> > BBs(tts.size(),vector<double>(NS*NI));
    for(int ii = 0; ii < tts.size(); ++ii) {
        Eigen::Map< NSxNS_type > AA_mat(AAs[ii].data(), NS, NS);
        pendcart.calc_AA(xxs[ii], uuzero, AA_mat);
        
        Eigen::Map< NSxNI_type > BB_mat(BBs[ii].data(), NS, NI);
        pendcart.calc_BB(xxs[ii], uuzero, BB_mat);
    }
    AA = new InterpVector(tts, AAs);
    BB = new InterpVector(tts, BBs);
    
    tts.clear();
    vector<ode_state_type > PPs;
    ode_state_type P1lqr_vec(NS*NS);
    Eigen::Map<NSxNS_type> P1lqr_mat(P1lqr_vec.data(), NS, NS);
    P1lqr_mat = P1lqr;

    if( !sys::integrate_LTV_Ricatti (AA, BB, QQlqr, RRlqr_inv, P1lqr_vec, t_h_max_upper, tts, PPs))
        cout << "ERROR in treevertex: Ricatti integration failed!";
//        InterpVector *PP = new InterpVector(tts, PPs);
    
    vector< vector<double> > KKs(tts.size(),vector<double>(NI*NS));

//        clock_t t_start3 = clock();
    for(int ii = 0; ii < tts.size(); ++ii) {
        Eigen::Map< NSxNS_type > PP_mat(PPs[ii].data(), NS, NS);

        vector<double> BB_vec(NS*NI);
        BB->pt(tts[ii], BB_vec);
        Eigen::Map<NSxNI_type> BB_mat(BB_vec.data(), NS, NI);

        Eigen::Map<NIxNS_type> KK_mat(KKs[ii].data(), NI, NS);
        KK_mat = RRinv * BB_mat.transpose() * PP_mat;
    }
//        cout << "avg time: " << ((float)(clock()-t_start3))/CLOCKS_PER_SEC << " seconds." << endl;
    KK = new InterpVector(tts, KKs);

    tts.clear();
    vector< ode_state_type > WWKs;
    if( !sys::integrate_WWK (AA, BB, KK, RRinv, t_h_max_upper,  tts, WWKs))
        cout << "ERROR in treevertex: WWK integration failed!";
    WWK = new InterpVector(tts, WWKs);

    tts.clear();
    vector< ode_state_type > SSKs;
    if( !sys::integrate_SSK (AA, BB, KK, WWK, RR, RRinv, t_h_max_upper,  tts, SSKs))
        cout << "ERROR in treevertex: SSK integration failed!";
    SSK = new InterpVector(tts, SSKs);

////        tts.clear();
////        vector< ode_state_type > Phis;
////        if( !sys::integrate_LTV_closedloop_STMslot2 (AA, BB, KK, t_h_MAX_UPPER, tts, Phis))
////            cout << "ERROR in treevertex: STM integration failed!";
////        Phi = new InterpVector(tts, Phis);
}


double TreeVertex::JJ_linear_steer (const NSx1_type & xxsamp, const double t_h)
{
    if (t_h > t_h_bar)
        return LARGENUM;
    
    t_h_cur = t_h;
    xxsamp_cur = xxsamp;
    
    vector<double> xxzero_vec(NS*1);
    xxzero->pt(t_h, xxzero_vec);
    Eigen::Map<NSx1_type> xxzero_mat(xxzero_vec.data(), NS, 1);
    
    vector<double> WWK_vec(NS*NS);
    WWK->pt(t_h, WWK_vec);
    Eigen::Map<NSxNS_type> WWK_mat(WWK_vec.data(), NS, NS);
    
    vector<double> SSK_vec(NS*NS);
    SSK->pt(t_h, SSK_vec);
    Eigen::Map<NSxNS_type> SSK_mat(SSK_vec.data(), NS, NS);
    
    Pth_cur = (WWK_mat*P1*WWK_mat + SSK_mat).inverse()*WWK_mat*P1;
    xxzeroerror_cur = xxzero_mat - xxsamp;
    eta_cur = Pth_cur * xxzeroerror_cur;
    
    JJapprox_cur = 0.5*eta_cur.dot(SSK_mat*eta_cur) + 0.5*(xxzeroerror_cur - WWK_mat*eta_cur).dot(P1*(xxzeroerror_cur - WWK_mat*eta_cur));
    return JJapprox_cur;
}

void TreeVertex::set_xxgoal (const NSx1_type & xxgoal_in)
{
    xxgoal = xxgoal_in; computedisttogoal = true;
}

bool TreeVertex::projection (const NSx1_type & xxsamp, const double t_h)
{
    if (xxsamp_cur!=xxsamp || t_h_cur!=t_h) // if linear steering computations have not already been made:
        JJ_linear_steer (xxsamp, t_h);
    
    vector<double> tts;
    vector< ode_state_type > Phis;
//        if( !sys::integrate_LTV_openloop_STMslot2 (AA, t_h, tts, Phis))
    if( !sys::integrate_LTV_closedloop_STMslot2 (AA, BB, KK, t_h, tts, Phis))
        cout << "ERROR in linear_steer: STM integration failed!";

//        if( tts.size() <= 2 ) {
//            cout << "PROJECTION FAILED ERROR in linear_steer: STM integration failed!";
//            return false;
//        }
    
    InterpVector Phi(tts, Phis);// = new InterpVector(tts, Phis);
    

    tts.clear();
    vector< ode_state_type > xxs;
    vector< ode_state_type > uus;
    
//        if (!sys::integrate_trajectory_linear_steering_projection_w_constraints (xxzero, BB, KK, KK, WWK, Phi, eta_cur, RRinv, x0, constraints, t_h, tts, xxs, uus)) {
    if (!sys::integrate_trajectory_linear_steering_projection_w_constraints_w_cost (xxsamp, xxzero, BB, KK, KK, WWK, &Phi, eta_cur, QQ, RR, P1, RRinv, x0, constraints, t_h, tts, xxs, uus, JJproj_cur)) {
//            cout << "WARNING in TreeVertex::projection integration did not finish!";
        return false;
    }

    if(computedisttogoal) {
        disttogoal_cur = LARGENUM;
        for (int ii = 0; ii < xxs.size(); ++ii) {
            Eigen::Map<NSx1_type> xx_mat(xxs[ii].data(), NS, 1);
            double tmp  = 0.5 * (xxgoal-xx_mat).dot(P1*(xxgoal-xx_mat));
            if ( tmp < disttogoal_cur )
                disttogoal_cur = tmp;
        }
    }
    
    //remove old if needed
    if ( xxproj_cur )
        delete xxproj_cur;
    if ( uuproj_cur)
        delete uuproj_cur;

    
    xxproj_cur = new InterpVector(tts, xxs);
    uuproj_cur = new InterpVector(tts, uus);
    
    
    
//        cout << "JJproj ------- " << JJproj_cur << endl;
    
    return true;
}

void TreeVertex::print_traj(double dt, const string & name, ostream & stream)
{
    xxproj_cur->print(dt, name+"xxproj", stream);
    uuproj_cur->print(dt, name+"uuproj", stream);
}

void TreeVertex::print_edge_traj(double dt, const string & name, ostream & stream)
{
    xxedge->print(dt, name+"xxedge", stream);
    uuedge->print(dt, name+"uuedge", stream);
}

TreeVertex TreeVertex::new_vertex()
{
    ode_state_type x0_new;
    xxproj_cur->pt(t_h_cur, x0_new);
    
    ///////////////////////////////
    //
    //SPECIAL CODE FOR THE PENDULUM FOR ANGLE WRAPPING
    //
    ///////////////////////////////
    for ( int ii = 0; ii < NS/2-1; ++ii ) {
        if ( x0_new[2*ii] < -M_PI && x0_new[2*ii+1] < 0 )
            x0_new[2*ii] += 2*M_PI;
        else if ( x0_new[2*ii] > M_PI && x0_new[2*ii+1] > 0 )
            x0_new[2*ii] -= 2*M_PI;
    }
    ///////////////////////////////
    //
    //END SPECIAL CODE FOR THE PENDULUM
    //
    ///////////////////////////////
    
    TreeVertex new_vert(x0_new, usezerotraj, t_h_max_upper, xxproj_cur, uuproj_cur, constraints, QQ, RR, P1, QQlqr, RRlqr, P1lqr, t_h_bar, save_edge);
    if ( computedisttogoal )
        new_vert.set_xxgoal(xxgoal);
    return new_vert;
}
