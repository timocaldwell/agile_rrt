//
//  treevertex.h
//  Agile_RRT
//
//  Created by Timothy Caldwell on 3/16/15.
//  Copyright (c) 2015 Timothy Caldwell. All rights reserved.
//


#ifndef __Agile_RRT__treevertex__
#define __Agile_RRT__treevertex__

#include "global.h"
#include "interp.h"
#include "system.h"
#include "pendcart_3link.h"

class TreeVertex {
private:
    ode_state_type x0;
    InterpVector *xxedge, *uuedge, *xxzero, *AA, *BB, *KK, *WWK, *SSK;
    const NSxNS_type QQlqr, P1lqr, QQ, P1;
    const NIxNI_type RR, RRinv, RRlqr, RRlqr_inv;
    bool usezerotraj, computedisttogoal, save_edge;
    double t_h_max_upper, t_h_bar;

    const kin_constraints constraints;
    
    double t_h_cur, JJapprox_cur, JJproj_cur, disttogoal_cur;
    NSxNS_type Pth_cur;
    NSx1_type xxzeroerror_cur, xxsamp_cur, eta_cur;
    InterpVector *xxproj_cur, *uuproj_cur;
    
    NSx1_type xxgoal;
   
    void initialize();
    
public:
    TreeVertex (ode_state_type x0_in, bool usezerotraj_in, double t_h_max_upper_in, InterpVector * xxedge_in, InterpVector * uuedge_in, const kin_constraints & constraints_in, const NSxNS_type & QQ_in, const NIxNI_type & RR_in, const NSxNS_type & P1_in, const NSxNS_type & QQlqr_in, const NIxNI_type & RRlqr_in, const NSxNS_type & P1lqr_in, double t_h_bar_in, bool save_edge_in);
    
    ~TreeVertex();
    
    double JJ_linear_steer (const NSx1_type & xxsamp, const double t_h);
    
    void get_xxzero_pt(double tt, ode_state_type & xxzero_pt) { xxzero->pt(tt, xxzero_pt); }
    
    void set_xxgoal (const NSx1_type & xxgoal_in);
    
    double get_disttogoal() { return disttogoal_cur; }
    
    bool projection (const NSx1_type & xxsamp, const double t_h);
    
    void print_traj(double dt, const string & name, ostream & stream);
    
    void print_edge_traj(double dt, const string & name, ostream & stream);
    
    TreeVertex new_vertex();    
};

#endif /* defined(__Agile_RRT__treevertex__) */
