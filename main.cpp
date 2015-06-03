//
//  main.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 9/19/14.
//  Copyright (c) 2014 Timothy Caldwell. All rights reserved.
//

#include "global.h"
#include "interp.h"
#include "system.h"
#include "pendcart_3link.h"
#include "treevertex.h"
#include "tree.h"
//#include "rrtreach.h"

int main(int argc, const char * argv[])
{
    // ARGUMENTS : [t_h_max_upper = 1.0] , [ usezerotraj = true] , [max_cnt = 1000], [max_miss = 1000], [max_samplemiss = 10000], [numbruns = 1] , [stopdist = 2.0] , [printskip = 1] , [nndelta = LARGENUM], [stats_name = "0" (use "0" or no arg to not make file)] , [traj_name = "0" (use "0" or no arg to not make file)] , [seed = 0]

    double t_h_max_upper = 1.0;
    bool usezerotraj = true;
    int max_cnt = 1000;//20000;
    int max_miss = 1000;//20000;
    int max_samplemiss = 10000;//30000;
    int numbruns = 1;//13;//100;
    double stopdist = 0.0;
    bool computexxgoaldist = false;
    int printskip = 1;//1000;
    double nndelta = LARGENUM;
    ofstream stats_f, traj_f;
    bool is_stats_out = false;
    bool is_traj_out = false;
    bool save_edges = false;
    double seed = 0.0;

    if (argc > 1)
        sscanf(argv[1], "%lf", &t_h_max_upper);
    if (argc > 2)
        if ( strcmp(argv[2],"false") == 0 ) // anything besides false is assumed true
            usezerotraj = false;
    if (argc > 3)
        sscanf(argv[3], "%d", &max_cnt);
    if (argc > 4)
        sscanf(argv[4], "%d", &max_miss);
    if (argc > 5)
        sscanf(argv[5], "%d", &max_samplemiss);
    if (argc > 6)
        sscanf(argv[6], "%d", &numbruns);
    if (argc > 7) {
        sscanf(argv[7], "%lf", &stopdist);
        if ( stopdist != 0 )
            computexxgoaldist = true;
    }
    if (argc > 8)
        sscanf(argv[8], "%d", &printskip);
    if (argc > 9) {
        sscanf(argv[9], "%lf", &nndelta);
    }
    if (argc > 10) {
        if ( strcmp(argv[10],"0") != 0 ) {
            stats_f.open(argv[10], ios::trunc);
            is_stats_out = true;
        }
    }
    if (argc > 11) {
        if ( strcmp(argv[11],"0") != 0 ) {
            traj_f.open(argv[11], ios::trunc);
            traj_f << "{" << argv[11];
            is_traj_out = true;
            save_edges = true;
        }
    }
    if (argc > 12)
        sscanf(argv[12], "%lf", &seed);
    

    // Initial Conditions
    vector<double> u0 = { 0.0 };
//    ode_state_type x0 = { M_PI, 0.0, 3.0, 0.0 };
    ode_state_type x0 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // Bounds
    sampling_bnds bnds;
    
    double bnds_delta = 10;
    bnds.xbnds << -2.0 * M_PI    , 2.0 * M_PI,
             -bnds_delta * M_PI * 0.5   , bnds_delta * M_PI * 0.5,
             -2.0 * M_PI    , 2.0 * M_PI,
             -bnds_delta * M_PI * 0.5   , bnds_delta * M_PI * 0.5,
             -2.0 * M_PI    , 2.0 * M_PI,
             -bnds_delta * M_PI * 0.5   , bnds_delta * M_PI * 0.5,
             -1.0                       , 8.0,
             -bnds_delta * 1.5          ,bnds_delta * 1.5;
    bnds.ubnds << -20.0 , 20.0;
    
    kin_constraints constraints;
    constraints.xbnds = bnds.xbnds;
    constraints.ubnds = bnds.ubnds;

    // Obstacles
    constraints.obs_radii = { 0.6, 0.6 };
    constraints.obs_Xs = { 3.0, 3.0 };
    constraints.obs_Ys = { 0.85, -0.85 };

    
    // QQ, RR, P1
    NSxNS_type QQ = NSxNS_type::Zero();
    for (int ii = 0; ii< NS; ++ii)
        QQ(ii, ii) = 1/pow((constraints.xbnds(ii,1) - constraints.xbnds(ii,0))/2,2);
    NIxNI_type invRR = NIxNI_type::Zero();
    NIxNI_type RR = NIxNI_type::Zero();
    for (int ii = 0; ii< NI; ++ii){
        invRR(ii, ii) = pow((constraints.ubnds(ii,1) - constraints.ubnds(ii,0))/2,2);
        RR(ii, ii) = 1/invRR(ii,ii);
    }
    NSxNS_type P1;
    P1 = QQ;
    


    // Root Vertex
    double tf = t_h_max_upper;
    vector<double> tt_vec;
    vector< ode_state_type > xx_vec;
    sys::integrate_free_dynamics(x0, tf, tt_vec, xx_vec);
    InterpVector xx_interp(tt_vec, xx_vec);
    vector< vector<double> > uu_vec(3,vector<double>(NI,0.0));
    tt_vec.clear();
    tt_vec.push_back(0.0);
    tt_vec.push_back(tf/2.0);
    tt_vec.push_back(tf);
    InterpVector uu_interp(tt_vec, uu_vec);

    clock_t t_start2 = clock();

    TreeVertex * init_vert = new TreeVertex(x0, usezerotraj, t_h_max_upper, &xx_interp, &uu_interp, constraints, QQ, RR, P1, QQ, RR, P1, tf, save_edges);
    Tree tree(init_vert, t_h_max_upper, bnds, max_cnt, max_miss, max_samplemiss, printskip, seed);
    
    if ( computexxgoaldist ) {
        NSx1_type xxgoal;
        xxgoal << 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 6.0 , 0.0;
        tree.set_xxgoal_stopdist (stopdist, xxgoal);
    }
    
    tree.set_nndelta(nndelta);
    
    double meancnt = 0;
    double meanmiss = 0;
    double meandist = 0;
    double meansamplemiss = 0;
    int numbxxgoalfails = 0;
    for (int ii = 0; ii < numbruns; ++ii) {
//        cout << " ########################### " << endl;
        cout << " ########################### " << endl;
        cout << "           ii = " << ii << endl;
        cout << " ########################### " << endl;
//        cout << " ########################### " << endl;
        
        tree.runrrt();
        
        meandist += tree.get_disttogoal()/numbruns;
        meancnt += double(tree.get_cnttot())/numbruns;
        meanmiss += double(tree.get_misstot())/numbruns;
        meansamplemiss += double(tree.get_samplemisstot())/numbruns;
        
        if ( tree.get_disttogoal() > stopdist )
            numbxxgoalfails++;
        
        if ( is_traj_out )
            tree.print( 0.02, "tree_"+to_string(ii), traj_f);
        
        tree.reset_tree();
    }
    cout << "mean dist is: " << meandist << endl;
    cout << "mean cnt is: " << meancnt << endl;
    cout << "mean miss is: " << meanmiss << endl;
    cout << "mean sample miss is: " << meansamplemiss << endl;
    cout << "mean run time: " << ((float)(clock()-t_start2))/CLOCKS_PER_SEC/numbruns << " seconds." << endl;// USE THIS FOR TIMING
    cout << "number of xxgoal fails: " << numbxxgoalfails << endl;
    

    if ( is_stats_out ) {
        stats_f << "meandist, " << meandist << endl;
        stats_f << "meancnt, " << meancnt << endl;
        stats_f << "meanmiss, " << meanmiss << endl;
        stats_f << "meansamplemiss, " << meansamplemiss << endl;
        stats_f << "meanruntime, " << ((float)(clock()-t_start2))/CLOCKS_PER_SEC/numbruns << endl;
        stats_f << "numbxxgoalsfails, " << numbxxgoalfails << endl;
        stats_f.close();
    }
    
    if ( is_traj_out ) {
        traj_f << "}";
        traj_f.close();
    }

}
