//
//  main.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell.
//  Copyright (c) 2015 Timothy Caldwell. Refer to license.txt
//
//  AgileRRT motion plans the highly dynamic 3-link pendulum on a cart. It executes the rapidly exploring random tree using LQR-based steering with linearization about the zero-control trajectory.
//
//  Arguments: (default value given right of '=' sign)
//    t_h_max_upper = 1.0       | Max time horizon.
//    usezerotraj = true        | Choose linearization for steering: true -> Linearization about zero-control trajectory.
//                                                                   false -> Linearization about vertex x0.
//    max_cnt = 1000            | Max number of successful insertions into the tree. RRT ends if max_cnt reached.
//    max_miss = 1000           | Max number of failed insertions into the tree. RRT ends if max_miss reached.
//    max_samplemiss = 10000    | Max number of failed attempts to find a
//    numbruns = 1              | Number of runs of the RRT. Run more to get more statistics.
//    stopdist = 2.0            | Stops RRT if a point is explored within a ball of radius stopdist from the goal state.
//    printskip = 1             | Number of loops of the RRT to skip between printing details to standard out.
//    nndelta = LARGENUM        | Nearest Neighbor ball radius for quickly checking
//    stats_name = "0"          | Filename to store statistics as a csv. Use "0" or no arg to not make file.
//    traj_name = "0"           | Filename to store tree trajectories easily readable by Mathematica. Use "0" or no arg to not make file
//    seed = 0                  | To choose the seed to the random number generator



#include "global.h"
#include "interp.h"
#include "system.h"
#include "pendcart_3link.h"
#include "treevertex.h"
#include "tree.h"
//#include "rrtreach.h"

// Struct for all options given as inputs through arguments.
struct inputs_struct {
  double t_h_max_upper = 1.0;
  bool usezerotraj = true;
  int max_cnt = 1000;
  int max_miss = 1000;
  int max_samplemiss = 10000;
  int numbruns = 1;
  double stopdist = 0.0;
  bool computexxgoaldist = false;
  int printskip = 1;
  double nndelta = kLARGENUM;
  ofstream stats_f, traj_f;
  bool is_stats_out = false;
  bool is_traj_out = false;
  bool save_edges = false;
  double seed = 0.0;
};

void SpecifyInputs(int argc, const char * argv[], inputs_struct * inputs) {
  if (argc > 1)
    sscanf(argv[1], "%lf", &inputs->t_h_max_upper);
  if (argc > 2)
    if ( strcmp(argv[2],"false") == 0 ) // anything besides false is assumed true
        inputs->usezerotraj = false;
  if (argc > 3)
    sscanf(argv[3], "%d", &inputs->max_cnt);
  if (argc > 4)
    sscanf(argv[4], "%d", &inputs->max_miss);
  if (argc > 5)
    sscanf(argv[5], "%d", &inputs->max_samplemiss);
  if (argc > 6)
    sscanf(argv[6], "%d", &inputs->numbruns);
  if (argc > 7) {
    sscanf(argv[7], "%lf", &inputs->stopdist);
    if (inputs->stopdist != 0)
      inputs->computexxgoaldist = true;
  }
  if (argc > 8)
    sscanf(argv[8], "%d", &inputs->printskip);
  if (argc > 9)
    sscanf(argv[9], "%lf", &inputs->nndelta);
  if (argc > 10) {
    if ( strcmp(argv[10],"0") != 0 ) {
      inputs->stats_f.open(argv[10], ios::trunc);
      inputs->is_stats_out = true;
    }
  }
  if (argc > 11) {
    if ( strcmp(argv[11],"0") != 0 ) {
      inputs->traj_f.open(argv[11], ios::trunc);
      inputs->traj_f << "{" << argv[11];
      inputs->is_traj_out = true;
      inputs->save_edges = true;
    }
  }
  if (argc > 12)
    sscanf(argv[12], "%lf", &inputs->seed);
}

// Initializes the system state and control bounds and circular obstacle constraints.
void SpecifyConstraints(sampling_bounds_struct * bounds, constraints_struct * constraints) {
  // Sets upper and lower boundaries of state and control
  double bounds_delta = 10;
  bounds->xbounds << -2.0 * M_PI    , 2.0 * M_PI,
                -bounds_delta * M_PI * 0.5   , bounds_delta * M_PI * 0.5,
                -2.0 * M_PI    , 2.0 * M_PI,
                -bounds_delta * M_PI * 0.5   , bounds_delta * M_PI * 0.5,
                -2.0 * M_PI    , 2.0 * M_PI,
                -bounds_delta * M_PI * 0.5   , bounds_delta * M_PI * 0.5,
                -1.0                       , 8.0,
                -bounds_delta * 1.5          ,bounds_delta * 1.5;
  bounds->ubounds << -20.0 , 20.0;
  
  constraints->xbounds = bounds->xbounds;
  constraints->ubounds = bounds->ubounds;

  // Circular obstacle details
  constraints->obs_radii = { 0.6, 0.6 };
  constraints->obs_Xs = { 3.0, 3.0 };
  constraints->obs_Ys = { 0.85, -0.85 };
}

// Initializes quadratic cost function gains QQ, RR, P1 based on the specified state and control bounds.
// Cost function has form (using Latex syntax) J = 0.5 \int_0^{th}xx^T*QQ*xx + uu^T*RR*uu dt + 0.5*xx(th)^T*P1*xx(th).
// RRinv is the matrix inverse of RR.
void SpecifyCostGains(const sampling_bounds_struct & bounds, NSxNS_type * QQ, NIxNI_type * RR, NIxNI_type * RRinv, NSxNS_type * P1) {
  *QQ = NSxNS_type::Zero();
  for (int ii = 0; ii< NS; ++ii)
    (*QQ)(ii, ii) = 1/pow((bounds.xbounds(ii,1) - bounds.xbounds(ii,0))/2,2);
  *RRinv = NIxNI_type::Zero();
  *RR = NIxNI_type::Zero();
  for (int ii = 0; ii< NI; ++ii) {
    (*RRinv)(ii, ii) = pow((bounds.ubounds(ii,1) - bounds.ubounds(ii,0))/2,2);
    (*RR)(ii, ii) = 1/(*RRinv)(ii,ii);
  }
  *P1 = *QQ;
}

void SpecifyRootTrajectory(const ode_state_type & x0, double t_h_max_upper, InterpVector * xx_interp, InterpVector * uu_interp) {
  vector<double> tt_vec;
  vector<ode_state_type> xx_vec;
  
  // Initializing to free dynamics
  sys::IntegrateFreeDynamics(x0, t_h_max_upper, &tt_vec, &xx_vec);
  xx_interp = new InterpVector(tt_vec, xx_vec);
  
  // Control input is set to 0.0
  vector< vector<double> > uu_vec(3,vector<double>(NI,0.0));
  tt_vec.clear();
  tt_vec.push_back(0.0);
  tt_vec.push_back(t_h_max_upper/2.0);
  tt_vec.push_back(t_h_max_upper);
  uu_interp = new InterpVector(tt_vec, uu_vec);
}


// ARGUMENTS : [t_h_max_upper = 1.0] , [ usezerotraj = true] , [max_cnt = 1000], [max_miss = 1000], [max_samplemiss = 10000], [numbruns = 1] , [stopdist = 2.0] , [printskip = 1] , [nndelta = LARGENUM], [stats_name = "0" (use "0" or no arg to not make file)] , [traj_name = "0" (use "0" or no arg to not make file)] , [seed = 0]
int main(int argc, const char * argv[])
{
  inputs_struct inputs;
  SpecifyInputs(argc, argv, &inputs);
  
  // The initial conditions
  ode_state_type x0 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  // Set bounds and constraints.
  sampling_bounds_struct bounds;
  constraints_struct constraints;
  SpecifyConstraints(&bounds, &constraints);
  
  // Set cost gains.
  NSxNS_type QQ, P1;
  NIxNI_type RR, RRinv;
  SpecifyCostGains(bounds, &QQ, &RR, &RRinv, &P1);

  // Set root trajectory to free dynamics.
  InterpVector * xx_interp = nullptr;
  InterpVector * uu_interp = nullptr;
  SpecifyRootTrajectory(x0, inputs.t_h_max_upper, xx_interp, uu_interp);
  
  // Create initial vertex.
  TreeVertex * init_vert = new TreeVertex(x0, inputs.usezerotraj, inputs.t_h_max_upper, *xx_interp, *uu_interp, constraints, QQ, RR, P1, QQ, RR, P1, inputs.t_h_max_upper, inputs.save_edges);

  clock_t t_start2 = clock();

//  Execute( )
  Tree tree(init_vert, inputs.t_h_max_upper, bounds, inputs.max_cnt, inputs.max_miss, inputs.max_samplemiss, inputs.printskip, inputs.seed);

  if (inputs.computexxgoaldist) {
    NSx1_type xxgoal;
    xxgoal << 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 6.0 , 0.0;
    tree.set_xxgoal_stopdist (inputs.stopdist, xxgoal);
  }
  
  tree.set_nndelta(inputs.nndelta);
  
  double meancnt = 0;
  double meanmiss = 0;
  double meandist = 0;
  double meansamplemiss = 0;
  int numbxxgoalfails = 0;
  for (int ii = 0; ii < inputs.numbruns; ++ii) {
    cout << " ########################### " << endl;
    cout << "           ii = " << ii << endl;
    cout << " ########################### " << endl;
    
    tree.RunRRT();
    
    meandist += tree.get_disttogoal()/inputs.numbruns;
    meancnt += double(tree.get_cnttot())/inputs.numbruns;
    meanmiss += double(tree.get_misstot())/inputs.numbruns;
    meansamplemiss += double(tree.get_samplemisstot())/inputs.numbruns;
    
    if ( tree.get_disttogoal() > inputs.stopdist )
      numbxxgoalfails++;
    
    if ( inputs.is_traj_out )
      tree.Print( 0.02, "tree_"+to_string(ii), &inputs.traj_f);
    
    tree.ResetTree();
  }
  cout << "mean dist is: " << meandist << endl;
  cout << "mean cnt is: " << meancnt << endl;
  cout << "mean miss is: " << meanmiss << endl;
  cout << "mean sample miss is: " << meansamplemiss << endl;
  cout << "mean run time: " << ((float)(clock()-t_start2))/CLOCKS_PER_SEC/inputs.numbruns << " seconds." << endl;
  cout << "number of xxgoal fails: " << numbxxgoalfails << endl;
  
  if ( inputs.is_stats_out ) {
    inputs.stats_f << "meandist, " << meandist << endl;
    inputs.stats_f << "meancnt, " << meancnt << endl;
    inputs.stats_f << "meanmiss, " << meanmiss << endl;
    inputs.stats_f << "meansamplemiss, " << meansamplemiss << endl;
    inputs.stats_f << "meanruntime, " << ((float)(clock()-t_start2))/CLOCKS_PER_SEC/inputs.numbruns << endl;
    inputs.stats_f << "numbxxgoalsfails, " << numbxxgoalfails << endl;
    inputs.stats_f.close();
  }
  if ( inputs.is_traj_out ) {
    inputs.traj_f << "}";
    inputs.traj_f.close();
  }
}
