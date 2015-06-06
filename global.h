//
//  global.h
//  Agile_RRT
//
//  Created by Timothy Caldwell on 9/20/14.
//  Copyright (c) 2014 Timothy Caldwell. All rights reserved.
//

#ifndef __Agile_RRT__global__
#define __Agile_RRT__global__

#include <fstream>
//#include <string>
//#include <sstream>

#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <math.h>
#include <random>

//boost libraries
#include <boost/numeric/odeint.hpp>

//gsl libraries
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

//Eigen libraries
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

//namespaces
using namespace std;
namespace odeint = boost::numeric::odeint;

extern const int NI;
extern const int NS;
extern const int numlinks;

extern const double kLARGENUM;
extern const double kSMALLNUM;
extern const int kMAX_NUM_INTEGRATION_STEPS;
extern const double kINTEGRATION_DT_START;

//integrate state types
extern const int kFREEDYNAMICS_TYPE;
extern const int kFEEDFORWARD_TYPE;
extern const int kTRACKPOINT_TYPE;
extern const int kTRACKTRAJECTORY_TYPE;
extern const int kLINEARSTEERINGPROJECTION_TYPE;
extern const int kLINEARSTEERINGPROJECTIONWCOST_TYPE;

//integrate linear state and ricatti types
extern const int kLTI_TYPE;
extern const int kLTV_TYPE;

//STMSlot2 types
extern const int kLTI_OPENLOOP_TYPE;
extern const int kLTI_CLOSEDLOOP_TYPE;
extern const int kLTV_OPENLOOP_TYPE;
extern const int kLTV_CLOSEDLOOP_TYPE;

//control type for linear state equation
extern const int kFREE_CONTROLTYPE;
extern const int kFEEDFORWARD_CONTROLTYPE;
extern const int kMINENERGY_OPENLOOP_CONTROLTYPE;
extern const int kMINENERGY_CLOSEDLOOP_CONTROLTYPE;

//Gramian types
extern const int kWW_GRAMIANTYPE;
extern const int kWWK_GRAMIANTYPE;
extern const int kSSK_GRAMIANTYPE;

typedef vector<double> ode_state_type; // odeint integrator state_type

///////////////////////////////
//Currently implemented as global. Needs to be updated for NS and NI
///////////////////////////////
typedef Eigen::Matrix< double , 8 , 8 > NSxNS_type;
typedef Eigen::Matrix< double , 8 , 1 > NSxNI_type;
typedef Eigen::Matrix< double , 1 , 1 > NIxNI_type;
typedef Eigen::Matrix< double , 1 , 8 > NIxNS_type;
typedef Eigen::Matrix< double , 8 , 1 > NSx1_type;
typedef Eigen::Matrix< double , 1 , 1 > NIx1_type;
struct kin_constraints {
  Eigen::Matrix<double, 8, 2> xbnds;
  Eigen::Matrix<double, 1, 2> ubnds;

  vector< double > obs_radii;
  vector< double > obs_Xs;
  vector< double > obs_Ys;
};
struct sampling_bnds {
  Eigen::Matrix<double, 8, 2> xbnds;
  Eigen::Matrix<double, 1, 2> ubnds;
};


//typedef odeint::adams_bashforth_moulton< 5, double > stepper_type;
//typedef odeint::runge_kutta_dopri5< double > stepper_type;
//typedef odeint::runge_kutta_dopri5< double[8] > stepper_type;
//typedef odeint::controlled_runge_kutta< odeint::runge_kutta_dopri5<ode_state_type> > stepper_type;
typedef odeint::controlled_runge_kutta< odeint::runge_kutta_cash_karp54< ode_state_type > > stepper_type;

extern void PrintVec(const vector<double> & vec, ostream & stream = std::cout);
extern void PrintVec(const vector<double> & vec, const string & name, ostream & stream);


#endif // __Agile_RRT__global__
