//
//  global.h
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
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

extern const int NI;                              // Constant; Number of inputs to the dynamic system
extern const int NS;                              // Constant; Number of states to the dynamic system

extern const double kLARGENUM;                    // Constant; A sufficiently large number

// Constants; For integration using odeint
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

typedef vector<double> ode_state_type;            // odeint integrator state_type

///////////////////////////////
//Currently implemented as global. Needs to be updated for NS and NI
///////////////////////////////
// Eigen Matrix types of size NS or NI by NS or NI
typedef Eigen::Matrix< double , 8 , 8 > NSxNS_type;
typedef Eigen::Matrix< double , 8 , 1 > NSxNI_type;
typedef Eigen::Matrix< double , 1 , 1 > NIxNI_type;
typedef Eigen::Matrix< double , 1 , 8 > NIxNS_type;
typedef Eigen::Matrix< double , 8 , 1 > NSx1_type;
typedef Eigen::Matrix< double , 1 , 1 > NIx1_type;

struct constraints_struct {
  Eigen::Matrix<double, 8, 2> xbounds;              // list of upper and lower bounds for each state
  Eigen::Matrix<double, 1, 2> ubounds;              // list of upper and lower bounds for each input

  // list of circle obstacles of radius obs_radii[.] at point (obs_Xs[.], obs_Ys[.])
  vector< double > obs_radii;
  vector< double > obs_Xs;
  vector< double > obs_Ys;
};
// list of upper and lower bounds for each state and input. Should be the same as those set in constraints_struct
struct sampling_bounds_struct {
  Eigen::Matrix<double, 8, 2> xbounds;
  Eigen::Matrix<double, 1, 2> ubounds;
};

// stepper type used by odeint for integration
//typedef odeint::controlled_runge_kutta< odeint::runge_kutta_dopri5<ode_state_type> > stepper_type; // Alternative stepper choice
typedef odeint::controlled_runge_kutta< odeint::runge_kutta_cash_karp54<ode_state_type> > stepper_type;

// Printing parameters
extern const int kPRINT_PRECISION;                // Printing precision number
extern const double kPRINT_CHOP_BOUND;            // Lower bound for chopping (i.e. if abs(x)<kPRINT_CHOP_BOUND, then 0 is returned)

// prints a vec to stream in a form easily readable by Mathematica.
extern void PrintVec(const vector<double> & vec, ostream * stream = &(std::cout));
extern void PrintVec(const vector<double> & vec, const string & name, ostream * stream);


#endif // __Agile_RRT__global__
