//
//  global.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 9/20/14.
//  Copyright (c) 2014 Timothy Caldwell. All rights reserved.
//

#include "global.h"

const size_t NI = 1;
const size_t NS = 8;
const size_t numlinks = 3;

//const double t_h_MAX_UPPER = .1;//1.0;
//const bool usezerotraj = true;//false; //true;


const double LARGENUM = 987654321;
const double SMALLNUM = 1/LARGENUM;

//const double gravity = 9.81, headmass = 0.1, cartmass = 1.0, pendlength = 1.0;

const int MAX_NUM_INTEGRATION_STEPS = 10000;
const double INTEGRATION_DT_START = 0.001;

// control types for state integration
const int FREEDYNAMICS_TYPE = 0;
const int FEEDFORWARD_TYPE = 1;
const int TRACKPOINT_TYPE = 2;
const int TRACKTRAJECTORY_TYPE = 3;
const int LINEARSTEERINGPROJECTION_TYPE = 4;
const int LINEARSTEERINGPROJECTIONWCOST_TYPE = 5;

const int LTI_TYPE = 0;
const int LTV_TYPE = 1;

const int LTI_OPENLOOP_TYPE = 0;
const int LTI_CLOSEDLOOP_TYPE = 1;
const int LTV_OPENLOOP_TYPE = 2;
const int LTV_CLOSEDLOOP_TYPE = 3;

const int FREE_CONTROLTYPE = 0;
const int FEEDFORWARD_CONTROLTYPE = 1;
const int MINENERGY_OPENLOOP_CONTROLTYPE = 2;
const int MINENERGY_CLOSEDLOOP_CONTROLTYPE = 2;


const int WW_GRAMIANTYPE = 0;
const int WWK_GRAMIANTYPE = 1;
const int SSK_GRAMIANTYPE = 2;

//struct kin_constraints {
//    Eigen::Matrix<double, 4, 2> xbnds;
//    Eigen::Matrix<double, 1, 2> ubnds;
//
//    vector< double > obs_radii;
//    vector< double > obs_Xs;
//    vector< double > obs_Ys;
//};


void print_vec(const vector<double> & vec, ostream & stream){
    for(vector<double>::const_iterator it = vec.begin(); it != vec.end()-1; ++it)
        stream << std::setprecision(6) << *it << ", ";
    stream << vec.back();
}
void print_vec(const vector<double> & vec, const string & name, ostream & stream){
    stream << ",{" << name << ",{";
    print_vec(vec, stream);
    stream << "}}";
}

