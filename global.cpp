//
//  global.cpp
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//

#include "global.h"

// see global.h for comments
const int NI = 1;
const int NS = 8;

const double kLARGENUM = 987654321;
const int kMAX_NUM_INTEGRATION_STEPS = 10000;
const double kINTEGRATION_DT_START = 0.001;

const int kFREEDYNAMICS_TYPE = 0;
const int kFEEDFORWARD_TYPE = 1;
const int kTRACKPOINT_TYPE = 2;
const int kTRACKTRAJECTORY_TYPE = 3;
const int kLINEARSTEERINGPROJECTION_TYPE = 4;
const int kLINEARSTEERINGPROJECTIONWCOST_TYPE = 5;

const int kLTI_TYPE = 0;
const int kLTV_TYPE = 1;

const int kLTI_OPENLOOP_TYPE = 0;
const int kLTI_CLOSEDLOOP_TYPE = 1;
const int kLTV_OPENLOOP_TYPE = 2;
const int kLTV_CLOSEDLOOP_TYPE = 3;

const int kFREE_CONTROLTYPE = 0;
const int kFEEDFORWARD_CONTROLTYPE = 1;
const int kMINENERGY_OPENLOOP_CONTROLTYPE = 2;
const int kMINENERGY_CLOSEDLOOP_CONTROLTYPE = 2;

const int kWW_GRAMIANTYPE = 0;
const int kWWK_GRAMIANTYPE = 1;
const int kSSK_GRAMIANTYPE = 2;

const int kPRINT_PRECISION = 6;
const double kPRINT_CHOP_BOUND = .0001;

void PrintVec(const vector<double> & vec, ostream * stream) {
  for(vector<double>::const_iterator it = vec.begin(); it != vec.end()-1; ++it)
    (*stream) << std::setprecision(kPRINT_PRECISION) << *it << ", ";
  (*stream) << vec.back();
}
void PrintVec(const vector<double> & vec, const string & name, ostream * stream) {
  (*stream) << ",{" << name << ",{";
  PrintVec(vec, stream);
  (*stream) << "}}";
}

