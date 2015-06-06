//
//  interp.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 9/20/14.
//  Copyright (c) 2014 Timothy Caldwell. All rights reserved.
//

#include "interp.h"

Interp::Interp(){
  xx_left_ = 0;
  xx_right_ = 0;
}
Interp::Interp(const vector<double> & xx_vec, const vector<double> & yy_vec){
  set(xx_vec, yy_vec);
}
Interp::Interp(const Interp & other){ //copy constructor
  set(other.xxraw_vec_, other.yyraw_vec_);
}
Interp::~Interp(){
  gsl_spline_free(spline_);
  gsl_interp_accel_free(acc_);
}
void Interp::set(const vector<double> & xx_vec, const vector<double> & yy_vec){
  xxraw_vec_ = xx_vec;
  yyraw_vec_ = yy_vec;

  if (xx_vec.size() <= 1) {
    cout << " ERROR --- Interp input vectors' size too small" << endl;
    abort();
  }

  xx_left_ = xx_vec[0];
  xx_right_ = xx_vec.back();

  if (xxraw_vec_.size() == 2) { // inserts a 3rd point through linear interpolation
    xxraw_vec_.insert(xxraw_vec_.begin()+1, (xx_right_ + xx_left_) * 0.5);
    yyraw_vec_.insert(yyraw_vec_.begin()+1, (yyraw_vec_[0] + yyraw_vec_.back()) * 0.5);
  }

  acc_ = gsl_interp_accel_alloc();
  spline_ = gsl_spline_alloc(gsl_interp_cspline, xxraw_vec_.size());
  
  gsl_spline_init(spline_, &xxraw_vec_[0], &yyraw_vec_[0], xxraw_vec_.size());
}

double Interp::pt(double xx) const {
// extrapolation holds to left or right bound.
  if (xx > xx_right_)
    return gsl_spline_eval(spline_, xx_right_, acc_);
  else if (xx < xx_left_)
    return gsl_spline_eval(spline_, xx_left_, acc_);
  return gsl_spline_eval(spline_, xx, acc_);
}

void Interp::Print(double dx) const {
  double xx = xx_left_;
  cout << "{";
  for (int ii = 0; ii < 1000; ii++ ) {
    cout << std::setprecision(6) << "{" << xx << "," << pt(xx) << "}";
    xx += dx;
    if (xx <= xx_right_)
      cout << ",";
    else {
      cout << "}" << endl;
      break;
    }
  }
}


InterpVector::InterpVector() {
  xx_pt_size_ = 0;
}

InterpVector::InterpVector(const vector<double> & tt_vec, const vector<vector<double>> & xx_vec) {
  set(tt_vec, xx_vec);
}

InterpVector::InterpVector(const InterpVector& other) { // copy constructor
  xx_pt_size_ = other.xx_pt_size_;
  xx_.clear();
  
  Interp * xx_copy;
  for (vector< Interp* >::const_iterator it = other.xx_.begin(); it!=other.xx_.end(); ++it) {
    xx_copy = new Interp(**it);
    xx_.push_back(xx_copy);
  }
}

void InterpVector::set(const vector<double> & tt_vec, const vector<vector<double>> & xx_vec) {
  Interp * interp_temp;
  xx_pt_size_ = xx_vec[0].size();
  for (int ii = 0; ii < xx_pt_size_; ++ii)
  {
    vector<double> xx_temp;
    for (vector<vector<double>>::const_iterator it = xx_vec.begin(); it!=xx_vec.end(); ++it)
        xx_temp.push_back((*it)[ii]);
    interp_temp = new Interp(tt_vec , xx_temp);
    xx_.push_back(interp_temp);
  }
}

InterpVector::~InterpVector() {
  if(xx_.size()!=0){
    for(unsigned ii = 0; ii<xx_.size(); ++ii)
      delete xx_.at(ii);
    xx_.clear();
  }
}

void InterpVector::pt(double tt, vector<double> * xx_pt) const {
  xx_pt->clear();
  for(unsigned ii = 0; ii < xx_pt_size_; ++ii)
    xx_pt->push_back(xx_[ii]->pt(tt));
}


//void InterpVector::print(double dt, ostream & stream){
//    double tt;
//    cout << "{";
//    for (unsigned jj = 0; jj < xx_pt_size; ++jj){
//        tt = xx[0]->begin();
//        cout << "{";
//        for (int ii = 0; ii < 1000; ++ii ){
//            cout << "{" << tt << "," << xx[jj]->pt(tt) << "}";
//            tt += dt;
//            if (tt <= xx[0]->end())
//                cout << ",";
//            else {
//                cout << "}";
//                break;
//            }
//        }
//        if (jj < xx_pt_size - 1)
//            cout << ",";
//        else
//            cout << "}" << endl;
//    }
//}

void InterpVector::Print(double dt, ostream & stream) const {
  double tt;
  if (xx_pt_size_ == 0) {
    stream << ",{}";
    return;
  }
  stream << ",{";
  for (unsigned jj = 0; jj < xx_pt_size_; ++jj) {
    tt = xx_[0]->begin();
    stream << "{";
    for (int ii = 0; ii < 1000; ++ii ){
      if ( xx_[jj]->pt(tt) > .0001 || xx_[jj]->pt(tt) < -.0001) // chop for printing
        stream << "{" << tt << "," << std::setprecision(6) <<  xx_[jj]->pt(tt) << "}";
      else
        stream << "{" << tt << "," << std::setprecision(6) <<  0 << "}";
      tt += dt;
      if (tt <= xx_[0]->end())
        stream << ",";
      else {
        stream << "}";
        break;
      }
    }
    if (jj < xx_pt_size_ - 1)
      stream << ",";
    else
      stream << "}" << endl;
  }
}

void InterpVector::Print(double dt, const string & name, ostream & stream) const {
  stream << ",{" << name ;
  Print(dt, stream);
  stream << "}";
}
