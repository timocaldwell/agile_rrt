//
//  interp.cpp
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//

#include "interp.h"

Interp::Interp(){
  tt_left_ = 0;
  tt_right_ = 0;
}
Interp::Interp(const vector<double> & tt_vec, const vector<double> & xx_vec){
  set(tt_vec, xx_vec);
}
// Copies by calling set() on the initially set points.
Interp::Interp(const Interp & other){
  set(other.ttraw_vec_, other.xxraw_vec_);
}
Interp::~Interp(){
  gsl_spline_free(spline_);
  gsl_interp_accel_free(acc_);
}
void Interp::set(const vector<double> & tt_vec, const vector<double> & xx_vec) {
  ttraw_vec_ = tt_vec;
  xxraw_vec_ = xx_vec;
  if (tt_vec.size() <= 1) {
    cout << " ERROR --- Interp input vectors' size too small" << endl;
    abort();
  }
  tt_left_ = tt_vec[0];
  tt_right_ = tt_vec.back();
  
  // GSL interpolation needs 3 or more points. When two are given, a 3rd is inserted directly between the two given.
  if (ttraw_vec_.size() == 2) {
    ttraw_vec_.insert(ttraw_vec_.begin()+1, (tt_right_ + tt_left_) * 0.5);
    xxraw_vec_.insert(xxraw_vec_.begin()+1, (xxraw_vec_[0] + xxraw_vec_.back()) * 0.5);
  }

  // Initialize acc_ and spline_ for GSL interpolation
  acc_ = gsl_interp_accel_alloc();
  spline_ = gsl_spline_alloc(gsl_interp_cspline, ttraw_vec_.size());
  
  gsl_spline_init(spline_, &ttraw_vec_[0], &xxraw_vec_[0], ttraw_vec_.size());
}
double Interp::pt(double tt) const {
  // extrapolation holds to left or right bound.
  if (tt > tt_right_)
    return gsl_spline_eval(spline_, tt_right_, acc_);
  else if (tt < tt_left_)
    return gsl_spline_eval(spline_, tt_left_, acc_);
  return gsl_spline_eval(spline_, tt, acc_);
}
void Interp::Print(double dt, ostream * stream) const {
  double tt = tt_left_;
  (*stream) << "{";
  for (int ii = 0; ii < 1000; ii++ ) {
    (*stream) << std::setprecision(kPRINT_PRECISION) << "{" << tt << "," << pt(tt) << "}";
    tt += dt;
    if (tt <= tt_right_)
      (*stream) << ",";
    else {
      (*stream) << "}" << endl;
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
  
  // copies each Interp object individually.
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
// Destructs each Interp object.
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
void InterpVector::Print(double dt, ostream * stream) const {
  double tt;
  if (xx_pt_size_ == 0) {
    (*stream) << ",{}";
    return;
  }
  (*stream) << ",{";
  for (unsigned jj = 0; jj < xx_pt_size_; ++jj) {
    tt = xx_[jj]->begin();
    (*stream) << "{";
    for (int ii = 0; ii < 1000; ++ii ){
      if ( xx_[jj]->pt(tt) > kPRINT_CHOP_BOUND || xx_[jj]->pt(tt) < -kPRINT_CHOP_BOUND ) // Chop for printing.
        (*stream) << "{" << tt << "," << std::setprecision(kPRINT_PRECISION) <<  xx_[jj]->pt(tt) << "}";
      else
        (*stream) << "{" << tt << "," << 0 << "}";
      tt += dt;
      (*stream) << ",";
      
      if (tt >= xx_[jj]->end()) {
        // Print final time.
        tt = xx_[jj]->end();
        if ( xx_[jj]->pt(tt) > kPRINT_CHOP_BOUND || xx_[jj]->pt(tt) < -kPRINT_CHOP_BOUND ) // Chop for printing.
          (*stream) << "{" << tt << "," << std::setprecision(kPRINT_PRECISION) <<  xx_[jj]->pt(tt) << "}";
        else
          (*stream) << "{" << tt << "," << 0 << "}";
        
        (*stream) << "}";
        break;
      }
    }
    if (jj < xx_pt_size_ - 1)
      (*stream) << ",";
    else
      (*stream) << "}" << endl;
  }
}
void InterpVector::Print(double dt, const string & name, ostream * stream) const {
  (*stream) << ",{" << name ;
  Print(dt, stream);
  (*stream) << "}";
}
void InterpVector::PrintEndPt(ostream * stream) const {
  (*stream) << "{{" << xx_[0]->end() << "},{";
  double tt = xx_[0]->end();
  for (unsigned jj = 0; jj < xx_pt_size_; ++jj) {
    if ( xx_[jj]->pt(tt) > kPRINT_CHOP_BOUND || xx_[jj]->pt(tt) < -kPRINT_CHOP_BOUND ) // Chop for printing.
      (*stream) << std::setprecision(kPRINT_PRECISION) <<  xx_[jj]->pt(tt);
    else
      (*stream) << 0;
    if (jj < xx_pt_size_ - 1)
      (*stream) << ",";
  }
  (*stream) << "}}";
}
