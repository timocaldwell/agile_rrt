//
//  interp.h
//  Agile_RRT
//
//  Created by Timothy Caldwell on 9/20/14.
//  Copyright (c) 2014 Timothy Caldwell. All rights reserved.
//

#ifndef __Agile_RRT__interp__
#define __Agile_RRT__interp__

#include "global.h"

class Interp {
 public:
  Interp();
  Interp(const vector<double> & xx_vec, const vector<double> & yy_vec);
  Interp(const Interp & other);
  ~Interp();
  void set(const vector<double> & xx_vec, const vector<double> & yy_vec);
  double pt(double xx) const; // Note extrapolation holds to left or right bound.
  double d_pt(double xx) const {  // Derivative at pt
    return gsl_spline_eval_deriv(spline_, xx, acc_);
  }
  double dd_pt(double xx) const { // 2nd Derivative at pt
    return gsl_spline_eval_deriv2(spline_, xx, acc_);
  }
  double begin() const {return xx_left_;};
  double end() const {return xx_right_;};
  void Print(double dx) const;
  
 private:
  gsl_interp_accel * acc_;
  gsl_spline * spline_;
  vector<double> xxraw_vec_;
  vector<double> yyraw_vec_;
  double xx_left_, xx_right_;
};

class InterpVector {
 public:
  InterpVector();
  InterpVector(const vector<double> & tt_vec, const vector<vector<double>> & xx_vec);
  InterpVector(const InterpVector & other);
  ~InterpVector();
  void set(const vector<double> & tt_vec, const vector<vector<double>> & xx_vec);
  void pt(double tt, vector<double> * xx_pt) const;
  void Print(double dt, ostream & stream = cout) const;
  void Print(double dt, const string & name, ostream & stream) const;

 private:
  vector< Interp* > xx_;
  size_t xx_pt_size_;
};

#endif // __Agile_RRT__interp__