//
//  interp.h
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//
//  Within are classes Interp and InterpVector.
//  Interp uses GSL interpolation to create a scalar valued interpolation function xx(tt) using cubic splines. The interpolation function is constructed from a list of points (tt_vec[.], xx_vec[.]) for the strictly monotonically increasing vector tt_vec and values xx_vec. The interpolated value xx(tt) for any double tt in the range [tt_vec[0], tt_vec.back()] is the point pt(tt). If tt is outside this range, the extrapolation holds to the boundaries values.
//  InterpVector is a vector of Interp objects. It is constructed from the vector of scalars tt_vec and the vector of vectors xx_vec. It is assumed that each xx_vec[ii] has equal size. The interpolation function xx(tt) returns a vector of this size. The interpolated value xx(tt) is xx_pt passed by pointer from pt(tt, xx_pt).


#ifndef __Agile_RRT__interp__
#define __Agile_RRT__interp__

#include "global.h"

class Interp {
 public:
  // Basic constructor. Public method set() must be called before using object constructed this way.
  Interp();
  // Constructs an interpolation function from the points (tt_vec[.], xx_vec[.]).
  Interp(const vector<double> & tt_vec, const vector<double> & xx_vec);
  Interp(const Interp & other);
  // frees acc_ and spline_
  ~Interp();
  // Sets up the interpolation function by populating private members for the points (tt_vec[.], xx_vec[.]).
  // tt_vec must be strictly monotonically increasing and tt_vec and xx_vec must be the same size greater than 1.
  void set(const vector<double> & tt_vec, const vector<double> & xx_vec);
  // Returns the interpolated value xx(tt).
  // Extrapolation holds to left or right bound.
  double pt(double tt) const;
  // Returns first derivative at tt.
  double d_pt(double tt) const {
    return gsl_spline_eval_deriv(spline_, tt, acc_);
  }
  // Returns second derivative at tt.
  double dd_pt(double tt) const {
    return gsl_spline_eval_deriv2(spline_, tt, acc_);
  }
  double begin() const {return tt_left_;};
  double end() const {return tt_right_;};
  // Prints to stream a list of interpolated points in a form easily readable by Mathematica
  //   starting from t0 = tt_left_, the output is
  //     {{t0,xx(t0)},{tt+dt,xx(tt+dt)},...,{t0+ii*dt,xx(t0+ii*dt)}}
  //   where ii is such that tt+ii*dt < t_right_ and tt+(ii+1)*dt > t_right_.
  void Print(double dt, ostream * stream = (&cout)) const;
  
 private:
  // Needed for GSL interpolation.
  gsl_interp_accel * acc_;
  gsl_spline * spline_;
  
  // The raw points initially set. They are used for copying.
  vector<double> ttraw_vec_;
  vector<double> xxraw_vec_;
  
  // The left and right bounds of tt_vec
  double tt_left_, tt_right_;
};

class InterpVector {
 public:
  // Member methods mirror those in Interp except where xx_vec is a vector of vectors such that xx_vec[.] are vectors of the equal size. The interplation function xx(tt) for double tt is called by pt(tt, xx_pt) which passes by pointer the interpolated vector xx_pt where xx_pt has the same size as each xx_vec[.].
  InterpVector();
  InterpVector(const vector<double> & tt_vec, const vector<vector<double>> & xx_vec);
  InterpVector(const InterpVector & other);
  ~InterpVector();
  void set(const vector<double> & tt_vec, const vector<vector<double>> & xx_vec);
  void pt(double tt, vector<double> * xx_pt) const;
  
  // Prints to stream a list of interpolated points in a form easily readable by Mathematica
  //   starting from t0 = tt_left_, the output is
  //    ,{{{t0,xx0(t0)},{t0,xx1(t0)},...,{t0,xxn(t0)}}
  //     {{{t0+dt,xx0(t0+dt)},{t0+dt,xx1(t0+dt)},...,{t0+dt,xxn(t0+dt)}}
  //                          .
  //                          .
  //                          .
  //     {{{t0+ii*dt,xx0(t0+ii*dt)},{t0+ii*dt,xx1(t0+ii*dt)},...,{t0+ii*dt,xxn(t0+ii*dt)}}}
  //   where n+1 is the size of the interpolated vector and ii is such that tt+ii*dt < t_right_ and tt+(ii+1)*dt > t_right_.
  void Print(double dt, ostream * stream = (&cout)) const;
  // Prepends previous output with ",{"name and appends "}"
  void Print(double dt, const string & name, ostream * stream) const;

 private:
  // Stores each interpolation as a vector of Interp objects.
  vector<Interp*> xx_;
  // Size of the interpolated vector.
  size_t xx_pt_size_;
};

#endif // __Agile_RRT__interp__