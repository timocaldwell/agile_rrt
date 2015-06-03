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
private:
    gsl_interp_accel * acc;
    gsl_spline * spline;
    vector<double> xxraw_vec;
    vector<double> yyraw_vec;
    double xx_left, xx_right;
//    vector<double> xx_int;
    
public:
    Interp();
    Interp(const vector<double> &, const vector<double> &);
    Interp(const Interp &);
//    Interp& operator= (const Interp& other)
//    {
//        Interp copy(other); // re-use copy-constructor
//        *this = std::move(copy); // re-use move-assignment
//        return *this;
//    }
    ~Interp();
    void set(const vector<double> &, const vector<double> &);
    void print(double) const;
    double pt(double) const; // Note extrapolation holds to left or right bound.
    double d_pt(double) const; //Derivative at pt
    double dd_pt(double) const; // 2nd Derivative at pt
//    double integ(double, double);
    double begin() const {return xx_left;};
    double end() const {return xx_right;};
};


class InterpVector {
private:
    vector< Interp* > xx;
    size_t xx_pt_size;
public:
    InterpVector();
    InterpVector(const vector<double> &, const vector<vector<double>> &);
    InterpVector(const InterpVector &);
//    InterpVector(InterpVector *);
//    InterpVector& operator= (const InterpVector& other)
//    {
//        InterpVector copy(other); // re-use copy-constructor
//        *this = std::move(copy); // re-use move-assignment
//        return *this;
//    }
    ~InterpVector();
    void set(const vector<double> &, const vector<vector<double>> &);
    void print(double, ostream & stream = cout) const;
    void print(double, const string &, ostream & stream) const;
    void pt(double, vector<double> &) const;
};

#endif
