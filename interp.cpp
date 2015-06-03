//
//  interp.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 9/20/14.
//  Copyright (c) 2014 Timothy Caldwell. All rights reserved.
//

//#include <stdio.h>
#include "interp.h"

Interp::Interp(){
//    acc = nullptr;
//    spline = nullptr;
    xx_left = 0;
    xx_right = 0;
}


Interp::Interp(const vector<double> & xx_vec, const vector<double> & yy_vec){
    set(xx_vec, yy_vec);
}

Interp::Interp(const Interp & other){ //copy constructor
    set(other.xxraw_vec, other.yyraw_vec);
}

Interp::~Interp(){
    
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

//Interp::Interp& operator= (const Interp& other)
//{
//    Interp copy(other); // re-use copy-constructor
//    *this = std::move(copy); // re-use move-assignment
//    return *this;
//}


void Interp::set(const vector<double> & xx_vec, const vector<double> & yy_vec){
//    xx_int.resize(0);
    xxraw_vec = xx_vec;
    yyraw_vec = yy_vec;

    if ( xx_vec.size() <= 1) {
        cout << " ERROR --- Interp input vectors' size too small" << endl;
        abort();
    }

    xx_left = xx_vec[0];
    xx_right = xx_vec.back();

    if (  xxraw_vec.size() == 2) { // inserts a 3rd point through linear interpolation
        xxraw_vec.insert(xxraw_vec.begin()+1, (xx_right + xx_left) * 0.5);
        yyraw_vec.insert(yyraw_vec.begin()+1, (yyraw_vec[0] + yyraw_vec.back()) * 0.5);
    }

    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, xxraw_vec.size());
    
    gsl_spline_init(spline, &xxraw_vec[0], &yyraw_vec[0], xxraw_vec.size());
}

double Interp::pt(double xx) const {
// extrapolation holds to left or right bound.
    if ( xx > xx_right )
        return gsl_spline_eval(spline, xx_right, acc);
    else if ( xx < xx_left )
        return gsl_spline_eval(spline, xx_left, acc);
    return gsl_spline_eval(spline, xx, acc);
}

double Interp::d_pt(double xx) const {
    return gsl_spline_eval_deriv(spline,xx,acc);
}

double Interp::dd_pt(double xx) const {
    return gsl_spline_eval_deriv2(spline,xx,acc);
}

void Interp::print(double dt) const {
    double xx = xx_left;
    cout << "{";
    for (int ii = 0; ii < 1000; ii++ ){
        cout << std::setprecision(6) << "{" << xx << "," << pt(xx) << "}";
        xx += dt;
        if (xx <= xx_right)
            cout << ",";
        else {
            cout << "}" << endl;
            break;
        }
    }
}


InterpVector::InterpVector(){
    xx_pt_size = 0;
}

InterpVector::InterpVector(const vector<double> & tt_vec, const vector<vector<double>> & xx_vec){
    set(tt_vec, xx_vec);
}

InterpVector::InterpVector(const InterpVector& other) { // copy constructor
    xx_pt_size = other.xx_pt_size;
    xx.clear();
    
    Interp * xx_copy;
    for (vector< Interp* >::const_iterator it = other.xx.begin(); it!=other.xx.end(); ++it) {
        xx_copy = new Interp(**it);
        xx.push_back(xx_copy);
    }
}

//InterpVector::InterpVector(InterpVector * other) { // copy constructor
//    xx_pt_size = other->xx_pt_size;
//    xx.clear();
//    Interp * xx_copy;
//    for (vector< Interp* >::const_iterator it = other->xx.begin(); it!=other->xx.end(); ++it) {
//        xx_copy = new Interp(**it);
//        xx.push_back(xx_copy);
//    }
//}


void InterpVector::set(const vector<double> & tt_vec, const vector<vector<double>> & xx_vec){
    Interp * interp_temp;
    xx_pt_size = xx_vec[0].size();
    for (size_t ii = 0; ii < xx_pt_size; ++ii)
    {
        vector<double> xx_temp;
        for (vector<vector<double>>::const_iterator it = xx_vec.begin(); it!=xx_vec.end(); ++it)
            xx_temp.push_back((*it)[ii]);
        interp_temp = new Interp(tt_vec , xx_temp);
        xx.push_back(interp_temp);
    }
}

InterpVector::~InterpVector(){
    if(xx.size()!=0){
        for(unsigned ii = 0; ii<xx.size(); ++ii)
            delete xx.at(ii);
        xx.clear();
    }
}

void InterpVector::pt(double tt, vector<double> & xx_pt) const {
    xx_pt.clear();
    for(unsigned ii = 0; ii<xx_pt_size; ++ii)
        xx_pt.push_back(xx[ii]->pt(tt));
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

void InterpVector::print(double dt, ostream & stream) const {
    double tt;
    if(xx_pt_size == 0) {
        stream << ",{}";
        return;
    }
    stream << ",{";
    for (unsigned jj = 0; jj < xx_pt_size; ++jj){
        tt = xx[0]->begin();
        stream << "{";
        for (int ii = 0; ii < 1000; ++ii ){
            if ( xx[jj]->pt(tt) > .0001 || xx[jj]->pt(tt) < -.0001) // CHOP FOR PRINTTING
                stream << "{" << tt << "," << std::setprecision(6) <<  xx[jj]->pt(tt) << "}";
            else
                stream << "{" << tt << "," << std::setprecision(6) <<  0 << "}";
            tt += dt;
            if (tt <= xx[0]->end())
                stream << ",";
            else {
                stream << "}";
                break;
            }
        }
        if (jj < xx_pt_size - 1)
            stream << ",";
        else
            stream << "}" << endl;
    }
}

void InterpVector::print(double dt, const string & name, ostream & stream) const {
    stream << ",{" << name ;
    print(dt, stream);
    stream << "}";
}
