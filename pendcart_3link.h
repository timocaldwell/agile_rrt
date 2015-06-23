//                                                                                                                            
//  pendcart_3link.h                                                                                                          
//  Agile_RRT                                                                                                                 
//                                                                                                                            
//  Created by Timothy M. Caldwell.                                                                                           
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt                                                              
//                                                                                                                            
//  Within is class PendCart which provides the system equations for a triple pendulum on a cart dynamic system of the form dxxdtt = f(xx, uu).  PendCart(...) constructs the object to be integrated by odeint (see system.h for integrations) through overloaded operator(). Different constructors correspond to differing control feedback schemes including those used for the Agile RRT. The methods CalcAA and CalcBB return the linearization about the state and control respectively for the dynamic system.
//                                                                                                                            
//  triple pendulum on a cart physical characteristics:                                                                       
//  totallen: 1.0 m;                                                                                                          
//  lenLink: 0.333333 m;                                                                                                      
//  gravity: 9.81 m/s/s;                                                                                                      
//  head mass: 0.1 Kg;                                                                                                        
//  body mass: 1.0 Kg;                                                                                                        
//

#ifndef __Agile_RRT__pendcart__
#define __Agile_RRT__pendcart__

#include "global.h"
#include "interp.h"

class PendCart {
 public:
  // For constructors:                                                                                                                                    
  //  1. The arguments to the constructor determine the control feedback form for integration by odeint. A point in control or state space is a vector<double>. A curve in control or state space is an InterpVector. All InterpVectors are interpolated over time tt (see interp.h).
  //  2. If constraints are passed, a boolean checking constraint satisfaction can be checked by calling IsConstraintSatisfied().                         
  //  3. If pointer uu_out is passed, the current point uu is stored at location uu_out. Note uu_out is needed for checking constraints.

  // free dynamics---i.e. zero input control---uu = 0
  PendCart();
  // feedforward control---dxxdtt = uu = uu_ff_interp
  PendCart(const InterpVector & uu_ff_interp);
  // track point---uu = kk_interp*(xx_ref_pt - xx)
  PendCart(const InterpVector & KK_interp, const ode_state_type & xx_ref_pt);
  // track trajectory with feedforward---uu = uu_ff_interp + kk_interp*(xx_ref_interp - xx)
  PendCart(const InterpVector & uu_ff_interp, const InterpVector & KK_interp, const InterpVector & xx_ref_interp);
  PendCart(const InterpVector & uu_ff_interp, const InterpVector & KK_interp, const InterpVector & xx_ref_interp, const constraints_struct & constraints, ode_state_type * uu_out);
  // For Agile RRT implemenation: Inexact linear steering with projection                                                                                  
  //   xxtilde = xxzero - WWK_mat*Phi^T*eta                                                                                                               
  //   uutilde = (KKlin*WWK - RRinv*BB^T)*Phi^T*eta                                                                                                       
  //   uu = uutilde + KKproj*(xxtilde - xx)
  PendCart(const InterpVector & xxzero, const InterpVector & BB, const InterpVector & KKlin, const InterpVector & KKproj, const InterpVector & WWK, const InterpVector & Phi, const NSx1_type & eta, const NIxNI_type & RRinv, const constraints_struct & constraints, ode_state_type * uu_out);
  // For Agile RRT implemenation: Inexact inear steering with projection (same as previous).                                                              
  //   Additionally computes the running cost ell = 0.5*xx^T*QQ*xx + 0.5*uu^T*RR*uu at position NS in dxxdtt---i.e. f(xx,uu) is length NS+1 now.
  PendCart(const InterpVector & xxzero, const InterpVector & BB, const InterpVector & KKlin, const InterpVector & KKproj, const InterpVector & WWK, const InterpVector & Phi, const NSx1_type & eta, const NSxNS_type & QQ, const NIxNI_type & RR, const NIxNI_type & RRinv, const constraints_struct & constraints, ode_state_type * uu_out);

  // Computes dxxdtt = f(xx,uu) for time tt. For use by odeint through method try_step (refer to system.h).
  void operator()( const ode_state_type xx , ode_state_type &dxxdtt , const double tt );

  // Passes the linearization of f(xx,uu) with respect to xx for points xx and uu to AA
  template<typename TT>
  void CalcAA(const vector<double> & xx, const vector<double> & uu, Eigen::MatrixBase<TT> * AA);
  
  // Passes the linearization of f(xx,uu) with respect to uu for points xx and uu to BB
  template<typename TT>
  void CalcBB(const vector<double> & xx, const vector<double> &uu, Eigen::MatrixBase<TT> * BB);

  // Returns false if linear interpoloation between xx_prev and xx_cur violates kinematic constraints specified by member variable constraints.           
  // Returns false if uu_cur is outside the bounds specified by member variable constraints.
  bool IsConstraintSatisfied(const vector<double> &, const vector<double> &, const vector<double> &) const;

 private:
  const double linklen_ = 0.333333;
  const int numlinks_ = 3;
  vector<double> uu_;
  int type_;
  const InterpVector * uu_ff_interp_;         // Feedforward trajectory.                                                                                  
  const InterpVector * KK_interp_;            // Feedback Gains.                                                                                          
  const InterpVector * xx_ref_interp_;        // Reference for trajectory tracking.                                                                       
  const ode_state_type * xx_ref_pt_;          // Feedback point.                                                                                          

  // for efficient inexact linear steering and projection for Agile RRT implementation (refer to RRT paper)                                               
  const InterpVector * xxzero_;               // Zero-control trajectory---i.e. xx when dxxdtt = f(xx,0).                                                 
  const InterpVector * BB_;                   // Linearization about uu.                                                                                  
  const InterpVector * KKlin_;                // Gains for linear steering.                                                                               
  const InterpVector * KKproj_;               // Trajectory tracking projection gains.                                                                    
  const InterpVector * WWK_;                  // Reachability Gramian.                                                                                    
  const InterpVector * Phi_;                  // State-transition matrix for a linearization AA.                                                          
  const NSx1_type * eta_;                     // Vector that determines linear steering location.                                                         

  // Weighting matrices QQ = QQ^T>=0, RR = RR^T>0, RRinv = RR^{-1}                                                                                        
  const NSxNS_type * QQ_;
  const NIxNI_type * RR_;
  const NIxNI_type * RRinv_;

  // Constraints on state and control to be checked through IsConstraintSatisfied(). When constraints_==nullptr, IsConstraintSatisfied() returns true.    
  const constraints_struct * constraints_;
  
  // The location of the control input used in integration. This pointer is needed since odeint does not provide a direct way to retain details of the integration.
  ode_state_type * uu_out_;

  // The dynamics for the system---dxxdtt = ff(xx,uu)
  void ff(const ode_state_type & xx, const vector<double> & uu, ode_state_type * dxxdtt);
  
  // Returns true if line segment between points (ppX, ppY) and (qqX, qqY) intersects circle at point (rrX, rrY) with squared radius radius_squared.
  bool IsLineCircleIntersect(double radius_squared, double ppX, double ppY, double qqX, double qqY, double rrX, double rrY) const;
};

// Hard coded linearization of ff(xx, uu) about xx for points xx and uu. The equations were exported from Mathematica.
template<typename TT>
void PendCart::CalcAA(const vector<double> & xx, const vector<double> & uu, Eigen::MatrixBase<TT> * AA) {
  EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(TT, 8, 8);

  double cx0 = cos(xx[0]);
  double sx0 = sin(xx[0]);
  double cx2 = cos(xx[2]);
  double sx2 = sin(xx[2]);
  double cx4 = cos(xx[4]);
  double sx4 = sin(xx[4]);
  
  double cx0x2 = cos(xx[0] - xx[2]);
  double sx0x2 = sin(xx[0] - xx[2]);
  double cx0x4 = cos(xx[0] - xx[4]);
  double sx0x4 = sin(xx[0] - xx[4]);
  double cx2x4 = cos(xx[2] - xx[4]);
  double sx2x4 = sin(xx[2] - xx[4]);
  
  double sqx1 = pow(xx[1],2);
  double sqx3 = pow(xx[3],2);
  double sqx5 = pow(xx[5],2);
  
  double sqcx4 = pow(cx4,2);
  
  double sqcx0x2 = pow(cx0x2,2);
  double sqcx0x4 = pow(cx0x4,2);
  double sqcx2x4 = pow(cx2x4,2);
  double sqsx2x4 = pow(sx2x4,2);
  
  double sqvar0 = pow(cx0x2/36450. - (cx0x4*cx2x4)/72900.,2);
  double sqvar1 = pow(-sqvar0 +
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2);
  double sqvar2 = pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2);
  double sqvar3 = pow(cx2/12150. - (cx2x4*cx4)/24300.,2);
  double sqvar5 = pow(-0.06666666666666667 + sqcx2x4/30.,2);
  double sqvar6 = pow((0.000027434842249657064 - sqcx2x4/72900.)*
          (cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2);
  double sqvar7 = pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
       (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2);
  double sqvar8 = pow(-sqvar7 +
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)),2);
  
  double sqvar9 = pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
        (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2);
  
  (*AA)( 0, 0 ) = 0;
  (*AA)( 0, 1 ) = 1;
  (*AA)( 0, 2 ) = 0;
  (*AA)( 0, 3 ) = 0;
  (*AA)( 0, 4 ) = 0;
  (*AA)( 0, 5 ) = 0;
  (*AA)( 0, 6 ) = 0;
  (*AA)( 0, 7 ) = 0;

  (*AA)( 1, 0 ) = (((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. -
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450 + (cx2x4*sx0x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. -
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  pow(-sqvar0 +
    (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2) - 
 (-((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)) - 
    (-sx0x2/36450. + (cx2x4*sx0x4)/72900.)*
     ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
       (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
    (0.000027434842249657064 - sqcx2x4/72900.)*
     ((sqx1*sqcx0x4)/72900. + ((-327*cx0)/1000. + (sqx3*cx0x2)/135. + (sqx1*cx0x4)/270.)/270. + 
       (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/
  (-sqvar0 + 
    (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.)) - 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  (sqvar1*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) + 
 ((-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
      (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) - 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    pow(-sqvar2 +
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)),2)) + 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    ((-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*(-(sqx1*cx0)/8100. + (sqx1*cx0x4*cx4)/24300.) - 
         (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)) - 
      (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*
          (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)) - 
         (-sx0x2/36450. + (cx2x4*sx0x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          ((sqx1*sqcx0x4)/72900. + ((-327*cx0)/1000. + (sqx3*cx0x2)/135. + (sqx1*cx0x4)/270.)/
             270. + (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
      ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 1, 1 ) = -(((xx[1]*cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.))/
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((xx[1]*cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(xx[1]*sx0)/4050. + (xx[1]*cx4*sx0x4)/12150.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 1, 2 ) = ((((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  pow(-sqvar0 +
    (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2) - 
 ((-(sqx3*cx0x2)/36450. + (sqx3*cx0x4*cx2x4)/72900.)*(0.000027434842249657064 - sqcx2x4/72900.) + 
    (cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
    (sx0x2/36450. + (cx0x4*sx2x4)/72900.)*
     ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
       (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
    (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
     ((sqx3*sqcx2x4)/72900. + ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx1*cx2x4)/270.)/270. + 
       (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/
  (-sqvar0 + 
    (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.)) - 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  (sqvar1*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) + 
 (((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
      (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) - 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       ((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((cx2x4*(0.004074074074074074 - sqcx4/8100.)*sx2x4)/36450. - 
         2*(cx2/12150. - (cx2x4*cx4)/24300.)*(-sx2/12150. + (cx4*sx2x4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    pow(-sqvar2 +
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)),2)) + 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    (-(((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
           (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       ((-(sqx3*cx0x2)/36450. + (sqx3*cx0x4*cx2x4)/72900.)*
          (0.000027434842249657064 - sqcx2x4/72900.) + 
         (cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
         (sx0x2/36450. + (cx0x4*sx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((sqx3*sqcx2x4)/72900. + ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx1*cx2x4)/270.)/
             270. + (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
      (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*(-(sqx3*cx2)/12150. + (sqx3*cx2x4*cx4)/24300.) - 
         (-sx2/12150. + (cx4*sx2x4)/24300.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*((sqx3*sqcx2x4)/72900. + 
            ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx1*cx2x4)/270.)/270. + 
            (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
         (cx2x4*sx2x4*(-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
              (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))/36450.)))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 1, 3 ) = -((-(xx[3]*cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450. + 
      (0.000027434842249657064 - sqcx2x4/72900.)*((xx[3]*sx0x2)/18225. + (xx[3]*cx0x4*sx2x4)/36450.))/
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-(xx[3]*cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450. + 
           (0.000027434842249657064 - sqcx2x4/72900.)*((xx[3]*sx0x2)/18225. + (xx[3]*cx0x4*sx2x4)/36450.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(xx[3]*cx2x4*(cx2/12150. - (cx2x4*cx4)/24300.)*sx2x4)/36450. + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(xx[3]*sx2)/6075. + (xx[3]*cx4*sx2x4)/12150.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 1, 4 ) = ((-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  pow(-sqvar0 +
    (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2) - 
 (-(cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
          (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
    (-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)*
     ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
       (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
    (0.000027434842249657064 - sqcx2x4/72900.)*
     (-(sqx1*cx0x4)/72900. - (cx0x4*((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/
        270. - (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
    (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
     (-(sqx1*cx2x4)/72900. - (cx2x4*((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/
        270. - (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/
  (-sqvar0 + 
    (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.)) - 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  (sqvar1*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) + 
 ((-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
      (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
      (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) - 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(cx2x4*(0.004074074074074074 - sqcx4/8100.)*sx2x4)/36450. + 
         ((0.000027434842249657064 - sqcx2x4/72900.)*cx4*sx4)/4050. - 
         2*(cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    pow(-sqvar2 +
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)),2)) + 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-(cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
                 (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
           (-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
           (0.000027434842249657064 - sqcx2x4/72900.)*
            (-(sqx1*cx0x4)/72900. - (cx0x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            (-(sqx1*cx2x4)/72900. - (cx2x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*
          (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))*
       (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)) + 
      (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            (-(sqx1*cx2x4)/72900. - (cx2x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) - 
         ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)*
          (-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(sqx1*cx4)/24300. - (((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.)*cx4)/90. + 
            ((-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.)*sx4)/90.) - 
         (cx2x4*sx2x4*(-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
              (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))/36450.)))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 1, 5 ) = -(((xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      (xx[5]*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450.)/
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
           (xx[5]*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(xx[5]*(cx2/12150. - (cx2x4*cx4)/24300.)*sx2x4)/36450. - 
         (xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx4)/12150.)))/
  ((-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
    (-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 1, 6 ) = 0;
  (*AA)( 1, 7 ) = 0;

  (*AA)( 2, 0 ) = 0;
  (*AA)( 2, 1 ) = 0;
  (*AA)( 2, 2 ) = 0;
  (*AA)( 2, 3 ) = 1;
  (*AA)( 2, 4 ) = 0;
  (*AA)( 2, 5 ) = 0;
  (*AA)( 2, 6 ) = 0;
  (*AA)( 2, 7 ) = 0;
  
  (*AA)( 3, 0 ) = (-3*((sqx1*cx0x2)/45. - (sqx1*cx0x4*cx2x4)/90.))/(-0.06666666666666667 + sqcx2x4/30.) + 
 ((-sx0x2/36450. + (cx2x4*sx0x4)/72900.)*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    sqvar1) + 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    (-((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
         (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)) - 
      (-sx0x2/36450. + (cx2x4*sx0x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
      (0.000027434842249657064 - sqcx2x4/72900.)*
       ((sqx1*sqcx0x4)/72900. + ((-327*cx0)/1000. + (sqx3*cx0x2)/135. + (sqx1*cx0x4)/270.)/
          270. + (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    ((-2*cx0x2*sx0x2)/225. + (cx0x4*cx2x4*sx0x2)/225. - (cx0x4*sx0x4)/225. + 
      (cx0x2*cx2x4*sx0x4)/225.)*(-(((0.000027434842249657064 - sqcx2x4/72900.)*
            (cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  (pow(-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. +
      sqcx2x4/300.,2)*(-sqvar6 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) + 
 (3*(-(cx0x2*sx0)/150. + (cx0x4*cx2x4*sx0)/300. - (cx0*sx0x2)/150. + 
      (cx0x4*cx4*sx0x2)/450. - (cx2*cx0x4*sx0x4)/225. + (cx0*cx2x4*sx0x4)/300. + 
      (cx0x2*cx4*sx0x4)/450.)*(-(((0.000027434842249657064 - sqcx2x4/72900.)*
            (cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) - 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*pow(-sqvar6 +
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)),2)) + 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    ((-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*(-(sqx1*cx0)/8100. + (sqx1*cx0x4*cx4)/24300.) - 
         (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)) - 
      (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*
          (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)) - 
         (-sx0x2/36450. + (cx2x4*sx0x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          ((sqx1*sqcx0x4)/72900. + ((-327*cx0)/1000. + (sqx3*cx0x2)/135. + (sqx1*cx0x4)/270.)/
             270. + (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
      ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 3, 1 ) = (-3*((2*xx[1]*sx0x2)/45. - (xx[1]*cx2x4*sx0x4)/45.))/(-0.06666666666666667 + sqcx2x4/30.) + 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    ((xx[1]*cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((xx[1]*cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(xx[1]*sx0)/4050. + (xx[1]*cx4*sx0x4)/12150.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 3, 2 ) = -(cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4*
     ((0.000027434842249657064 - sqcx2x4/72900.)*
        (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
          (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
       (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
        ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
          (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  (36450.*pow(0.000027434842249657064 - sqcx2x4/72900.,2)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 ((sx0x2/36450. + (cx0x4*sx2x4)/72900.)*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    sqvar1) + 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    ((-(sqx3*cx0x2)/36450. + (sqx3*cx0x4*cx2x4)/72900.)*
       (0.000027434842249657064 - sqcx2x4/72900.) + 
      (cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
           (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
      (sx0x2/36450. + (cx0x4*sx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((sqx3*sqcx2x4)/72900. + ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx1*cx2x4)/270.)/270. + 
         (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 (cx2x4*sx2x4*((sqx1*sx0x2)/45. + (327*sx2)/500. - (sqx1*cx2x4*sx0x4)/90. - 
      (sqx1*sx2x4)/90. - (sqx3*cx2x4*sx2x4)/90. - (327*cx2x4*sx4)/1000.))/
  (5.*sqvar5) - 
 (3*(-(sqx1*cx0x2)/45. + (327*cx2)/500. - (sqx1*cx2x4)/90. - (sqx3*sqcx2x4)/90. + 
      (sqx1*sx0x4*sx2x4)/90. + (sqx3*sqsx2x4)/90. + (327*sx2x4*sx4)/1000.))/
  (-0.06666666666666667 + sqcx2x4/30.) - (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - 
      (cx0*cx0x4*cx2x4)/300. - (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    ((2*cx0x2*sx0x2)/225. - (cx0x4*cx2x4*sx0x2)/225. + (cx0x2*cx0x4*sx2x4)/225. - 
      (cx2x4*sx2x4)/150.)*(-(((0.000027434842249657064 - sqcx2x4/72900.)*
            (cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  (pow(-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. +
      sqcx2x4/300.,2)*(-sqvar6 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) + 
 (3*((cx0*sx0x2)/150. - (cx0x4*cx4*sx0x2)/450. + sx2/150. - (sqcx0x4*sx2)/450. + 
      (cx0*cx0x4*sx2x4)/300. - (cx4*sx2x4)/300.)*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) - 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       ((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((cx2x4*(0.004074074074074074 - sqcx4/8100.)*sx2x4)/36450. - 
         2*(cx2/12150. - (cx2x4*cx4)/24300.)*(-sx2/12150. + (cx4*sx2x4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*pow(-sqvar6 +
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)),2)) + 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    (-(((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
           (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       ((-(sqx3*cx0x2)/36450. + (sqx3*cx0x4*cx2x4)/72900.)*
          (0.000027434842249657064 - sqcx2x4/72900.) + 
         (cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
         (sx0x2/36450. + (cx0x4*sx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((sqx3*sqcx2x4)/72900. + ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx1*cx2x4)/270.)/
             270. + (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
      (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*(-(sqx3*cx2)/12150. + (sqx3*cx2x4*cx4)/24300.) - 
         (-sx2/12150. + (cx4*sx2x4)/24300.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*((sqx3*sqcx2x4)/72900. + 
            ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx1*cx2x4)/270.)/270. + 
            (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
         (cx2x4*sx2x4*(-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
              (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))/36450.)))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 3, 3 ) = (xx[3]*cx2x4*sx2x4)/(15.*(-0.06666666666666667 + sqcx2x4/30.)) + 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    (-(xx[3]*cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450. + 
      (0.000027434842249657064 - sqcx2x4/72900.)*((xx[3]*sx0x2)/18225. + (xx[3]*cx0x4*sx2x4)/36450.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-(xx[3]*cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450. + 
           (0.000027434842249657064 - sqcx2x4/72900.)*((xx[3]*sx0x2)/18225. + (xx[3]*cx0x4*sx2x4)/36450.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(xx[3]*cx2x4*(cx2/12150. - (cx2x4*cx4)/24300.)*sx2x4)/36450. + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(xx[3]*sx2)/6075. + (xx[3]*cx4*sx2x4)/12150.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 3, 4 ) = (cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  (36450.*pow(0.000027434842249657064 - sqcx2x4/72900.,2)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 ((-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    sqvar1) + 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    (-(cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
      (-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
      (0.000027434842249657064 - sqcx2x4/72900.)*
       (-(sqx1*cx0x4)/72900. - (cx0x4*
            ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
         (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       (-(sqx1*cx2x4)/72900. - (cx2x4*
            ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
         (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (cx2x4*sx2x4*((sqx1*sx0x2)/45. + (327*sx2)/500. - (sqx1*cx2x4*sx0x4)/90. - 
      (sqx1*sx2x4)/90. - (sqx3*cx2x4*sx2x4)/90. - (327*cx2x4*sx4)/1000.))/
  (5.*sqvar5) - 
 (3*((sqx1*cx2x4)/90. + (sqx1*cx0x4*cx2x4)/90. + (sqx3*sqcx2x4)/90. - 
      (327*cx2x4*cx4)/1000. - (sqx1*sx0x4*sx2x4)/90. - (sqx3*sqsx2x4)/90. - 
      (327*sx2x4*sx4)/1000.))/(-0.06666666666666667 + sqcx2x4/30.) - 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    ((cx0x4*sx0x4)/225. - (cx0x2*cx2x4*sx0x4)/225. - (cx0x2*cx0x4*sx2x4)/225. + 
      (cx2x4*sx2x4)/150.)*(-(((0.000027434842249657064 - sqcx2x4/72900.)*
            (cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  (pow(-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. +
      sqcx2x4/300.,2)*(-sqvar6 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) + 
 (3*((cx2*cx0x4*sx0x4)/225. - (cx0*cx2x4*sx0x4)/300. - (cx0x2*cx4*sx0x4)/450. - 
      (cx0*cx0x4*sx2x4)/300. + (cx4*sx2x4)/300. + (cx0x2*cx0x4*sx4)/450. - (cx2x4*sx4)/300.
      )*(-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)))) - 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(cx2x4*(0.004074074074074074 - sqcx4/8100.)*sx2x4)/36450. + 
         ((0.000027434842249657064 - sqcx2x4/72900.)*cx4*sx4)/4050. - 
         2*(cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*pow(-sqvar6 +
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)),2)) + 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-(cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
                 (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
           (-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
           (0.000027434842249657064 - sqcx2x4/72900.)*
            (-(sqx1*cx0x4)/72900. - (cx0x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            (-(sqx1*cx2x4)/72900. - (cx2x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*
          (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx1*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))*
       (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)) + 
      (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            (-(sqx1*cx2x4)/72900. - (cx2x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) - 
         ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx1*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)*
          (-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(sqx1*cx4)/24300. - (((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.)*cx4)/90. + 
            ((-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.)*sx4)/90.) - 
         (cx2x4*sx2x4*(-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
              (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx1*sx4)/90.)/270.))/36450.)))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 3, 5 ) = (xx[5]*sx2x4)/(15.*(-0.06666666666666667 + sqcx2x4/30.)) + 
 ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
    ((xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      (xx[5]*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450.))/
  ((0.000027434842249657064 - sqcx2x4/72900.)*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
      (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
           (xx[5]*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(xx[5]*(cx2/12150. - (cx2x4*cx4)/24300.)*sx2x4)/36450. - 
         (xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx4)/12150.)))/
  ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
      sqcx2x4/300.)*(-sqvar2 + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 3, 6 ) = 0;
  (*AA)( 3, 7 ) = 0;
  (*AA)( 4, 0 ) = 0;
  (*AA)( 4, 1 ) = 0;
  (*AA)( 4, 2 ) = 0;
  (*AA)( 4, 3 ) = 0;
  (*AA)( 4, 4 ) = 0;
  (*AA)( 4, 5 ) = 1;
  (*AA)( 4, 6 ) = 0;
  (*AA)( 4, 7 ) = 0;
  
  (*AA)( 5, 0 ) = sqx1*cx0x4 + (3*cx2x4*((sqx1*cx0x2)/45. - (sqx1*cx0x4*cx2x4)/90.))/
  (-0.06666666666666667 + sqcx2x4/30.) + (3*(-(cx2x4*sx0x2)/45. + sx0x4/45.)*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/
       36450. - 2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*pow(-sqvar0 +
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2)) + 
 (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*(-((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
         (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)) - 
      (-sx0x2/36450. + (cx2x4*sx0x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
      (0.000027434842249657064 - sqcx2x4/72900.)*
       ((sqx1*sqcx0x4)/72900. + ((-327*cx0)/1000. + (sqx3*cx0x2)/135. + (sqx5*cx0x4)/270.)/
          270. + (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (((-3*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-(cx2x4*sx0x2)/45. + sx0x4/45.))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
           2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         pow(-sqvar0 +
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2)) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
           (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))) - 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  sqvar8 + 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    ((-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*(-(sqx1*cx0)/8100. + (sqx1*cx0x4*cx4)/24300.) - 
         (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)) - 
      (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*
          (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)) - 
         (-sx0x2/36450. + (cx2x4*sx0x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          ((sqx1*sqcx0x4)/72900. + ((-327*cx0)/1000. + (sqx3*cx0x2)/135. + (sqx5*cx0x4)/270.)/
             270. + (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
      ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
  (*AA)( 5, 1 ) = 2*xx[1]*sx0x4 + (3*cx2x4*((2*xx[1]*sx0x2)/45. - (xx[1]*cx2x4*sx0x4)/45.))/
  (-0.06666666666666667 + sqcx2x4/30.) + (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
    ((xx[1]*cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((xx[1]*cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(xx[1]*sx0)/4050. + (xx[1]*cx4*sx0x4)/12150.))))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
  (*AA)( 5, 2 ) = sqx3*cx2x4 + (cx2x4*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*sx2x4*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  (5.*sqvar5*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (3*((cx2x4*sx0x2)/45. - (cx0x2*sx2x4)/45.)*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*(((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/
       36450. - 2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*pow(-sqvar0 +
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2)) + 
 (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*((-(sqx3*cx0x2)/36450. + (sqx3*cx0x4*cx2x4)/72900.)*
       (0.000027434842249657064 - sqcx2x4/72900.) + 
      (cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
           (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
      (sx0x2/36450. + (cx0x4*sx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((sqx3*sqcx2x4)/72900. + ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx5*cx2x4)/270.)/270. + 
         (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (sqcx2x4*sx2x4*((sqx1*sx0x2)/45. + (327*sx2)/500. - (sqx1*cx2x4*sx0x4)/90. - 
      (sqx5*sx2x4)/90. - (sqx3*cx2x4*sx2x4)/90. - (327*cx2x4*sx4)/1000.))/
  (5.*sqvar5) - 
 (3*sx2x4*((sqx1*sx0x2)/45. + (327*sx2)/500. - (sqx1*cx2x4*sx0x4)/90. - 
      (sqx5*sx2x4)/90. - (sqx3*cx2x4*sx2x4)/90. - (327*cx2x4*sx4)/1000.))/
  (-0.06666666666666667 + sqcx2x4/30.) + (3*cx2x4*
    (-(sqx1*cx0x2)/45. + (327*cx2)/500. - (sqx5*cx2x4)/90. - (sqx3*sqcx2x4)/90. + 
      (sqx1*sx0x4*sx2x4)/90. + (sqx3*sqsx2x4)/90. + (327*sx2x4*sx4)/1000.))/
  (-0.06666666666666667 + sqcx2x4/30.) + (((cx2x4*((cx2*cx2x4)/15. - cx4/15.)*sx2x4)/
       (5.*sqvar5) - 
      (cx2x4*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*sx2x4)/
       (5.*sqvar5*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
      (3*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((cx2x4*sx0x2)/45. - (cx0x2*sx2x4)/45.))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
      (3*(-(cx2x4*sx2)/15. - (cx2*sx2x4)/15.))/(-0.06666666666666667 + sqcx2x4/30.) + 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
           2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         pow(-sqvar0 +
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2)) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
           (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))) - 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       ((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((cx2x4*(0.004074074074074074 - sqcx4/8100.)*sx2x4)/36450. - 
         2*(cx2/12150. - (cx2x4*cx4)/24300.)*(-sx2/12150. + (cx4*sx2x4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  sqvar8 + 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    (-(((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
           (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       ((-(sqx3*cx0x2)/36450. + (sqx3*cx0x4*cx2x4)/72900.)*
          (0.000027434842249657064 - sqcx2x4/72900.) + 
         (cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
         (sx0x2/36450. + (cx0x4*sx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((sqx3*sqcx2x4)/72900. + ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx5*cx2x4)/270.)/
             270. + (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
      (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*(-(sqx3*cx2)/12150. + (sqx3*cx2x4*cx4)/24300.) - 
         (-sx2/12150. + (cx4*sx2x4)/24300.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*((sqx3*sqcx2x4)/72900. + 
            ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx5*cx2x4)/270.)/270. + 
            (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
         (cx2x4*sx2x4*(-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
              (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))/36450.)))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
  (*AA)( 5, 3 ) = 2*xx[3]*sx2x4 - (xx[3]*sqcx2x4*sx2x4)/(15.*(-0.06666666666666667 + sqcx2x4/30.)) + 
 (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*(-(xx[3]*cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          sx2x4)/36450. + (0.000027434842249657064 - sqcx2x4/72900.)*
       ((xx[3]*sx0x2)/18225. + (xx[3]*cx0x4*sx2x4)/36450.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-(xx[3]*cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450. + 
           (0.000027434842249657064 - sqcx2x4/72900.)*((xx[3]*sx0x2)/18225. + (xx[3]*cx0x4*sx2x4)/36450.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(xx[3]*cx2x4*(cx2/12150. - (cx2x4*cx4)/24300.)*sx2x4)/36450. + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(xx[3]*sx2)/6075. + (xx[3]*cx4*sx2x4)/12150.))))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 +
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
  (*AA)( 5, 4 ) = -3*((sqx1*cx0x4)/3. + (sqx3*cx2x4)/3. - (981*cx4)/100.) - 
 (cx2x4*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*sx2x4*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  (5.*sqvar5*
    (-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (3*(-sx0x4/45. + (cx0x2*sx2x4)/45.)*((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*(-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/
       36450. - ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
      2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.))*
    ((0.000027434842249657064 - sqcx2x4/72900.)*
       (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
         (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*pow(-sqvar0 +
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2)) + 
 (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*(-(cx2x4*sx2x4*
          (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
      (-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)*
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
         (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
      (0.000027434842249657064 - sqcx2x4/72900.)*
       (-(sqx5*cx0x4)/72900. - (cx0x4*
            ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
         (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
      (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
       (-(sqx5*cx2x4)/72900. - (cx2x4*
            ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
         (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
 (sqcx2x4*sx2x4*((sqx1*sx0x2)/45. + (327*sx2)/500. - (sqx1*cx2x4*sx0x4)/90. - 
      (sqx5*sx2x4)/90. - (sqx3*cx2x4*sx2x4)/90. - (327*cx2x4*sx4)/1000.))/
  (5.*sqvar5) + 
 (3*sx2x4*((sqx1*sx0x2)/45. + (327*sx2)/500. - (sqx1*cx2x4*sx0x4)/90. - 
      (sqx5*sx2x4)/90. - (sqx3*cx2x4*sx2x4)/90. - (327*cx2x4*sx4)/1000.))/
  (-0.06666666666666667 + sqcx2x4/30.) + (3*cx2x4*
    ((sqx5*cx2x4)/90. + (sqx1*cx0x4*cx2x4)/90. + (sqx3*sqcx2x4)/90. - 
      (327*cx2x4*cx4)/1000. - (sqx1*sx0x4*sx2x4)/90. - (sqx3*sqsx2x4)/90. - 
      (327*sx2x4*sx4)/1000.))/(-0.06666666666666667 + sqcx2x4/30.) + 
 ((-(cx2x4*((cx2*cx2x4)/15. - cx4/15.)*sx2x4)/(5.*sqvar5) + 
      (cx2x4*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*sx2x4)/
       (5.*sqvar5*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
      (3*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-sx0x4/45. + (cx0x2*sx2x4)/45.))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
           ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
           2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         pow(-sqvar0 +
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.),2)) + 
      (3*((cx2*sx2x4)/15. + sx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
           (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
           (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))) - 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    ((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(cx2x4*(0.004074074074074074 - sqcx4/8100.)*sx2x4)/36450. + 
         ((0.000027434842249657064 - sqcx2x4/72900.)*cx4*sx4)/4050. - 
         2*(cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  sqvar8 + 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         (-(cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
                 (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
           (-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
           (0.000027434842249657064 - sqcx2x4/72900.)*
            (-(sqx5*cx0x4)/72900. - (cx0x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            (-(sqx5*cx2x4)/72900. - (cx2x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) - 
      ((0.000027434842249657064 - sqcx2x4/72900.)*
          (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))*
       (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)) + 
      (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            (-(sqx5*cx2x4)/72900. - (cx2x4*
                 ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
              (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) - 
         ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)*
          (-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(sqx5*cx4)/24300. - (((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.)*cx4)/90. + 
            ((-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.)*sx4)/90.) - 
         (cx2x4*sx2x4*(-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
              (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))/36450.)))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
  (*AA)( 5, 5 ) = -(xx[5]*cx2x4*sx2x4)/(15.*(-0.06666666666666667 + sqcx2x4/30.)) + 
 (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*((xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
      (xx[5]*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450.))/
  ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
 (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
      (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
         ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
       ((-0.06666666666666667 + sqcx2x4/30.)*
         (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
           (xx[5]*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(xx[5]*(cx2/12150. - (cx2x4*cx4)/24300.)*sx2x4)/36450. - 
         (xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx4)/12150.)))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));

  (*AA)( 5, 6 ) = 0;
  (*AA)( 5, 7 ) = 0;
  (*AA)( 6, 0 ) = 0;
  (*AA)( 6, 1 ) = 0;
  (*AA)( 6, 2 ) = 0;
  (*AA)( 6, 3 ) = 0;
  (*AA)( 6, 4 ) = 0;
  (*AA)( 6, 5 ) = 0;
  (*AA)( 6, 6 ) = 0;
  (*AA)( 6, 7 ) = 1;
  
  (*AA)( 7, 0 ) = (((-sqvar3 +
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  sqvar8 - 
 ((-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     ((0.000027434842249657064 - sqcx2x4/72900.)*(-(sqx1*cx0)/8100. + (sqx1*cx0x4*cx4)/24300.) - 
       (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)) - 
    (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.)) + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(-sx0/8100. + (cx4*sx0x4)/24300.))*
     ((0.000027434842249657064 - sqcx2x4/72900.)*
        (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
          (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
       (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
        ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
          (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) - 
    ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
       (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
     (-((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          (-(sqx1*cx0x2)/36450. + (sqx1*cx0x4*cx2x4)/72900.)) - 
       (-sx0x2/36450. + (cx2x4*sx0x4)/72900.)*
        ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
          (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
       (0.000027434842249657064 - sqcx2x4/72900.)*
        ((sqx1*sqcx0x4)/72900. + ((-327*cx0)/1000. + (sqx3*cx0x2)/135. + (sqx5*cx0x4)/270.)/
           270. + (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
    ((cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
       2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx0x2/36450. + (cx2x4*sx0x4)/72900.))*
     (-((cx2/12150. - (cx2x4*cx4)/24300.)*((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/
             270. - (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
       (0.000027434842249657064 - sqcx2x4/72900.)*
        (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
          (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.)))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
  (*AA)( 7, 1 ) = -((-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
          (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
        ((xx[1]*cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
          (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.))) + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-(xx[1]*sx0x2)/18225. + (xx[1]*cx2x4*sx0x4)/36450.)) + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(-(xx[1]*sx0)/4050. + (xx[1]*cx4*sx0x4)/12150.)))/
   (-sqvar9 + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-sqvar3 + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 7, 2 ) = (((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       ((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       ((cx2x4*(0.004074074074074074 - sqcx4/8100.)*sx2x4)/36450. - 
         2*(cx2/12150. - (cx2x4*cx4)/24300.)*(-sx2/12150. + (cx4*sx2x4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  sqvar8 - 
 (-(((cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-sx2/12150. + (cx4*sx2x4)/24300.))*
       ((0.000027434842249657064 - sqcx2x4/72900.)*
          (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) - 
    ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
       (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
     ((-(sqx3*cx0x2)/36450. + (sqx3*cx0x4*cx2x4)/72900.)*
        (0.000027434842249657064 - sqcx2x4/72900.) + 
       (cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
            (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
       (sx0x2/36450. + (cx0x4*sx2x4)/72900.)*
        ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
          (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
       (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
        ((sqx3*sqcx2x4)/72900. + ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx5*cx2x4)/270.)/
           270. + (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
    (((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
       2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(sx0x2/36450. + (cx0x4*sx2x4)/72900.))*
     (-((cx2/12150. - (cx2x4*cx4)/24300.)*((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/
             270. - (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
       (0.000027434842249657064 - sqcx2x4/72900.)*
        (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
          (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.)) + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     ((0.000027434842249657064 - sqcx2x4/72900.)*(-(sqx3*cx2)/12150. + (sqx3*cx2x4*cx4)/24300.) - 
       (-sx2/12150. + (cx4*sx2x4)/24300.)*((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/
           270. - (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
       (cx2/12150. - (cx2x4*cx4)/24300.)*((sqx3*sqcx2x4)/72900. + 
          ((sqx1*cx0x2)/135. - (109*cx2)/500. + (sqx5*cx2x4)/270.)/270. + 
          (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
       (cx2x4*sx2x4*(-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))/36450.))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
  (*AA)( 7, 3 ) = -((-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
          (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
        (-(xx[3]*cx2x4*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450. + 
          (0.000027434842249657064 - sqcx2x4/72900.)*((xx[3]*sx0x2)/18225. + (xx[3]*cx0x4*sx2x4)/36450.))) + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-(xx[3]*cx2x4*(cx2/12150. - (cx2x4*cx4)/24300.)*sx2x4)/36450. + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(-(xx[3]*sx2)/6075. + (xx[3]*cx4*sx2x4)/12150.)))/
   (-sqvar9 + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-sqvar3 + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 7, 4 ) = (((-sqvar3 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))*
       (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
         ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
         2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)) - 
      2*((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
         (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-(cx2x4*(0.004074074074074074 - sqcx4/8100.)*sx2x4)/36450. + 
         ((0.000027434842249657064 - sqcx2x4/72900.)*cx4*sx4)/4050. - 
         2*(cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)))*
    (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
         ((0.000027434842249657064 - sqcx2x4/72900.)*
            (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
              (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-((cx2/12150. - (cx2x4*cx4)/24300.)*
            ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
              (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))))/
  sqvar8 - 
 (-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
       (-(cx2x4*sx2x4*(((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
               (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/36450. - 
         (-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.)*
          ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
            (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) + 
         (0.000027434842249657064 - sqcx2x4/72900.)*
          (-(sqx5*cx0x4)/72900. - (cx0x4*
               ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
            (sx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
          (-(sqx5*cx2x4)/72900. - (cx2x4*
               ((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
            (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))) - 
    ((0.000027434842249657064 - sqcx2x4/72900.)*
        (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
          (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
       (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
        ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
          (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))*
     (-(cx2x4*(cx0/8100. - (cx0x4*cx4)/24300.)*sx2x4)/36450. - 
       (cx2/12150. - (cx2x4*cx4)/24300.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.) + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(-(cx4*sx0x4)/24300. + (cx0x4*sx4)/24300.) - 
       (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.)) + 
    (-(cx0x4*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
       ((0.0000411522633744856 - sqcx0x4/72900.)*cx2x4*sx2x4)/36450. - 
       2*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(-(cx2x4*sx0x4)/72900. - (cx0x4*sx2x4)/72900.))*
     (-((cx2/12150. - (cx2x4*cx4)/24300.)*((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/
             270. - (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) + 
       (0.000027434842249657064 - sqcx2x4/72900.)*
        (-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
          (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.)) + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-((cx2/12150. - (cx2x4*cx4)/24300.)*(-(sqx5*cx2x4)/72900. - 
            (cx2x4*((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.))/270. - 
            (sx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)) - 
       ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
          (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)*
        (-(cx4*sx2x4)/24300. + (cx2x4*sx4)/24300.) + 
       (0.000027434842249657064 - sqcx2x4/72900.)*
        (-(sqx5*cx4)/24300. - (((sqx1*cx0x4)/270. + (sqx3*cx2x4)/270. - (109*cx4)/1000.)*cx4)/90. + 
          ((-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.)*sx4)/90.) - 
       (cx2x4*sx2x4*(-(cx4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/90. + 
            (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.))/36450.))/
  (-sqvar7 + 
    (-sqvar0 + 
       (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
     (-sqvar3 + 
       (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
  (*AA)( 7, 5 ) = -((-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
          (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
        ((xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx0x4)/36450. - 
          (xx[5]*(cx0x2/36450. - (cx0x4*cx2x4)/72900.)*sx2x4)/36450.)) + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-(xx[5]*(cx2/12150. - (cx2x4*cx4)/24300.)*sx2x4)/36450. - 
        (xx[5]*(0.000027434842249657064 - sqcx2x4/72900.)*sx4)/12150.))/
   (-pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
        (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2) + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-sqvar3 + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*AA)( 7, 6 ) = 0;
  (*AA)( 7, 7 ) = 0;

}

// Hard coded linearization of ff(xx, uu) about uu for points xx and uu. The equations were exported from Mathematica.
template<typename TT>
void PendCart::CalcBB(const vector<double> & xx, const vector<double> & uu, Eigen::MatrixBase<TT> * BB) {
  EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(TT, 8, 1);

  double cx0 = cos(xx[0]);
  double cx2 = cos(xx[2]);
  double cx4 = cos(xx[4]);


  double cx0x2 = cos(xx[0] - xx[2]);
  double cx0x4 = cos(xx[0] - xx[4]);
  double cx2x4 = cos(xx[2] - xx[4]);
  
  double sqcx4 = pow(cx4,2);
  
  double sqcx0x2 = pow(cx0x2,2);
  double sqcx0x4 = pow(cx0x4,2);
  double sqcx2x4 = pow(cx2x4,2);

  
  double sqvar0 = pow(cx0x2/36450. - (cx0x4*cx2x4)/72900.,2);
  double sqvar3 = pow(cx2/12150. - (cx2x4*cx4)/24300.,2);
  double sqvar9 = pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
        (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2);
  
  (*BB)( 0, 0 ) = 0;
  (*BB)( 1, 0 ) = ((-0.000027434842249657064 + sqcx2x4/72900.)*((0.000027434842249657064 - sqcx2x4/72900.)*
      (cx0/8100. - (cx0x4*cx4)/24300.) - (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
      (cx2/12150. - (cx2x4*cx4)/24300.)))/
 (270.*(-sqvar9 + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-sqvar3 + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*BB)( 2, 0 ) = 0;
  (*BB)( 3, 0 ) = ((-0.000027434842249657064 + sqcx2x4/72900.)*(-sqvar0 +
     (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
   ((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - (cx0*cx0x4*cx2x4)/300. - 
     (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.))/
 (90.*(-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - (cx0x2*cx0x4*cx2x4)/225. + 
     sqcx2x4/300.)*(-sqvar9 + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-sqvar3 + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*BB)( 4, 0 ) = 0;
  (*BB)( 5, 0 ) = ((-0.000027434842249657064 + sqcx2x4/72900.)*(-sqvar0 +
     (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
   ((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
     (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
        ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
          (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
      ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
          (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.)))))/
 (270.*(-sqvar9 + 
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-sqvar3 + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
  (*BB)( 6, 0 ) = 0;
  (*BB)( 7, 0 ) = -((-0.000027434842249657064 + sqcx2x4/72900.)*(-sqvar0 +
      (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.)))/
 (270.*(-sqvar9 +
     (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-sqvar3 + 
        (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
}
#endif // __Agile_RRT__pendcart__
