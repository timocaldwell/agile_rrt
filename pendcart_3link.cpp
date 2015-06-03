//
//  pendcart.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 5/9/15.
//  Copyright (c) 2015 Timothy Caldwell. All rights reserved.
//

#include "pendcart_3link.h"

void PendCart::operator()( const ode_state_type xx , ode_state_type &dxxdtt , const double tt ) {
//    xx = {theta1 , theta1_dot , theta2 , theta2_dot , theta3 , theta3_dot , X , X_dot}
    if (type == FEEDFORWARD_TYPE)
        uu_ff_interp->pt(tt, uu);
    else if (type == TRACKPOINT_TYPE)
    {
        vector<double> KK_vec(NI*NS);
        KK_interp->pt(tt, KK_vec);
        
        Eigen::Map<const NIxNS_type> KK_mat(KK_vec.data(), NI, NS);
        Eigen::Map<const NSx1_type> xx_ref_mat(xx_ref_pt.data(), NS, 1);
        Eigen::Map<const NSx1_type> xx_mat(xx.data(), NS, 1);
        uu.resize(NI);
        Eigen::Map<NIx1_type> uu_mat(uu.data(), NI, 1);
        uu_mat = KK_mat * (xx_ref_mat - xx_mat);
    }
    else if (type == TRACKTRAJECTORY_TYPE)
    {
        vector<double> KK_vec(NI*NS);
        KK_interp->pt(tt, KK_vec);
        Eigen::Map<const NIxNS_type> KK_mat(KK_vec.data(), NI, NS);

        vector<double> xx_ref_vec(NS);
        xx_ref_interp->pt(tt, xx_ref_vec);
        Eigen::Map<const NSx1_type> xx_ref_mat(xx_ref_vec.data(), NS, 1);
        
        vector<double> uu_ff_vec(NI);
        uu_ff_interp->pt(tt, uu_ff_vec);
        Eigen::Map<const NIx1_type> uu_ff_mat(uu_ff_vec.data(), NI, 1);

        Eigen::Map<const NSx1_type> xx_mat(xx.data(), NS, 1);
        uu.resize(NI);
        Eigen::Map<NIx1_type> uu_mat(uu.data(), NI, 1);

        uu_mat = uu_ff_mat + KK_mat * (xx_ref_mat - xx_mat);
        
        *uu_out = uu;
    }
    else if (type == LINEARSTEERINGPROJECTION_TYPE || type == LINEARSTEERINGPROJECTIONWCOST_TYPE)
    {
        vector<double> xxzero_vec(NS*1);
        xxzero->pt(tt, xxzero_vec);
        Eigen::Map<const NSx1_type> xxzero_mat(xxzero_vec.data(), NS, 1);
        
        vector<double> BB_vec(NS*NI);
        BB->pt(tt, BB_vec);
        Eigen::Map<const NSxNI_type> BB_mat(BB_vec.data(), NS, NI);
    
        vector<double> KKlin_vec(NI*NS);
        KKlin->pt(tt, KKlin_vec);
        Eigen::Map<const NIxNS_type> KKlin_mat(KKlin_vec.data(), NI, NS);
        
        vector<double> KKproj_vec(NI*NS);
        KKproj->pt(tt, KKproj_vec);
        Eigen::Map<const NIxNS_type> KKproj_mat(KKproj_vec.data(), NI, NS);
        
        vector<double> WWK_vec(NS*NS);
        WWK->pt(tt, WWK_vec);
        Eigen::Map<const NSxNS_type> WWK_mat(WWK_vec.data(), NS, NS);

        vector<double> Phi_vec(NS*NS);
        Phi->pt(tt, Phi_vec);
        Eigen::Map<const NSxNS_type> Phi_mat(Phi_vec.data(), NS, NS);
        
        NSx1_type xxtilde_mat = xxzero_mat - WWK_mat*Phi_mat.transpose()*eta;
        NIx1_type uutilde_mat = (KKlin_mat*WWK_mat - RRinv*BB_mat.transpose())*Phi_mat.transpose()*eta;
        
        Eigen::Map<const NSx1_type> xx_mat(xx.data(), NS, 1);
        uu.resize(NI);
        Eigen::Map<NIx1_type> uu_mat(uu.data(), NI, 1);

        uu_mat = uutilde_mat - KKproj_mat*(xx_mat - xxtilde_mat);
        *uu_out = uu;
        
        if ( type == LINEARSTEERINGPROJECTIONWCOST_TYPE ) // integrate ell
            dxxdtt.back() = 0.5*xx_mat.dot(QQ*xx_mat) + 0.5*uu_mat.dot(RR*uu_mat);
    
    }
    
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

    double sqvar0 = pow(cx0x2/36450. - (cx0x4*cx2x4)/72900.,2);
    double sqvar1 = pow(cx2/12150. - (cx2x4*cx4)/24300.,2);
    
    dxxdtt[0] = xx[1];
    dxxdtt[1] = -(((0.000027434842249657064 - sqcx2x4/72900.)*
         (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
           (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
        (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
         ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
           (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.))/
      (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) + 
   (((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
        (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.))*
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
    ((-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
      (-pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2) + 
        (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
         (-sqvar1 + 
           (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
    dxxdtt[2] = xx[3];
    dxxdtt[3] = ((cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
      ((0.000027434842249657064 - sqcx2x4/72900.)*
         (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
           (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
        (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
         ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
           (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
    ((0.000027434842249657064 - sqcx2x4/72900.)*
      (-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
   (3*((sqx1*sx0x2)/45. + (327*sx2)/500. - (sqx1*cx2x4*sx0x4)/90. - 
        (sqx5*sx2x4)/90. - (sqx3*cx2x4*sx2x4)/90. - (327*cx2x4*sx4)/1000.))/
    (-0.06666666666666667 + sqcx2x4/30.) + 
   (3*((cx0*cx0x2)/150. - cx2/150. + (cx2*sqcx0x4)/450. - 
        (cx0*cx0x4*cx2x4)/300. - (cx0x2*cx0x4*cx4)/450. + (cx2x4*cx4)/300.)*
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
    ((-0.006666666666666667 + sqcx0x2/225. + sqcx0x4/450. - 
        (cx0x2*cx0x4*cx2x4)/225. + sqcx2x4/300.)*
      (-pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
           (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2) + 
        (-sqvar0 + 
           (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
         (-sqvar1 + 
           (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));
    dxxdtt[4] = xx[5];
    dxxdtt[5] = (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
      ((0.000027434842249657064 - sqcx2x4/72900.)*
         (((-327*sx0)/1000. + (sqx3*sx0x2)/135. + (sqx5*sx0x4)/270.)/270. - 
           (cx0x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.) - 
        (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*
         ((-(sqx1*sx0x2)/135. - (109*sx2)/500. + (sqx5*sx2x4)/270.)/270. - 
           (cx2x4*(-(sqx1*sx0x4)/270. - (sqx3*sx2x4)/270. - (109*sx4)/1000.))/270.)))/
    ((-0.06666666666666667 + sqcx2x4/30.)*(-sqvar0 + 
        (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))) - 
   3*(-(sqx1*sx0x4)/3. - (sqx3*sx2x4)/3. - (981*sx4)/100.) + 
   (3*cx2x4*((sqx1*sx0x2)/45. + (327*sx2)/500. - (sqx1*cx2x4*sx0x4)/90. - 
        (sqx5*sx2x4)/90. - (sqx3*cx2x4*sx2x4)/90. - (327*cx2x4*sx4)/1000.))/
    (-0.06666666666666667 + sqcx2x4/30.) + 
   (((3*((cx2*cx2x4)/15. - cx4/15.))/(-0.06666666666666667 + sqcx2x4/30.) - 
        (3*(-cx0x4/45. + (cx0x2*cx2x4)/45.)*
           ((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
             (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.)))/
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
    (-pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
         (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2) + 
      (-sqvar0 + 
         (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
       (-sqvar1 + 
         (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.)));
    dxxdtt[6] = xx[7];
    dxxdtt[7] = -((-(((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) - 
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
             (-uu[0] - (sqx1*sx0)/30. - (sqx3*sx2)/45. - (sqx5*sx4)/90.)/270.)))/
     (-pow((0.000027434842249657064 - sqcx2x4/72900.)*(cx0/8100. - (cx0x4*cx4)/24300.) -
          (cx0x2/36450. - (cx0x4*cx2x4)/72900.)*(cx2/12150. - (cx2x4*cx4)/24300.),2) + 
       (-sqvar0 + 
          (0.0000411522633744856 - sqcx0x4/72900.)*(0.000027434842249657064 - sqcx2x4/72900.))*
        (-sqvar1 + 
          (0.000027434842249657064 - sqcx2x4/72900.)*(0.004074074074074074 - sqcx4/8100.))));

}

bool PendCart::is_constraintsatisfied( vector<double> & xx_prev, vector<double> & xx_cur, vector<double> & uu_cur ) const {
    for (int ii = 0; ii < NS; ++ii)
    {
        if ( xx_cur[ii] < constraints.xbnds( ii, 0 ) )
            return false;
        else if ( xx_cur[ii] > constraints.xbnds( ii, 1 ) )
            return false;
    }

    for (int ii = 0; ii < NI; ++ii)
    {
        if ( uu_cur[ii] < constraints.ubnds( ii, 0 ))
            return false;
        if ( uu_cur[ii] > constraints.ubnds( ii, 1 ))
            return false;
    }

    double qqX = xx_prev[NS-2];
    double rrX = xx_cur[NS-2];
    double qqY = 0;
    double rrY = 0;
    for (int jj = 0; jj < numlinks; ++jj) {
        qqX += linklen * sin(xx_prev[jj*2]);
        qqY += linklen * cos(xx_prev[jj*2]);
        rrX += linklen * sin(xx_cur[jj*2]);
        rrY += linklen * cos(xx_cur[jj*2]);

        for (int ii = 0; ii < constraints.obs_radii.size(); ++ii)
        { // is_constraintsatisfied = false when the line segment between the pendulum head of state xx (i.e. point rr) and xx_prev (i.e. point qq) breaks disk constraint of being obs_radii from point pp.
            double sqradius = pow(constraints.obs_radii[ii],2);
            double ppX = constraints.obs_Xs[ii];
            double ppY = constraints.obs_Ys[ii];
            if ( is_linecircleintersect(sqradius, ppX, ppY, qqX, qqY, rrX, rrY) )
                return false;
//            cout << "here " << qqX << " " << qqY << " " << is_linecircleintersect(sqradius, ppX, ppY, qqX, qqY, rrX, rrY) << endl;
        }
    }
//    cout << endl;
    return true;
}

bool PendCart::is_linecircleintersect(double sqradius, double ppX, double ppY, double qqX, double qqY, double rrX, double rrY) const {
    double rrqq_distance_sqaured = ((rrX-qqX)*(rrX-qqX) + (rrY-qqY)*(rrY-qqY));
    double ss;
    if (rrqq_distance_sqaured == 0) // e.g. same point
        ss = 0; // i.e. pick one.
    else
        ss = ((ppX-qqX)*(rrX-qqX) + (ppY-qqY)*(rrY-qqY)) / rrqq_distance_sqaured;
    if ( ss <= 0.0 && pow(qqX-ppX,2) + pow(qqY-ppY,2) < sqradius )
        return true;
    if ( ss >= 1.0 && pow(rrX-ppX,2) + pow(rrY-ppY,2) < sqradius )
        return true;
    double llX = (rrX - qqX)*ss + qqX;
    double llY = (rrY - qqY)*ss + qqY;
    if ( ss >= 0.0 && ss <= 1.0 && pow(llX-ppX,2) + pow(llY-ppY,2) < sqradius )
        return true;
    return false;
}


