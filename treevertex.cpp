//
//  treevertex.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 5/10/15.
//  Copyright (c) 2015 Timothy Caldwell. All rights reserved.
//

#include "treevertex.h"

TreeVertex::TreeVertex (ode_state_type x0, bool usezerotraj, double t_h_max_upper, InterpVector * xxedge, InterpVector * uuedge, const kin_constraints & constraints, const NSxNS_type & QQ, const NIxNI_type & RR, const NSxNS_type & P1, const NSxNS_type & QQlqr, const NIxNI_type & RRlqr, const NSxNS_type & P1lqr, double t_h_bar, bool save_edge)
    : usezerotraj_(usezerotraj), t_h_max_upper_(t_h_max_upper), QQlqr_(QQlqr), QQ_(QQ), RR_(RR), RRinv_(RR.inverse()), RRlqr_(RRlqr), RRlqr_inv_(RRlqr.inverse()), P1_(P1), P1lqr_(P1lqr), constraints_(constraints), t_h_bar_(t_h_bar), computedisttogoal_(false), xxgoal_(NSx1_type::Zero()), disttogoal_cur_(kLARGENUM), save_edge_(save_edge)
{
  t_h_cur_ = kLARGENUM;
  JJapprox_cur_ = kLARGENUM;
  JJproj_cur_ = kLARGENUM;
  x0_ = x0;
  
  xxproj_cur_ = nullptr;
  uuproj_cur_ = nullptr;

  if (save_edge) {
    xxedge_ = new InterpVector(*xxedge); // deep copy
    uuedge_ = new InterpVector(*uuedge);
  } else {
    xxedge_ = nullptr;
    uuedge_ = nullptr;
  }
  Initialize();
}

TreeVertex::~TreeVertex()
{
  if ( xxproj_cur_ )
    delete xxproj_cur_;
  if ( uuproj_cur_)
    delete uuproj_cur_;
  if ( xxedge_ )
    delete xxedge_;
  if ( uuedge_ )
    delete uuedge_;
  delete xxzero_;
  delete AA_;
  delete BB_;
  delete KK_;
  delete WWK_;
  delete SSK_;
}

void TreeVertex::Initialize()
{
  vector<double> tts;
  vector< ode_state_type > xxs;
//  cout << "a" <<endl;
  if (usezerotraj_) {
    //------------- FOR Linearizing about zero trajectory!--------------
    if( !sys::IntegrateFreeDynamics(x0_, t_h_max_upper_, &tts, &xxs) )
      cout << "ERROR in treevertex: Zero dynamics integration failed!";
    xxzero_ = new InterpVector(tts, xxs);
  } else {
    //------------- FOR Linearizing about point!--------------
    tts.push_back(0.0);
    tts.push_back(t_h_max_upper_/2.0);
    tts.push_back(t_h_max_upper_);
    xxs.push_back(x0_);
    xxs.push_back(x0_);
    xxs.push_back(x0_);
    xxzero_ = new InterpVector(tts,xxs);
  }
//  cout << "b" <<endl;
  PendCart pendcart;
  vector< double > uuzero(NI, 0.0);
  vector< vector<double> > AAs(tts.size(),vector<double>(NS*NS));
  vector< vector<double> > BBs(tts.size(),vector<double>(NS*NI));
  for(int ii = 0; ii < tts.size(); ++ii) {
    Eigen::Map< NSxNS_type > AA_mat(AAs[ii].data(), NS, NS);
    pendcart.CalcAA(xxs[ii], uuzero, &AA_mat);
    
    Eigen::Map< NSxNI_type > BB_mat(BBs[ii].data(), NS, NI);
    pendcart.CalcBB(xxs[ii], uuzero, &BB_mat);
  }
  AA_ = new InterpVector(tts, AAs);
  BB_ = new InterpVector(tts, BBs);
  
  tts.clear();
  vector<ode_state_type > PPs;
  ode_state_type P1lqr_vec(NS*NS);
  Eigen::Map<NSxNS_type> P1lqr_mat(P1lqr_vec.data(), NS, NS);
  P1lqr_mat = P1lqr_;

  if( !sys::IntegrateLTVRicatti(*AA_, *BB_, QQlqr_, RRlqr_inv_, P1lqr_vec, t_h_max_upper_, &tts, &PPs) )
    cout << "ERROR in treevertex: Ricatti integration failed!";
//  cout << "c" <<endl;
  vector<vector<double>> KKs(tts.size(),vector<double>(NI*NS));

  for(int ii = 0; ii < tts.size(); ++ii) {
    Eigen::Map< NSxNS_type > PP_mat(PPs[ii].data(), NS, NS);

    vector<double> BB_vec(NS*NI);
    BB_->pt(tts[ii], &BB_vec);
    Eigen::Map<NSxNI_type> BB_mat(BB_vec.data(), NS, NI);

    Eigen::Map<NIxNS_type> KK_mat(KKs[ii].data(), NI, NS);
    KK_mat = RRinv_ * BB_mat.transpose() * PP_mat;
  }
  KK_ = new InterpVector(tts, KKs);
//  cout << "d" <<endl;
  tts.clear();
  vector< ode_state_type > WWKs;
  if( !sys::IntegrateWWK(*AA_, *BB_, *KK_, RRinv_, t_h_max_upper_, &tts, &WWKs) )
    cout << "ERROR in treevertex: WWK integration failed!";
  WWK_ = new InterpVector(tts, WWKs);
//  cout << "e" <<endl;
  tts.clear();
  vector< ode_state_type > SSKs;
  if( !sys::IntegrateSSK(*AA_, *BB_, *KK_, *WWK_, RR_, RRinv_, t_h_max_upper_, &tts, &SSKs) )
    cout << "ERROR in treevertex: SSK integration failed!";
  SSK_ = new InterpVector(tts, SSKs);
//  cout << "f" <<endl;
}


double TreeVertex::JJLinearSteer (const NSx1_type & xxsamp, const double t_h)
{
  if (t_h > t_h_bar_)
    return kLARGENUM;
//  cout << "g" <<endl;
  t_h_cur_ = t_h;
  xxsamp_cur_ = xxsamp;
  
  vector<double> xxzero_vec(NS*1);
  xxzero_->pt(t_h, &xxzero_vec);
  Eigen::Map<NSx1_type> xxzero_mat(xxzero_vec.data(), NS, 1);
//  cout << "h" <<endl;
  vector<double> WWK_vec(NS*NS);
  WWK_->pt(t_h, &WWK_vec);
  Eigen::Map<NSxNS_type> WWK_mat(WWK_vec.data(), NS, NS);
//  cout << "i" <<endl;
  vector<double> SSK_vec(NS*NS);
  SSK_->pt(t_h, &SSK_vec);
  Eigen::Map<NSxNS_type> SSK_mat(SSK_vec.data(), NS, NS);
//  cout << "j" <<endl;
  Pth_cur_ = (WWK_mat*P1_*WWK_mat + SSK_mat).inverse()*WWK_mat*P1_;
  xxzeroerror_cur_ = xxzero_mat - xxsamp;
  eta_cur_ = Pth_cur_ * xxzeroerror_cur_;
//  cout << "k" <<endl;
  JJapprox_cur_ = 0.5*eta_cur_.dot(SSK_mat*eta_cur_) + 0.5*(xxzeroerror_cur_ - WWK_mat*eta_cur_).dot(P1_*(xxzeroerror_cur_ - WWK_mat*eta_cur_));
  return JJapprox_cur_;
}

void TreeVertex::set_xxgoal (const NSx1_type & xxgoal)
{
  xxgoal_ = xxgoal;
  computedisttogoal_ = true;
}

bool TreeVertex::Projection(const NSx1_type & xxsamp, const double t_h)
{
  if(xxsamp_cur_!=xxsamp || t_h_cur_!=t_h) // if linear steering computations have not already been made:
    JJLinearSteer(xxsamp, t_h);
//  cout << "l" <<endl;
  vector<double> tts;
  vector< ode_state_type > Phis;
  if ( !sys::IntegrateLTVClosedLoopSTMSlot2(*AA_, *BB_, *KK_, t_h, &tts, &Phis) )
    cout << "ERROR in linear_steer: STM integration failed!";
//  cout << "m" <<endl;
  InterpVector Phi(tts, Phis);// = new InterpVector(tts, Phis);

  tts.clear();
  vector< ode_state_type > xxs;
  vector< ode_state_type > uus;
//  cout << "n" <<endl;
  if ( !sys::IntegrateTrajectoryLinearSteeringProjectionWithConstraintsWithCost(xxsamp, *xxzero_, *BB_, *KK_, *KK_, *WWK_, Phi, eta_cur_, QQ_, RR_, P1_, RRinv_, x0_, constraints_, t_h, &tts, &xxs, &uus, &JJproj_cur_) ) {
//    cout << "printing shit" << endl;
//    PrintVec(uus[0]);
    return false;
  }
//  cout << endl;
  if(computedisttogoal_) {
    disttogoal_cur_ = kLARGENUM;
    for (int ii = 0; ii < xxs.size(); ++ii) {
      Eigen::Map<NSx1_type> xx_mat(xxs[ii].data(), NS, 1);
      double tmp  = 0.5 * (xxgoal_-xx_mat).dot(P1_*(xxgoal_-xx_mat));
      if ( tmp < disttogoal_cur_ )
        disttogoal_cur_ = tmp;
    }
  }
//  cout << "o" <<endl;
  //remove old projections if needed.
  if (xxproj_cur_)
    delete xxproj_cur_;
  if (uuproj_cur_)
    delete uuproj_cur_;
//  cout << "p" <<endl;
  xxproj_cur_ = new InterpVector(tts, xxs);
  uuproj_cur_ = new InterpVector(tts, uus);

//  uuproj_cur_->Print(.1);
//  cout << "q" <<endl;
  return true;
}

void TreeVertex::PrintTraj(double dt, const string & name, ostream & stream)
{
  xxproj_cur_->Print(dt, name+"xxproj", stream);
  uuproj_cur_->Print(dt, name+"uuproj", stream);
}

void TreeVertex::PrintEdgeTraj(double dt, const string & name, ostream & stream)
{
  xxedge_->Print(dt, name+"xxedge", stream);
  uuedge_->Print(dt, name+"uuedge", stream);
}

TreeVertex TreeVertex::NewVertex()
{
  ode_state_type x0_new;
  xxproj_cur_->pt(t_h_cur_, &x0_new);
  
  //TODO (Tim): make a angle notifier for angle wrapping general for problems
  ///////////////////////////////
  //
  //SPECIAL CODE FOR THE PENDULUM FOR ANGLE WRAPPING
  //
  ///////////////////////////////
  for ( int ii = 0; ii < NS/2-1; ++ii ) {
    if ( x0_new[2*ii] < -M_PI && x0_new[2*ii+1] < 0 )
      x0_new[2*ii] += 2*M_PI;
    else if ( x0_new[2*ii] > M_PI && x0_new[2*ii+1] > 0 )
      x0_new[2*ii] -= 2*M_PI;
  }
  ///////////////////////////////
  //
  //END SPECIAL CODE FOR THE PENDULUM
  //
  ///////////////////////////////
  
  TreeVertex new_vert(x0_new, usezerotraj_, t_h_max_upper_, xxproj_cur_, uuproj_cur_, constraints_, QQ_, RR_, P1_, QQlqr_, RRlqr_, P1lqr_, t_h_bar_, save_edge_);
  if (computedisttogoal_)
    new_vert.set_xxgoal(xxgoal_);
  return new_vert;
}
