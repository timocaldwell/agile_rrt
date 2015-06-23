//
//  tree.cpp
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//

#include "tree.h"


void Tree::PrintTree(const treenode_ & node, double dt, const string & name, ostream * stream) {
  node.vert->PrintEdgeTraj(dt, name, stream);
  for (vector< treenode_ *>::const_iterator it = node.children.begin(); it!=node.children.end(); ++it)
    PrintTree(**it, dt, name, stream);
}

void Tree::PrintBranch(const treenode_ & node, double dt, const string & name, ostream * stream) {
  node.vert->PrintEdgeTraj(dt, name, stream);
  if (node.parent)
    PrintBranch(*node.parent, dt, name, stream);
}

bool Tree::IsNearer(treenode_ & node, double & best_dist, const NSx1_type & xxsamp, const double t_h) {
  double dist = node.vert->JJLinearSteer(xxsamp, t_h);
  if (dist < best_dist) {
    best_dist = dist;
    return true;
  }
  return false;
}

bool Tree::IsxxZeroNear(treenode_ & node, const NSx1_type & xxsamp, double delta, const double t_h) {
  ode_state_type xxzero_pt;
//  node.vert->get_xxzero_pt(t_h, xxzero_pt);  
  node.vert->get_xxzero_pt(t_h, &xxzero_pt);
  Eigen::Map<NSx1_type> xxzero_mat(xxzero_pt.data(), NS, 1);
  return sqrt((xxsamp-xxzero_mat).dot(xxsamp-xxzero_mat)) < delta;
}

void Tree::NearestNeighbor(treenode_ & node, const NSx1_type & xxsamp, const double delta, const double t_h, treenode_ ** best_node_ptr, double & best_dist) {
  if (IsxxZeroNear(node, xxsamp, delta, t_h))
    if (IsNearer(node, best_dist, xxsamp, t_h)) // updates best_dist
      *best_node_ptr = &node;
  
  //iterate through children
  for (vector< treenode_ *>::iterator it = node.children.begin(); it!=node.children.end(); ++it)
    NearestNeighbor(**it, xxsamp, delta, t_h, best_node_ptr, best_dist);
}

bool Tree::Extend(treenode_ & node, const NSx1_type & xxsamp, const double t_h ) {
  if (node.vert->Projection(xxsamp, t_h)) {
    // successful projection!
    TreeVertex * new_vert = new TreeVertex(node.vert->NewVertex());
    treenode_ * new_node = new treenode_;
    new_node->vert = new_vert;
    new_node->children.clear();
    new_node->parent = &node;
    node.children.push_back(new_node);
    return true;
  }
  return false;
}

void Tree::DeleteTree(treenode_ * node) {
  for (vector<treenode_*>::const_iterator it = node->children.begin(); it!=node->children.end(); ++it)
    DeleteTree(*it);
  delete node->vert;
  delete node;
}

Tree::Tree (TreeVertex * init_vert, double t_h_max_upper, sampling_bounds_struct bounds, int max_cnt, int max_miss, int max_samplemiss, int printskip, double seed)
    : t_h_max_upper_(t_h_max_upper), bounds_(bounds), max_cnt_(max_cnt), max_miss_(max_miss), max_samplemiss_(max_samplemiss), computedisttogoal_(false), printskip_(printskip), disttogoal_(kLARGENUM), cnt_(0), miss_(0), samplemiss_(0), isstopdist_(false), stopdist_(kLARGENUM), nndelta_(kLARGENUM), seed_(seed) {
  root_ = new treenode_;
  root_->vert = init_vert;
  root_->children.clear();
  root_->parent = nullptr;
  leastdist_node_ = root_;
}

void Tree::set_xxgoal_params (const NSx1_type & xxgoal) {
  computedisttogoal_ = true;
  root_->vert->set_xxgoal(xxgoal);
}

void Tree::set_xxgoal_stopdist (double stopdist, const NSx1_type & xxgoal) {
  set_xxgoal_params(xxgoal);
  stopdist_ = stopdist;
  isstopdist_ = true;
}

void Tree::Print(double dt, const string & name, ostream * stream) {
  PrintTree(*root_, dt, name, stream);
  if (computedisttogoal_)
    PrintBranch(*leastdist_node_, dt, "trajs to goal", stream);
}

void Tree::ResetTree () { // delete all nodes but root
  for (vector<treenode_*>::const_iterator it = root_->children.begin(); it!=root_->children.end(); ++it)
    DeleteTree(*it);
  root_->children.clear();
  disttogoal_ = kLARGENUM;
  cnt_ = 0;
  miss_ = 0;
  samplemiss_ = 0;
}

void Tree::PrintUpdate(double t_h, const NSx1_type & xxsamp) {
  cout << "-------------------------------------------------" << endl;
  cout << "    cnt " << cnt_ << " | misses " << miss_ << " | sample misses " << samplemiss_ << " | t_h " << t_h << endl;
  cout << "    xxsamp " << xxsamp.transpose() << endl;
//  if (computedisttogoal_)
//    cout << "--- best distance to goal is ---- " << disttogoal_ << endl;
}

void Tree::RunRRT() {
  ode_state_type xxsamp_vec(NS*1);
  Eigen::Map<NSx1_type> xxsamp(xxsamp_vec.data(), NS, 1);
  double t_h;
  vector<double> t_hs;

  double initialseed = seed_;
  if (seed_==0)
    initialseed = (double)time(0);
  default_random_engine rand_engine(initialseed); // TO SEED TO THE TIME
//        default_random_engine rand_engine((double)time(0)); // TO SEED TO THE TIME
//        default_random_engine rand_engine(100.1); // To seed fixed
  
  while (cnt_ < max_cnt_ && miss_ < max_miss_ && samplemiss_ < max_samplemiss_) {
    xxsamp_vec.clear();
    for (int ii = 0; ii < NS; ++ii)
    {
      uniform_real_distribution<double> unif(bounds_.xbounds(ii, 0), bounds_.xbounds(ii, 1));
      xxsamp_vec.push_back(unif(rand_engine));
    }
    uniform_real_distribution<double> unif(0,t_h_max_upper_);
    t_h = unif(rand_engine);

    treenode_ * best_node = nullptr;
    double best_dist = kLARGENUM;
    NearestNeighbor(*root_, xxsamp, nndelta_, t_h, &best_node, best_dist);
    if (best_node) { // check that a best_node is found. If not, continue
      if (Extend(*best_node, xxsamp, t_h)) {
        cnt_++;
        t_hs.push_back(t_h);
        if (computedisttogoal_) {
          if (best_node->vert->get_disttogoal() < disttogoal_) {
            disttogoal_ = best_node->vert->get_disttogoal();
            // new leastdist_node_ is the newly created node from Extend.
            cout << "           *** new best distance to goal is *** " << disttogoal_ << endl;
            leastdist_node_ = best_node->children.back();

//            ode_state_type * foundxx = nullptr;
//            leastdist_node_->vert->get_x0(foundxx);
//            cout << "new best state : ";
//            PrintVec(*foundxx);
//            cout << endl;

            if (isstopdist_ && disttogoal_ < stopdist_)
              break;
          }
        }
        if ((cnt_) % printskip_ == 0) {//( (miss + cnt) % printskip == 0) {
          PrintUpdate(t_h, xxsamp);
        }
      }
      else miss_++;
    }
    else samplemiss_++;
  }
  PrintUpdate(t_h, xxsamp);
  if (computedisttogoal_)
    cout << "--- best distance to goal is ---- " << disttogoal_ << endl;
}

