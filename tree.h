//
//  tree.h
//  Agile_RRT
//
//  Created by Timothy M. Caldwell.
//  Copyright (c) 2015 Timothy M. Caldwell. Refer to license.txt
//

#ifndef __Agile_RRT__tree__
#define __Agile_RRT__tree__

#include "global.h"
#include "treevertex.h"


class Tree {
 public:
  Tree(TreeVertex * init_vert, double t_h_max_upper, sampling_bounds_struct bounds, int max_cnt = 10, int max_miss = 10, int max_samplemiss = 10, int printskip = 1, double seed = 0);
  ~Tree() { DeleteTree(root_); }
  void set_xxgoal_params(const NSx1_type & xxgoal);
  void set_xxgoal_stopdist(double stopdist, const NSx1_type & xxgoal);
  void set_nndelta(double nndelta) { nndelta_ = nndelta; }
  void Print(double dt, const string & name = "tree", ostream * stream = (&cout));
  int get_cnttot() {return cnt_;}
  int get_misstot() {return miss_;}
  int get_samplemisstot() {return samplemiss_;}
  double get_disttogoal() {return disttogoal_;}

  void ResetTree();
  void PrintUpdate(double t_h, const NSx1_type & xxsamp);
  void RunRRT();
  
 private:
  struct treenode_
  {
   TreeVertex * vert;
   vector<treenode_ *> children;
   treenode_ * parent;
  };

  sampling_bounds_struct bounds_;
  int cnt_, miss_, samplemiss_, max_cnt_, max_miss_, max_samplemiss_, printskip_;
  bool isprint_, computedisttogoal_, isstopdist_;
  double disttogoal_, stopdist_, t_h_max_upper_;
  
  treenode_ * root_;
  treenode_ * leastdist_node_;
  
  double nndelta_;
  double seed_;
  
  void PrintTree(const treenode_ & node, double dt, const string & name, ostream * stream);
  void PrintBranch(const treenode_ & node, double dt, const string & name, ostream * stream);
  
  bool IsNearer(treenode_ & node, double & best_dist, const NSx1_type & xxsamp, const double t_h);
  bool IsxxZeroNear(treenode_ & node, const NSx1_type & xxsamp, double delta, const double t_h);
  void NearestNeighbor(treenode_ & node, const NSx1_type & xxsamp, const double delta, const double t_h, treenode_ ** best_node_ptr, double & best_dist);
  bool Extend(treenode_ & node, const NSx1_type & xxsamp, double t_h);
  void DeleteTree(treenode_ * node);
};

#endif // __Agile_RRT__tree__
