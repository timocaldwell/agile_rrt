//
//  tree.h
//  Agile_RRT
//
//  Created by Timothy Caldwell on 10/1/14.
//  Copyright (c) 2014 Timothy Caldwell. All rights reserved.
//

#ifndef __Agile_RRT__tree__
#define __Agile_RRT__tree__

#include "global.h"
#include "treevertex.h"


class Tree {
private:
    struct treenode
    {
       TreeVertex * vert;
       vector<treenode *> children;
       treenode * parent;
    };

    sampling_bnds bnds;
    int cnt, miss, samplemiss, max_cnt, max_miss, max_samplemiss, printskip;
    bool isprint, computedisttogoal, isstopdist;
    double disttogoal, stopdist, t_h_max_upper;
    
    treenode * root;
    treenode * leastdist_node;
    
    double nndelta;
    double seed;
    
    void print_tree (treenode & node, double dt, const string & name, ostream & stream);
    void print_branch (treenode & node, double dt, const string & name, ostream & stream);
    
    bool isnearer(treenode & node, double & best_dist, const NSx1_type & xxsamp, const double t_h);
    bool isxxzeronear(treenode & node, const NSx1_type & xxsamp, double delta, const double t_h);
    void nearestneighbor(treenode & node, treenode ** best_node_ptr, double & best_dist, const NSx1_type & xxsamp, const double delta, const double t_h);
    bool extend(treenode & node, const NSx1_type & xxsamp, const double t_h );
    void delete_tree( treenode * node );
    
public:
    Tree ( TreeVertex * init_vert, double t_h_max_upper_in, sampling_bnds bnds_in, int max_cnt_in = 10, int max_miss_in = 10, int max_samplemiss_in = 10, int printskip_in = 1, double seed_in = 0 );
    ~Tree() { delete_tree(root); }
    void set_xxgoal_params ( const NSx1_type & xxgoal );
    void set_xxgoal_stopdist (double stopdist_in, const NSx1_type & xxgoal);
    void set_nndelta (double nndelta_in) { nndelta = nndelta_in; }
    void print( double dt, const string & name = "tree", ostream & stream = cout);
    int get_cnttot() { return cnt; }
    int get_misstot() { return miss; }
    int get_samplemisstot() { return samplemiss; }
    double get_disttogoal() { return disttogoal; }

    void reset_tree ();
    void print_update(double t_h, const NSx1_type & xxsamp);
    void runrrt();
};

#endif /* defined(__Agile_RRT__tree__) */
