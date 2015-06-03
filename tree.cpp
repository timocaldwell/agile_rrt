//
//  tree.cpp
//  Agile_RRT
//
//  Created by Timothy Caldwell on 5/10/15.
//  Copyright (c) 2015 Timothy Caldwell. All rights reserved.
//

#include "tree.h"


void Tree::print_tree (treenode & node, double dt, const string & name, ostream & stream) {
    node.vert->print_edge_traj(dt, name, stream);
    for (vector< treenode *>::const_iterator it = node.children.begin(); it!=node.children.end(); ++it)
        print_tree(**it, dt, name, stream);
}

void Tree::print_branch (treenode & node, double dt, const string & name, ostream & stream) {
    node.vert->print_edge_traj(dt, name, stream);
    if ( node.parent )
        print_branch(*node.parent, dt, name, stream);
}

bool Tree::isnearer(treenode & node, double & best_dist, const NSx1_type & xxsamp, const double t_h) {
    double dist = node.vert->JJ_linear_steer(xxsamp, t_h);
    if ( dist < best_dist ) {
        best_dist = dist;
        return true;
    }
    return false;
}

bool Tree::isxxzeronear(treenode & node, const NSx1_type & xxsamp, double delta, const double t_h) {
    ode_state_type xxzero_pt;
    node.vert->get_xxzero_pt(t_h, xxzero_pt);
    Eigen::Map<NSx1_type> xxzero_mat(xxzero_pt.data(), NS, 1);
    return sqrt((xxsamp-xxzero_mat).dot(xxsamp-xxzero_mat)) < delta;
}

void Tree::nearestneighbor(treenode & node, treenode ** best_node_ptr, double & best_dist, const NSx1_type & xxsamp, const double delta, const double t_h) {
    if ( isxxzeronear(node, xxsamp, delta, t_h) )
        if ( isnearer(node, best_dist, xxsamp, t_h) ) // updates best_dist
            *best_node_ptr = &node;
    
    
    //iterate through children
    for (vector< treenode *>::iterator it = node.children.begin(); it!=node.children.end(); ++it)
        nearestneighbor(**it, best_node_ptr, best_dist, xxsamp, delta, t_h);
}

bool Tree::extend(treenode & node, const NSx1_type & xxsamp, const double t_h ) {
    if ( node.vert->projection(xxsamp, t_h) ) {
        // successful projection!
        TreeVertex * new_vert = new TreeVertex(node.vert->new_vertex());
        treenode * new_node = new treenode;
        new_node->vert = new_vert;
        new_node->children.clear();
        new_node->parent = &node;
        node.children.push_back(new_node);
        return true;
    }
    return false;
}

void Tree::delete_tree( treenode * node ) {
    for (vector< treenode *>::const_iterator it = node->children.begin(); it!=node->children.end(); ++it)            delete_tree(*it);
    delete node->vert;
    delete node;
}

Tree::Tree ( TreeVertex * init_vert, double t_h_max_upper_in, sampling_bnds bnds_in, int max_cnt_in, int max_miss_in, int max_samplemiss_in, int printskip_in, double seed_in ) : t_h_max_upper(t_h_max_upper_in), bnds(bnds_in), max_cnt(max_cnt_in), max_miss(max_miss_in), max_samplemiss(max_samplemiss_in), computedisttogoal(false), printskip(printskip_in), disttogoal(LARGENUM), cnt(0), miss(0), samplemiss(0), isstopdist(false), stopdist(LARGENUM), nndelta(LARGENUM), seed(seed_in) {
    
    root = new treenode;
    root->vert = init_vert;
    root->children.clear();
    root->parent = nullptr;
    leastdist_node = root;
}

void Tree::set_xxgoal_params ( const NSx1_type & xxgoal ) {
    computedisttogoal = true;
    root->vert->set_xxgoal(xxgoal);
}

void Tree::set_xxgoal_stopdist (double stopdist_in, const NSx1_type & xxgoal) {
    set_xxgoal_params(xxgoal);
    stopdist = stopdist_in;
    isstopdist = true;
}

void Tree::print( double dt, const string & name, ostream & stream ) {
    print_tree (*root, dt, name, stream);
    if ( computedisttogoal )
        print_branch ( *leastdist_node, dt, "trajs to goal", stream);
}

void Tree::reset_tree () { // delete all nodes but root
    for (vector< treenode *>::const_iterator it = root->children.begin(); it!=root->children.end(); ++it)            delete_tree(*it);
    root->children.clear();
    disttogoal = LARGENUM;
    cnt = 0;
    miss = 0;
    samplemiss = 0;
}

void Tree::print_update(double t_h, const NSx1_type & xxsamp) {
    cout << "-------------------------------------------------" << endl;
    cout << "    cnt " << cnt << " | misses " << miss << " | sample misses " << samplemiss << " | t_h " << t_h << endl;
    cout << "    xxsamp " << xxsamp.transpose() << endl;
    if ( computedisttogoal )
        cout << "--- best distance to goal is ---- " << disttogoal << endl;
}

void Tree::runrrt() {
    ode_state_type xxsamp_vec(NS*1);
    Eigen::Map<NSx1_type> xxsamp(xxsamp_vec.data(), NS, 1);
    double t_h;
    vector<double> t_hs;
    
    time_t cur_time = time(0);
    default_random_engine rand_engine((double)cur_time * seed); // TO SEED TO THE TIME
//        default_random_engine rand_engine((double)time(0)); // TO SEED TO THE TIME
//        default_random_engine rand_engine(100.1); // To seed fixed
    
//        disttogoal = LARGENUM;
//        cnt = 0;
//        miss = 0;
//        samplemiss = 0;
    
    while (cnt < max_cnt && miss < max_miss && samplemiss < max_samplemiss) {//(cnt + miss < max_cnt) {
        //random state and time horizon
        xxsamp_vec.clear();
        for ( int ii = 0; ii < NS; ++ii)
        {
            uniform_real_distribution<double> unif(bnds.xbnds(ii, 0), bnds.xbnds(ii, 1));
            xxsamp_vec.push_back( unif(rand_engine) );
        }
        uniform_real_distribution<double> unif(0,t_h_max_upper);
        t_h = unif(rand_engine);
    
        treenode * best_node = nullptr;
        double best_dist = LARGENUM;
//        double DELTA = LARGENUM;//10;
        nearestneighbor(*root, &best_node, best_dist, xxsamp, nndelta, t_h);
        if ( best_node ) { // check that a best_node is found. If not, continue
            if ( extend(*best_node, xxsamp, t_h) ) {
                cnt++;
                t_hs.push_back(t_h);
                if ( computedisttogoal ) {
                    if ( best_node->vert->get_disttogoal() < disttogoal ) {
                        disttogoal = best_node->vert->get_disttogoal();
                        leastdist_node = best_node;
                        if ( isstopdist && disttogoal < stopdist )
                            break;
                    }
                }
                if ( (cnt) % printskip == 0 ) {//( (miss + cnt) % printskip == 0) {
                    print_update(t_h, xxsamp);
                }
            }
            else miss++;
        }
        else samplemiss++;
    }
    print_update(t_h, xxsamp);
    if ( computedisttogoal )
        cout << "--- best distance to goal is ---- " << disttogoal << endl;
}

