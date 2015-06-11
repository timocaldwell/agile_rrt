# AgileRRT for 3-Link Pendulum on a Cart
AgileRRT motion plans the highly dynamic 3-link pendulum on a cart. It executes the rapidly exploring random tree using LQR-based steering with linearization about the zero-control trajectory. 

Zero-control steering is presented in caldwell_isrr2015.pdf submitted to ISRR. The code implements Algorithm 4 in the paper, efficient fixed time horizon inexact linear steering, as well as Equation 26, the projection operation.

## Dependencies
Requires GSL Interpolation, Boost odeint, and Eigen.
* [http://www.gnu.org/software/gsl/](http://www.gnu.org/software/gsl/)
* [http://www.boost.org/doc/libs/1_58_0/libs/numeric/odeint/doc/html/index.html](http://www.boost.org/doc/libs/1_58_0/libs/numeric/odeint/doc/html/index.html)
* [http://eigen.tuxfamily.org/](http://eigen.tuxfamily.org/)

## Arguments
Arguments are specified in order. If one is missing, the default value is used (default value given right of '=' sign)
* t_h_max_upper = 1.0       | Max time horizon.
* usezerotraj = true        | Choose linearization for steering: true -> Linearization about zero-control trajectory; false -> Linearization about vertex x0.
* max_cnt = 1000            | Max number of successful insertions into the tree. RRT ends if max_cnt reached.
* max_miss = 1000           | Max number of failed insertions into the tree. RRT ends if max_miss reached.
* max_samplemiss = 10000    | Max number of failed attempts to find a
* numbruns = 1              | Number of runs of the RRT. Run more to get more statistics.
* stopdist = 2.0            | Stops RRT if a point is explored within a ball of radius stopdist from the goal state.
* printskip = 1             | Number of loops of the RRT to skip between printing details to standard out.
* nndelta = LARGENUM        | Nearest Neighbor ball radius for quickly checking
* stats_name = "0"          | Filename to store statistics as a csv. Use "0" or no arg to not make file.
* traj_name = "0"           | Filename to store tree trajectories easily readable by Mathematica. Use "0" or no arg to not make file
* seed = 0                  | To choose the seed to the random number generator
