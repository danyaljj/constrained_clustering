CODE FOR CONSTRAINED 1-SPECTRAL CLUSTERING
------------------------------------------

This archive contains a Matlab implementation of Constrained 1-Spectral Clustering 
as described in the paper
 
   Syama Sundar Rangapuram, Matthias Hein
   Constrained 1-Spectral Clustering   
   In International conference on Artificial Intelligence and Statistics (AISTATS), 2012. 
   Available here: http://jmlr.csail.mit.edu/proceedings/papers/v22/sundar12/sundar12.pdf
   Long version:

Current version: V1.1


INSTALLATION
------------

The implementation uses two mex-files to solve the (constrained) inner convex problem of Algorithm 1 (see paper). 
Compile them on the matlab command prompt: 

mex -largeArrayDims mex_solve_inner_problem.cpp
mex -largeArrayDims mex_solve_cnstr_inner_problem.cpp


QUICK DEMO
----------

Run test_cosc.m

This file loads the two moons dataset and runs cosc for k=2 with default options.
Here you will also get help on how to call cosc with various options, in particular how to use a faster version of cosc for both k=2, k>2; check the commented lines in this file.
You can load your dataset by replacing the first line.

DOCUMENTATION
-------------

Usage:
  [cut, clusters, viols, clusters_intermediate] = cosc(W, deg, k, ML, CL, start_flags, nInitializationsForBinarySplit, nRunsOfRecursiveSplit, perturbation, MAX_ITERS, verbosity)

Input
  W:      n x n similarity matrix, where n is the number of data points.
          It has to be a sparse symmetric matrix, with zeros on the diagonal.
  deg:    n x 1 degree vector (can also be seen as vertex weights)
          for Normalized cut it is sum(W,2)
          for Ratio cut, it is ones(size(W,1),1)
  k:      number of clusters
  ML:     m x 2 martix, specifying m must-link pairs
  CL:     p x 2 matrix, specifying p cannot-link pairs

  OPTIONAL arguments:
  start_flags:    0, don't use any other special initializations 
                  >=1, initialize the method with the second eigenvector of
                     the standard Graph Laplacian (Default for k>2)
                     We used this option in the paper for both k=2 and k>2
                  >=2, initialize with the solution of the method: Fast Normalized
                     Cuts with Linear Constraints (Default for k=2)
                     Note that this option is valid only for the binary parititioning case, i.e., k=2)
  nInitializationsForBinarySplit: 
                  how many total initializations for each binary split? We start the method from (nInitializationsForBinarySplit - start_flags) random
                  initializations.
                  Default: 10 for k=2, 2 for k>2.
  nRunsOfRecursiveSplit: 
                  Total number of re-runs of the complete recursive
                  split. This is valid only for k>2.
                  Default: 5 for k>2; for k=2, value passed in this argument is ignored.
  perturabation:  true/false. Perturbing the solution will lead to better
                  results but at the expense of higher runtime. (Currently
                  available only for binary  partitioning and by default true for k=2.)
  MAX_ITERS:      Maximum number of FISTA iterations for solving the inner problem.
                  Recommendation: 1000 - 10000, depending on the size of the problem.                 
  verbosity:      Display Level
                  0: Display nothing. (Default option)
                  1: Track the progress by displaying the main steps.
                  2: Display intermediate results (useful for preliminary sanity checks).
                  3: Display all details. Use this option to check if the maximum number
                     of FISTA iterations reached before the inner problem
                     is solved optimally.

Output
  cut:        cut value of the constrained cut problem
  clusters:   vector containing the clustering of the data points.
  viols:      number of violated constraints: zero for the binary
              partitioning problem
  clusters_intermediate:
              matrix of size n times k-1, containing in each column the intermediate clusters found in the recursive
              split. For example, clusters_intermediate(:, end) is same
              as clusters, clusters_intermedaite(:, 1) is the first
              2-way split and clusters_intermediate(:, i) gives
              clustering into i+1 components.
 
 
ADDITIONAL NOTES
----------------

1. You can compute the clustering error using the function:

     error = cluster_err(clusters, Y)

   where Y is the true label vector. Note that this is the clustering error obtained by majority voting!

2. The balanced cut of a (multi-) partition specified by clusters can be computed using:

     cut = bal_cut(W, deg, clusters)

   where 'deg' is the degree vector (Normalized cut), vector of ones (Ratio cut) or can be any other vertex weight vector.

3. At the moment we always assume that the constraints are consistent. 
   Some assert statements in the code would throw errors if the constraints are not consistent!


LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

If you use this code for your publication, please include a reference 
to the paper: S. Rangapuram and M. Hein. Constrained 1-Spectral Clustering, AISTATS 2012. 
 
 
CONTACT
srangapu@mpi-inf.mpg.de
hein@cs.uni-saarland.de

(C)2012 Syama Sundar Rangapuram, Matthias Hein
Machine Learning Group, Saarland University, Germany
http://www.ml.uni-saarland.de
