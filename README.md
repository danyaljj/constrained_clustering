# Constrained Clustering 

This is a MATLAB code containing a set of clustering algorithms. 

Part of this code is used to simulate experiments in [this](http://arxiv.org/abs/1508.06235) work.

Also there is a list of constrained clustering algorithms with available codes [here](http://web.engr.illinois.edu/~khashab2/files/2015_constrained_clustering/constrainedClustering.html). 

## How to run: 
To see output on toy data, go to the directory `experiment`, and run the script `experiment_toy.m`. You should be able to see the following output, followed by some other outputs: 
![alt text](https://github.com/danyaljj/constrained_clustering/blob/master/experiment/Gaussian-Mixtures_iter=1.tif)

You can run the script `experiment_uci.m` to see the output of the algorithms on the UCI dataset as well.  

## Structure of this package
Here is how the code structured: 
- `algorithms` contains a the algorithms we have studied / experimented with, at some point. Many of these codes are downloaded from somewhere, and included directly (or with small modifications). Some of these algorithms contain a `README.md` inside their folder, which explain where they are downloaded, and possible modifications / extensions on them.  
Note that not all of these algorithms are used in the evaluation script (either due to instability, being slow, or not being compatible with our purposes). That said, you can always add these to the script and use them. 
- `data`: UCI data + toy data
- `distance`: some of the distance measures we have used across multiple algorithms. 
- `metrics`: contains the evaluation metrics we have used.    

## Questions / Comments / Suggestions  
Email Daniel: [http://web.engr.illinois.edu/~khashab2/](http://web.engr.illinois.edu/~khashab2/)
