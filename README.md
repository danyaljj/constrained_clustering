# Constrained Clustering 

This is a MATLAB code containing a set of clustering algorithms. 

Part of this code is used to simulate the experiments in [this](http://arxiv.org/abs/1508.06235) work.

BTW, there is a list of constrained clustering algorithms with available codes [here](http://web.engr.illinois.edu/~khashab2/files/2015_constrained_clustering/constrainedClustering.html). 

## How to run: 

## Structure of this package
Here is how the code structured: 
- `algorithms` contains a the algorithms we have studied / experimented with, at some point. Many of these codes are downloaded from somewhere, and included directly (or with small modifications). In case there is a big modification, it is mentioned in the `README.md` inside their folder. 
Note that not all of these algorithms are used in the evaluation script (either due to instability, being slow, or not being compatible with our purposes). That said, you can always add these to the script and use them. 
- `data`: UCI data + toy data
- `distance`: some of the distance measures we have used across multiple algorithms here. 
- `metrics`: contains the evaluation metrics we have used.    

## Questions? 
Ask Daniel: [http://web.engr.illinois.edu/~khashab2/](http://web.engr.illinois.edu/~khashab2/)
