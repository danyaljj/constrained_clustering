% TEXTCLUST toolbox
% A model-based document clustering toolbox.
%
% Author: Shi Zhong, ECE Dept., The University of Texas at Austin
%
% Version 1.0, June 2003
%
%
% TEXTCLUST is a free software and comes with ABSOLUTELY NO WARRANTY.
% You are free to use and redistribute it except for commercial purpose.
%
%
% Algorithms implemented:
%	k-means, stochastic k-means, EM, balanced k-means, and
%	deterministic annealing
%
% Models implemented:
%	Bernoulli, multinomial, and von Mises-Fisher (simplified)
%
% Bernoulli model-based clustering
%	kberns - Bernoulli-based k-means 
%	skberns - Bernoulli-based stochastic k-means
%	mixberns - Bernoulli-based EM
%	bkberns - balanced Bernoulli-based k-means
%	daberns - Bernoulli-based deterministic annealing
%
% Multinomial model-based clustering
%	kmns - multinomial-based k-means
%	skmns - multinomial-based stochastic k-means
%	mixmns - multinomial-based EM
%	bkmns - balanced multinomial-based k-means
%	damns - multinomial-based deterministic annealing
%
% von Mises-Fisher model-based clustering
%	kvmfs - vMF-based k-means
%	skvmfs - vMF-based stochastic k-means
%	mixvmfs - vMF-based soft clustering
%	bkvmfs - balanced vMF-based k-means
%	davmfs - vMF-based deterministic annealing
%
% Utility functions
%	cm - computes confusion matrix
%	mi - computes normalized mutual information
%	logidf - transforms dtm matrix using log(IDF) weighting
%	unitnorm - normalizes each row or column of a matrix into unit length
%	entro - calculates entropy (base-2) of a non-negative vector
%	entroa - average entropy of a matrix (average over rows)
%	puritya - average purity of a confusion matrix (average over rows)
%	bpart - balanced partition of a log-likelihood matrix
%
% Demo, example uses
%	test - example test code (needs enclosed data tr11.mat)
%	compare - example comparison code (needs enclosed data tr11.mat)

