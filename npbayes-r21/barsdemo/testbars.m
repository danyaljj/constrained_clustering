% bars data.  starts with number of class = number of dimensions

% description of generated data set.
imsize       = 5;
noiselevel   = .01;
numbarpermix = [0 ones(1,3)]/3;
numgroup     = 50;
numdata      = 50;

hh = ones(imsize*imsize,1)/imsize;

echo on

% This script demonstrates the HDP working on the topic bars which first 
% appears in Blei et al (Nested CRP; NIPS 2004).  Training set consists 
% of 50 groups of data, each with 50 data points arranged into a number 
% of horizontal and vertical bars (5x5 image). Number of bars is between 
% 1 and 3 for each group. 

trainss      = genbars(imsize,noiselevel,numbarpermix,numgroup,numdata);
pause

% Test set consists of 25 groups, one consisting of just one data point 
% corresponding to one pixel.  This is simply to test that our method to 
% estimate predictive probabilities is accurate.  

testss       = num2cell(1:25);
pause

% Now we initialize and run the HDP on the training set, and estimate
% probabilities of test groups given training groups.  

[hdp sample lik predlik] = hdp2Multinomial_run(hh,[1 1],[1 1],15,...
    trainss,testss,1000,10,100,500,1000,15,1,1);
pause

% Plot of the likelihood of training data over time.

plot(lik);
pause

% A sample of the final estimated clusters can be visualized as images.

imlayout(-hdp.base.classqq,[5 5 1 hdp.base.numclass+1]);
pause

% The predictive probabilities are returned in predlik.  If the estimation
% procedure is accurate, the sum of the predictive probabilities should be
% close to 1:

sum(exp(meanlik(predlik)))
pause

% Finally we show the MCMC iterations running:

hdp = hdp2Multinomial_init(hh,[1 1],[1 1],15,trainss,testss);
for iter=1:100
  hdp = hdp_iterate(hdp,10,15,0,0);
  imlayout(-hdp.base.classqq,[5 5 1 hdp.base.numclass+1]);
  drawnow;
end
pause

% And done!
echo off

