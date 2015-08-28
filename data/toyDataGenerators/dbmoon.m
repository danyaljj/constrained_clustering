function data=dbmoon(N,d,r,w)
% Usage: data=dbmoon(N,d,r,w)
% doublemoon.m - genereate the double moon data set in Haykin's book titled
% "neural networks and learning machine" third edition 2009 Pearson
% Figure 1.8 pp. 61
% The data set contains two regions A and B representing 2 classes
% each region is a half ring with radius r = 10, width = 6, one is upper
% half and the other is lower half
% d: distance between the two regions
% will generate region A centered at (0, 0) and region B is a mirror image
% of region A (w.r.t. x axis) with a (r, d) shift of origin
% N: # of samples each class, default = 1000
% d: seperation of two class, negative value means overlapping (default=1)
% r: radius (default=10), w: width of ring (default=6)
% 
% (C) 2010 by Yu Hen Hu
% Created: Sept. 3, 2010

% clear all; close all;
if nargin<4, w=6;
elseif nargin<3, r=10;
elseif nargin<2, d=1;
elseif nargin < 1, N=1000;
end

% generate region A:
% first generate a uniformly random distributed data points from (-r-w/2, 0)
% to (r+w/2, r+w/2)
N1=10*N;  % generate more points and select those meet criteria
w2=w/2; 
done=0; data=[]; tmp1=[];
while ~done, 
    tmp=[2*(r+w2)*(rand(N1,1)-0.5) (r+w2)*rand(N1,1)];
    % 3rd column of tmp is the magnitude of each data point
    tmp(:,3)=sqrt(tmp(:,1).*tmp(:,1)+tmp(:,2).*tmp(:,2)); 
    idx=find([tmp(:,3)>r-w2] & [tmp(:,3)<r+w2]);
    tmp1=[tmp1;tmp(idx,1:2)];
    if length(idx)>= N, 
        done=1;
    end
    % if not enough data point, generate more and test
end
% region A data and class label 0
% region B data is region A data flip y coordinate - d, and x coordinate +r
data=[tmp1(1:N,:) zeros(N,1);
    [tmp1(1:N,1)+r -tmp1(1:N,2)-d ones(N,1)]];

% plot(data(1:N,1),data(1:N,2),'.r',data(N+1:end,1),data(N+1:end,2),'.b');
% title(['Fig. 1.8 Double moon data set, d = ' num2str(d)]),
% axis([-r-w2 2*r+w2 -r-w2-d r+w2])

%save dbmoon N r w d data;
