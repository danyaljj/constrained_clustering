function el = ellipse(xx,var,varargin);

[v,d] = eig(var);
d = sqrt(d);
theta = -.05 :.05: 2*pi;
xy = repmat(xx,1,length(theta))+...
     d(1,1)*v(:,1)*sin(theta)+d(2,2)*v(:,2)*cos(theta);
el = plot(xy(1,:),xy(2,:),varargin{:});

