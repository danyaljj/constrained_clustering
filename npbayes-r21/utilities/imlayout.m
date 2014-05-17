function handle = imlayout(w,shape,s,varargin);
if nargin<2,
 shape = [size(w,1) size(w,2) size(w,3) size(w,4)];
 s = [min(w(:)) max(w(:))];
elseif nargin<3,
 w = reshape(w,shape);
 s = [min(w(:)) max(w(:))];
else
 w = reshape(w,shape);
 w(:) = min(s(2),max(s(1),w(:)));
end

if length(shape)<4,
 shape(end+1:4) = 1;
end
w = (w-s(1))/(s(2)-s(1));

v = shape(1)*shape(3);
h = shape(2)*shape(4);
w = reshape(permute(w,[1 3 2 4]),[v h]);
w = w(:,:,ones(1,3));
handle = image(w);
hold on; axis equal; axis off;
for i=0:shape(4),
 plot([i*shape(2)+.5 i*shape(2)+.5],[+.5 v+.5],varargin{:});
end;
for j=0:shape(3),
 plot([+.5 h+.5],[j*shape(1)+.5 j*shape(1)+.5],varargin{:});
end;
hold off
