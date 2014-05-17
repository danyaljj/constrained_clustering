function showtopics(qq,words,fid);

if nargin==2
  fid = 1;
  stop = 1;
elseif fid==0
  fid = 1;
  stop = 1;
else
  stop = 0;
end
for cc = 1:size(qq,2)
  [mm ii] = sort(-qq(:,cc));
  for jj = 1:10
    fprintf(fid,'%s ',words{ii(jj)});
    %fprintf(1,'%s ',words{ii(jj)})
  end
  fprintf(fid,'\n');
  if stop==1
    pause;
  end
end
