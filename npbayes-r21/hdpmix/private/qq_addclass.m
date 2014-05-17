function hdp = qq_addclass(hdp,newclass);

HELDOUT = 0;
ACTIVE = 2;

oldclass = hdp.base.numclass;
numclass = oldclass + newclass;
hdp.base.numclass = numclass;
hdp.base.classqq(:,oldclass+1:numclass+1) = ...
    hdp.base.classqq(:,oldclass+ones(1,1+newclass));

for jj = 1:hdp.numdp
  if hdp.dpstate(jj) ~= HELDOUT
    hdp.dp{jj}.classnd(numclass+1) = 0;
    hdp.dp{jj}.classnt(numclass+1) = 0;
  end
  if hdp.dpstate(jj) == ACTIVE
    hdp.dp{jj}.beta(oldclass+1:numclass+1) = hdp.dp{jj}.beta(oldclass+1)*...
        randstick(hdp.dp{jj}.alpha,newclass+1);
  end 
end
