function hdp = hdp_randbeta(hdp);

ACTIVE = 2;

basebeta = [zeros(1,hdp.base.numclass) 1];

for jj = 1:hdp.numdp
  if hdp.dpstate(jj) == ACTIVE
    pp = hdp.ppindex(jj);
    if pp == 0
      beta = basebeta;
    else
      beta = hdp.dp{pp}.beta;
    end
    hdp.dp{jj}.beta = randdir(hdp.dp{jj}.alpha*beta);
  end
end
