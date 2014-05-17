function ni = qq_numitems(hdp,cc);

ni = feval(hdp.func.numitems,hdp.base.hh,hdp.base.classqq(:,cc));
