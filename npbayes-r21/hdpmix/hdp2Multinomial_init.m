function hdp = hdp2Multinomial_init(hh,alphaa,alphab,...
    numclass,trainss,testss);

func = hdpMultinomial_func;

hdp = hdp_init(func,0,1,hh,alphaa,alphab); 
[hdp trainindex] = hdp_adddp(hdp,length(trainss),1,2);
[hdp testindex]  = hdp_adddp(hdp,length(testss),1,2);

hdp = hdp_setdata(hdp,trainindex,trainss);
hdp = hdp_setdata(hdp,testindex,testss);

hdp = dp_activate(hdp,[1 trainindex],numclass);

