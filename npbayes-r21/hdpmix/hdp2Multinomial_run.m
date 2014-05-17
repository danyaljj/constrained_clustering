function [hdp,sample,lik,predlik] = hdp2Multinomial_run(hh,alphaa,alphab,...
    numclass,trainss,testss,...
    trainnumburnin,trainnumsample,trainnumspace,...
    testnumburnin,testnumsample,trainconparam,testconparam,fid);

func = hdpMultinomial_func;

hdp = hdp_init(func,0,1,hh,alphaa,alphab); 
[hdp trainindex] = hdp_adddp(hdp,length(trainss),1,2);
[hdp testindex]  = hdp_adddp(hdp,length(testss),1,2);


hdp = hdp_setdata(hdp,trainindex,trainss);
hdp = hdp_setdata(hdp,testindex,testss);

hdp = dp_activate(hdp,[1 trainindex],numclass);

[sample hdp lik] = hdp_posterior(hdp,trainnumburnin,trainnumsample,...
                   trainnumspace,trainconparam,1,0,fid);
hdp = dp_freeze(hdp,trainindex); 

predlik = hdp_predict(hdp,sample,trainindex,testindex,...
          testnumburnin,testnumsample,testconparam,0,fid);
