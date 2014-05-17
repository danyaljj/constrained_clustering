function [hdp, ll] = hdp_iterate(hdp,numiter,doconparam,dolik,dodebug);

[hdp ll] = feval(hdp.func.iterate,hdp,numiter,doconparam,dolik,dodebug);
