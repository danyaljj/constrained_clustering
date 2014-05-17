
cd hdpmix
try
  mex hdpMultinomial_iterate.c -lm % add debugging options? -DDEBUG??
  mex hdpMultinomial_predict.c -lm 
catch
end
cd ..
