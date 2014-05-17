
cd utilities
try
  mex randgamma.c  -lm % add debugging options? -DDEBUG??
  mex randnumtable.c -lm
catch
end
cd ..
