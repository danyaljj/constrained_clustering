function [A, b, c] = formulateSDP(S, D, bb)
[F0, FI, c] = localformulateSDP(S, D, bb);
[A, b, c] = sdpToSeDuMi(F0, FI, c);




