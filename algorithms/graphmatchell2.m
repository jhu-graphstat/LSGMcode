function corr = graphmatchell2(A,B,s)
% Graph match non-seeds to non-seeds only

corr = seedgraphmatchell2(A(s+1:end,s+1:end), B(s+1:end,s+1:end), 0);
corr = [1:s corr+s];
