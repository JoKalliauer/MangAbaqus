function [idxs,xPoints,yHalfPoints,yDiffPoints] = markWende(lambdaplot,y4)
  lambdaplotHalf=(lambdaplot(1:end-1)+lambdaplot(2:end))/2;
  y4Half=(y4(1:end-1)+y4(2:end))/2;
  y4Diff=y4(2:end)-y4(1:end-1);
  idxsMin=islocalmin(y4Diff,'FlatSelection', 'all');
  idxsMax=islocalmax(y4Diff,'FlatSelection', 'all');
  idxs=logical(idxsMin+idxsMax);
  xPoints=lambdaplotHalf(idxs);
  yHalfPoints=y4Half(idxs);
  yDiffPoints=y4Diff(idxs);
  lambdaplotHalf(idxs)
end

