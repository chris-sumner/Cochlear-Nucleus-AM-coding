function WRTraces(Group_S,Group_T,a)
% WRTraces

% Based on code by Robert Mill <rob.mill.uk@gmail.com>

N = length(Group_S);

% N.B. Does not to every combination - for display purposes only.
for ri = 1:N-1
    op{ri} = WRTrace(Group_S{ri+1},Group_T{ri},a,[20 100])
    diffSqMat(:,ri) = op{ri}.diffSqTrace;
end;

meanDiffSq = mean(diffSqMat,2);

figure; plot(meanDiffSq); hold on;




