function bestwroutputs = getBestWROutputs(unitoutputs)
% Selects out the best performing classifier (tau) for each set.
inds = arrayfun(@(x) isfield(x.bestClassifier,'index'),unitoutputs);
bestwroutputs = unitoutputs;
for i =1:length(unitoutputs)
    if inds(i)~=0
        bestwroutputs(i).wroutput = unitoutputs(i).wroutput(unitoutputs(i).bestClassifier.index);
    end;
end;
    