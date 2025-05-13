function bestwroutputs = getWROutputsWithTau(unitoutputs,tau)
% Selects out the specified classifier (tau) for each set.

bestwroutputs = unitoutputs;
for i =1:length(unitoutputs)
    % Get all the taus for this dataset (all the same but hey).  
    if ~isempty(unitoutputs(i).wroutput)
        taus = arrayfun(@(x) x.classifier.tau,unitoutputs(i).wroutput);
        ind = find(taus==tau);        
        bestwroutputs(i).wroutput = unitoutputs(i).wroutput(ind);
    end;
end;
    
   
