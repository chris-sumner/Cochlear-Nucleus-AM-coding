function err = mtfErrFn(f_fn,tmp)

% ---- Boils down to this for each iteration ---. 

%f_fn = meanfModFn - min(meanfModFn);  % This is the weights for the model which are optimising. 

logitfn = @(g) 1./(1 + exp(-(g*f_fn-tmp.x0))); 
pcfn = @(p0,lt) p0 +(1-p0).*lt;
logittable = logitfn(tmp.gains');  % Generic table of functions for all neurons 

% Go through all neurons and comute the error. 
for ui = 1:length(tmp.unitlist)   
    tmp.pctable = pcfn(tmp.p0{ui},logittable(:,tmp.amfreq_inds{ui}));        
    tmp.difference = nanmean((ones(101,1)*tmp.measure{ui} - tmp.pctable).^2,2);
    tmp.best_ind(ui) =find(tmp.difference == min(tmp.difference) );
    tmp.fitted_err(ui) = tmp.difference(tmp.best_ind(ui));
end;
 
err = sum( tmp.fitted_err );

