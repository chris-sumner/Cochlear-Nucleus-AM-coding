function tmpout =  postTheoreticalMaxFit(tmp)



% We make a normalised version. Or renormalised, because the fit did not
% constrain the final values. 
% This will increase the fitted gains in proportion, but otherwise make no difference. 
tmp.normalised_fitted_fn = tmp.fitted_fn/max(tmp.fitted_fn);

% Define the functions.
logitfn = @(g) 1./(1 + exp(-(g*tmp.normalised_fitted_fn-tmp.x0)));         % Calcuates a logit function for a specified gain. 
pcfn = @(p0,lt) p0 +(1-p0).*lt;

% We don't know what the gains were for each record so
% recompute the fits for the indidual records.
tmp.logittable = logitfn(max(tmp.fitted_fn)*tmp.gains');                       % Generic table of functions - precomputed for all gains.

% Go through all neurons and compute the gains error. 
for ui = 1:length(tmp.unitlist)   
    % The table for looking up the correct gain.
    tmp.pctable = pcfn(tmp.p0{ui},tmp.logittable(:,tmp.amfreq_inds{ui}));   
    % Differences between data and each row of lookup table.
    tmp.difference = nanmean((ones(101,1)*tmp.measure{ui} - tmp.pctable).^2,2);
    % Find the gain with the minimum error.
    tmp.best_ind(ui) = find(tmp.difference == min(tmp.difference) );
    tmp.fitted_err(ui) = tmp.difference(tmp.best_ind(ui));
    tmp.fitted_normalised_gains(ui) = max(tmp.fitted_fn)*tmp.gains(tmp.best_ind(ui));
    tmp.fitted_c_metric{ui} = tmp.fitted_normalised_gains(ui).*tmp.normalised_fitted_fn( tmp.amfreq_inds{ui} );
    tmp.fitted_measures{ui} = tmp.pctable(tmp.best_ind(ui),:);
end;

tmp.measureV = [tmp.measure{:}];
tmp.predmeasureV = [tmp.fitted_measures{:}];
tmp.nonnans = ~isnan(tmp.measureV) &  ~isnan(tmp.predmeasureV);
tmp.totalCorrCoef = corrcoef(tmp.measureV(tmp.nonnans),tmp.predmeasureV(tmp.nonnans));
tmp.overallR = tmp.totalCorrCoef(2);

tmpout = tmp;