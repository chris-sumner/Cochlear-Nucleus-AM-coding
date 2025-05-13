function [tmp, p0] = prepTheoreticalMaxFit_logistic(stats_stacked_4models,gains,x0,initModFn)
% prepTheoreticalMaxFit pre-compute all that we can before fitting. 

% Function to invert the logistic function and compute the underlying metric from
% pecentage correct. It is not as good as fitting, because it cannot
% account fully for the compression. 
c_metric_fn = @(p,p0,x0)  -log( (1 -p0)./(p-p0)  - 1 ) + x0;
 
% Prepare the data structure for fitting.
tmp.unitlist = unique([stats_stacked_4models.unitid' stats_stacked_4models.modLvl_rationalised'  stats_stacked_4models.depthMod'], 'rows' ); % Loop through each dataset.
tmp.gains = gains;
tmp.x0 = x0;
tmp.initModFn = initModFn;  % Resonable to initialise to the mean normalised modulation transfer function.
flist = unique(stats_stacked_4models.usedamfreqs_dprime_rationalised);
tmp.flist = flist;
stats_stacked_4models.p0 = nan(size(stats_stacked_4models.dprimes));

% We pre-compute as much as possible - assumble info for each record.
for ui = 1:length(tmp.unitlist)
    % The indexes in stats_stacked_4models where this record is.
    tmp.inds{ui} = find(stats_stacked_4models.unitid == tmp.unitlist(ui,1) & ...
                    stats_stacked_4models.modLvl_rationalised == tmp.unitlist(ui,2) & ...
                    stats_stacked_4models.depthMod == tmp.unitlist(ui,3) );
    % The outcome measure as a function of modulation frequency for this record.
    tmp.measure{ui} = stats_stacked_4models.dprimes(tmp.inds{ui});
    % Which AM frequencies this record has 
    tmp.amfreqs{ui} = stats_stacked_4models.usedamfreqs_dprime_rationalised(tmp.inds{ui});
    % Indexes speed the computation for each record. 
    tmp.amfreq_inds{ui} =( tmp.amfreqs{ui}+50)/100;
    % Chance is 1/no of am freqs for this record.
    tmp.pchance{ui} = 1/length( tmp.measure{ui} );
    % p0 is a correction factor so the logistic function floor is at chance.
    tmp.p0{ui} = (tmp.pchance{ui}*(1+exp(tmp.x0))-1)/exp(tmp.x0);

    % Add p0 to stats_stacked_4models for fitting.
    stats_stacked_4models.p0(tmp.inds{ui}) =   tmp.p0{ui};

    % We can also compute the metric for a given dataset.  
    % Fitting to this is not a good idea because we are trying to predict the data, no the metric.  
    tmp.c_metric{ui} =  c_metric_fn( tmp.measure{ui}, tmp.p0{ui}, tmp.x0 );
end;

p0 = stats_stacked_4models.p0;