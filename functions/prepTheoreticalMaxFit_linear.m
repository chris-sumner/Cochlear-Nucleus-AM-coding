function [tmp] = prepTheoreticalMaxFit_linear(stats_stacked_4models,gains,initModFn)
% prepTheoreticalMaxFit_linear pre-compute all that we can before fitting. 

% Function which does nothing but keeps the code similar to the logistic fit.
c_metric_fn = @(p) p;
 
% Prepare the data structure for fitting.
tmp.unitlist = unique([stats_stacked_4models.unitid' stats_stacked_4models.modLvl_rationalised' ...
    stats_stacked_4models.depthMod' stats_stacked_4models.fullDataUnitIndex'], 'rows' ); % Loop through each dataset.
tmp.gains = gains;
tmp.initModFn = initModFn;  % Resonable to initialise to the mean normalised modulation transfer function.
flist = unique(stats_stacked_4models.usedamfreqs_dprime_rationalised);
tmp.flist = flist;

% We pre-compute as much as possible - assumble info for each record.
for ui = 1:length(tmp.unitlist)
    % The indexes in stats_stacked_4models where this record is.
    tmp.inds{ui} = find(stats_stacked_4models.unitid == tmp.unitlist(ui,1) & ...
                stats_stacked_4models.modLvl_rationalised == tmp.unitlist(ui,2) & ...
                stats_stacked_4models.depthMod == tmp.unitlist(ui,3) & ...
                stats_stacked_4models.fullDataUnitIndex == tmp.unitlist(ui,4) );
    
    % We need to check if the measure has any inf or nans which would screw
    % up the calclulation. We remove those frequencies. 
    thismeasure = stats_stacked_4models.dprimes(tmp.inds{ui});
    goodones = (~isnan(thismeasure ) & ~isinf( thismeasure));  
    tmp.inds{ui} = tmp.inds{ui}(goodones);
                        
    % The outcome measure as a function of modulation frequency for this record.
    tmp.measure{ui} = stats_stacked_4models.dprimes(tmp.inds{ui});
    
    % Which AM frequencies this record has 
    tmp.amfreqs{ui} = stats_stacked_4models.usedamfreqs_dprime_rationalised(tmp.inds{ui});
    % Indexes speed the computation for each record. 
    tmp.amfreq_inds{ui} =( tmp.amfreqs{ui}+50)/100;
    
    % We can also compute the metric for a given dataset.  
    % Does nothing for a linear fitting anyway. 
    tmp.c_metric{ui} =  c_metric_fn( tmp.measure{ui});
end;
