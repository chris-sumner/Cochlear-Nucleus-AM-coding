function [tmp, ClassIndex, DataIndex] = prepTheoreticalMaxFit_softmax(stats_stacked_4models,gains,initModFn)
% prepTheoreticalMaxFit pre-compute all that we can before fitting. 


% Prepare the data structure for fitting.
tmp.unitlist = unique([stats_stacked_4models.unitid' stats_stacked_4models.modLvl_rationalised' ...
    stats_stacked_4models.depthMod' stats_stacked_4models.fullDataUnitIndex'], 'rows' ); % Loop through each dataset.
tmp.gains = gains;
tmp.initModFn = initModFn;  % Resonable to initialise to the mean normalised modulation transfer function.
flist = unique(stats_stacked_4models.usedamfreqs_dprime_rationalised);
tmp.flist = flist;

% In order to calculate hitrates from the fitted softmax models, we need
% the B parameter vector. This is not stored in the main datastructure and is
% complicated to do because we need ALL the B values from teh data set.
% This means we can just look them up (on the fly during fitting).
stats_stacked_4models.ClassIndex = nan(size(stats_stacked_4models.dprimes));
stats_stacked_4models.DataIndex = nan(size(stats_stacked_4models.dprimes));

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

    % Add DataIndex to stats_stacked_4models for fitting.
    tmp.DataIndex{ui} = stats_stacked_4models.fullDataUnitIndex(tmp.inds{ui});
    %stats_stacked_4models.DataIndex(tmp.inds{ui}) =   tmp.DataIndex{ui};

    % Add "ClassIndex" which is just the index of the frequency within the
    % set used for the MTF.
    % Any not included are left as nans. 
    stats_stacked_4models.ClassIndex(tmp.inds{ui}) = 1:length(tmp.amfreq_inds{ui});

end;

% Passing these out seaprately in the right (long) format for adding back into the stacked data. 
ClassIndex = stats_stacked_4models.ClassIndex;  % should not need this. 
tmp.ClassIndex_stacked = stats_stacked_4models.ClassIndex;