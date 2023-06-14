function [popn, stats, inputs] = runPopn_WR(statsmat,unitCIandVSoutputs,unitoutputs,allSpikeTimes88and91,options)
% Runs a population classifier for a specific set of units. 


%options.sampleMode = 'rand_replace';     % 'best','worst'
% switch(options.sampleMode)                  
%     case {'best'}
%         % For best we pick the ones with the highest mean dprimes of the
%         % the 1st 7 frequencies. 
%         % N.B. This will only run if there is a "standardised set" of
%         % classifier results in statsmat. 
%         meandprimes = cellfun(@(x) mean(x(1:min(length(x),7))),statsmat.dprimes_standardset);
%         [~, options.sampleOrder] = sort(meandprimes,'descend');
%         
%         if options.sampleIter > 1
%             warning('sampleIter>1 will just repeat processing needlessly');
%         end;
%     otherwise
% end;
% 

% Initialise outputs to empty.
popn = []; 
stats = [];
inputs = [];

% Options can contain options to pass on to the the "combineUnits"
% command. 
if isfield(options,'combineoptions')
    combineoptions = options.combineoptions;
else
    combineoptions = {};
end;


inputs.fullDataMUIndexes = [statsmat.fullDataUnitIndex{:}]; 
noofUnits = length(inputs.fullDataMUIndexes);

options.sampleP = nan; % How many to select as a percentage.
if noofUnits>1
    options.sampleN_Margin = 5;
else
    options.sampleN_Margin = 0;
end;
if ~isfield(options,'minNumFreqs')
    options.minNumFreqs = 3;
end;

% Set up the classifier. 
classifier.store_distances = 0;
classifier.iterations = 500;
classifier.sweeppermute = false;
classifier.sweeppermutations = 1;

% Get those individual units.
[inputs.MUspktrainsets inputs.MUunitinfos] =  ...
    findSpikeTimes(inputs.fullDataMUIndexes,unitCIandVSoutputs,allSpikeTimes88and91);

for iter = 1:options.sampleIter
    fprintf('population sampling iteration: %d\n',iter);
    
    % Loop around until you have an MU smple with a minimum number of
    % frequencies
    sampleGen_iter = 1;  % Counter to make sure we don't look forever.
    while true 
        % Random sub-sampling is done each iteration 
        switch(options.sampleMode)    
            case 'rand_replace'
               % N.B. We increase the sample slightly to account for the likelihood
               % of having to exclude some. 
               unitSample = randi(noofUnits,[options.sampleN+options.sampleN_Margin 1]);
               
            case {'best'}
                % Very simple sampling - picks the best ones with no random
                % element. 
                if options.sampleN+options.sampleN_Margin<= length(options.sampleOrder)
                    unitSample = options.sampleOrder([1:options.sampleN+options.sampleN_Margin]);
                else
                    unitSample = options.sampleOrder;
                    warning('Not enough to run with requested margin');
                end;
                
            case {'ordered'}
                % Very simple sampling - picks the best ones with no random
                % element. 
                if options.sampleN+options.sampleN_Margin<= length(options.sampleOrder)
                    unitSample = options.sampleOrder([1:options.sampleN+options.sampleN_Margin]);
                else
                    unitSample = options.sampleOrder;
                    warning('Not enough to run with requested margin');
                end;

            case {'one'}
                unitSample = 1;

        end;   

        % Put all the spike trains together - N.B. This can result in a very small
        % number of modulation frequencies 
        [MUspktrainset MUunitinfo{iter}] = ...
            combineUnits(inputs.MUspktrainsets(unitSample),inputs.MUunitinfos(unitSample), ...
                    'userationalised','ignorelevel','ignoreunits',combineoptions{:});

        if length( MUunitinfo{iter}.usedamfreqs ) >= options.minNumFreqs ...
               && length(MUspktrainset)>=options.sampleN
            break;
        end;
        if sampleGen_iter>10
            warning('Cannot generate requested sample of neurons');        
            return;
        end;
        sampleGen_iter = sampleGen_iter + 1;
    end;
      
    % Now we reduce the population again to the size we want.
    switch(options.sampleMode)   
        % This option takes only the last neuron in the set (so the worst
        % performer of the n best).
        % Fairly sure this is dead code.
%         case 'nth_best_only'
%             unitSample = unitSample(options.sampleN);
%             MUspktrainset = MUspktrainset(options.sampleN);
%             fullDataMUIndexes = popn.fullDataMUIndexes(unitSample);
%             MUunitinfo{iter} = MUunitinfo{iter}(options.sampleN);
        otherwise
            unitSample = unitSample(1:options.sampleN);
            MUspktrainset = MUspktrainset(1:options.sampleN);
            fullDataMUIndexes = inputs.fullDataMUIndexes(unitSample);
            
            % Also shrink MUunitinfo
            MUunitinfo{iter}.unitid = MUunitinfo{iter}.unitid(1:options.sampleN);
            MUunitinfo{iter}.carrFreq = MUunitinfo{iter}.carrFreq(1:options.sampleN);
            MUunitinfo{iter}.cf = MUunitinfo{iter}.cf(1:options.sampleN);
            MUunitinfo{iter}.cf_thr = MUunitinfo{iter}.cf_thr(1:options.sampleN);
            MUunitinfo{iter}.TypeNum = MUunitinfo{iter}.TypeNum(1:options.sampleN);
            MUunitinfo{iter}.TypeName = MUunitinfo{iter}.TypeName(1:options.sampleN);
            MUunitinfo{iter}.wholeListInds = MUunitinfo{iter}.wholeListInds(1:options.sampleN);
            MUunitinfo{iter}.scholesWholeListFirstEntry = MUunitinfo{iter}.scholesWholeListFirstEntry(1:options.sampleN);
            MUunitinfo{iter}.fullDataMUIndexes = fullDataMUIndexes;
    end;
    
    % Select the right time constants from the individual data.
    classifier.tau = arrayfun( @(x) x.bestClassifier.tau, ...
        unitoutputs(fullDataMUIndexes));

    % Run the classifier
    popn.WRi{iter} = WR_MUClassifier(MUspktrainset, ...
        ones(size(MUspktrainset{1},1),1)*MUunitinfo{iter}.usedamfreqs, ...
        classifier);

    % Store whichi neurons were used in this iteration. 
    popn.WRi{iter}.fullDataMUIndexes = fullDataMUIndexes;
    popn.WRi{iter}.sampleOrder = unitSample;
end;    

% Store the unit info for each separate iteration.
popn.MUunitinfo = MUunitinfo;
%popn  = cellfun(@(x) rmfield(x,{'MUspktrainsets','MUunitinfos'}),popn,'uni',false);

% Store an overall summary of the runs. 
stats.WR = popn.WRi{1};

% These may vary so remove. 
stats.WR.classifier.tau = nan;
stats.WR.classifier.alpha = nan;

% Derive the full set of frequencies tested. 
allamflists = cellfun(@(x) x.stimulusconditions',popn.WRi,'Uni',false);
amflist = unique([allamflists{:}]);
stats.WR.stimulusconditions = amflist;

WRhitrate = nan(length(amflist),options.sampleIter);
WRfalsealarmrate = nan(length(amflist),options.sampleIter);
WRdprime = nan(length(amflist),options.sampleIter);
WRsweeps = nan(length(amflist),options.sampleIter);
for iter = 1:options.sampleIter    
    [ism, inds] = ismember(amflist,popn.WRi{iter}.stimulusconditions);
    inds = inds(ism);
    WRhitrate(inds,iter) = popn.WRi{iter}.hitrate;
    WRfalsealarmrate(inds,iter) = popn.WRi{iter}.falsealarmrate;
    WRdprime(inds,iter) = popn.WRi{iter}.dprime;
    WRsweeps(inds,iter) = popn.WRi{iter}.sweepspercondition;
end;

% dprime statistics.
stats.WR = rmfield(stats.WR,{'dprime'});
stats.WR.dprime.mean = nanmean(WRdprime,2);
stats.WR.dprime.sd = nanstd(WRdprime,0,2);
stats.WR.dprime.n = nansum(~isnan(WRdprime),2);

% hitrate statistics.
stats.WR = rmfield(stats.WR,{'hitrate'});
stats.WR.hitrate.mean = nanmean(WRhitrate,2);
stats.WR.hitrate.sd = nanstd(WRhitrate,0,2);
stats.WR.hitrate.n = nansum(~isnan(WRhitrate),2);

% false alarm rate statistics.
stats.WR = rmfield(stats.WR,{'falsealarmrate'});
stats.WR.falsealarmrate.mean = nanmean(WRfalsealarmrate,2);
stats.WR.falsealarmrate.sd = nanstd(WRfalsealarmrate,0,2);
stats.WR.falsealarmrate.n = nansum(~isnan(WRfalsealarmrate),2);

% sweep statistics.
stats.WR = rmfield(stats.WR,{'sweepspercondition'});
stats.WR.sweepspercondition.mean = nanmean(WRsweeps,2);
stats.WR.sweepspercondition.sd = nanstd(WRsweeps,0,2);
stats.WR.sweepspercondition.n = nansum(~isnan(WRsweeps),2);

