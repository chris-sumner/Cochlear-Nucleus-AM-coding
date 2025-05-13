function [spktrainsets, unitinfoout] = findSpikeTimes(unitList,unitCIandVS,allSpikeTimes88and91)
% findSpikeTimes   retreives spike times for a list of units. 
%
% This function works from the analysis of CI/VS, to produce the original
% spike times. This is so that we can go back to the data for additional
% analysis, such as to look at population coding. 
%
% [spktrainset, thisunitinfo] = findSpikeTimes(i,unitCIandVSoutputs,allSpikeTimes88and91)
% 
% i is the unit index in unitCIandVSoutputs.
% The 2nd and 3rd parameters are not optional event though they are relatively fixed. 

load('rawdata\unitList');
colfn =  @(s) find(strcmp(EXPLOGLIST,s));

% Can only use "wholeListInds" list if you first remove any frequencies>2000
keepConds = find(wholeList171212(:,find(strcmp(EXPLOGLIST,'modFreq')))<=2000); 
% Shorten the data to exclude MTFs we will not analyse. 
wholeList171212 = wholeList171212(keepConds,:);
allAMtypes88and91 = allAMtypes88and91(keepConds);
allSpikeTimes88and91 = allSpikeTimes88and91(keepConds);


% -------------------- Check the conditions match up -------------------

% We check the stimulus conditions in unitinfo match up with those in the
% "wholeList...." array which is row matched to the spike times cell
% array.

for ui=1:length(unitList)

    % The index in that dataset for the unit.
    i = unitList(ui);
    
    % Local copy of the unitinfo for this data set for convenience.
    thisunitinfo = unitCIandVS(i).unitinfo;

    % --------------- retrieve the AM type ----------------
    datasetChk = 1;
    if ~isempty( thisunitinfo.wholeListInds )
        AMType = unique( allAMtypes88and91( thisunitinfo.wholeListInds ) );

         % Store: give warning if ambiguous.
         if length(AMType)>1
             error('AMType is ambiguous for this dataset: %d\n',i);
             AMType = 'BAD';
         end;

        % Double check the unit info against to the conditions in whole unit list.
        % Extremely unlikely to be wrong but...

        % The Scholes' wholeList for this dataset.
        dataSetScholesList = wholeList171212(thisunitinfo.wholeListInds,:); 

        % Check all these things agains the unitinfo.
        % N.B. Special treatment of carrier frequency in case of nans (ie.
        % noise carrier, cf and cf_thr can all be nan)
        datasetChk = prod( ...
                dataSetScholesList(:,colfn('Unit')) == rem(thisunitinfo.unitid,1000)  & ...
                dataSetScholesList(:,colfn('Expt')) == round(thisunitinfo.unitid/1000)  & ...
                ( dataSetScholesList(:,colfn('cf')) == thisunitinfo.cf | ...
                      isnan(dataSetScholesList(:,colfn('cf'))) & isnan(thisunitinfo.cf) ) & ...
                ( dataSetScholesList(:,colfn('cf_thr')) == thisunitinfo.cf_thr | ...
                      isnan(dataSetScholesList(:,colfn('cf_thr')))  & isnan( thisunitinfo.cf_thr ) ) &...
                ( dataSetScholesList(:,colfn('carrFreq')) == thisunitinfo.carrFreq | ...
                    isnan(dataSetScholesList(:,colfn('carrFreq'))) == isnan( thisunitinfo.carrFreq ) ) & ...
                dataSetScholesList(:,colfn('amToneDur')) == thisunitinfo.amToneDur & ...
                dataSetScholesList(:,colfn('modLevel')) == thisunitinfo.modLevel & ...
                dataSetScholesList(:,colfn('depthMod')) == thisunitinfo.depthMod ...
                );

        % If this fails then the AM type cannot be correct.    
        if ~datasetChk
            error('Failure of correspondence of unitinfo with wholeList171212 for dataset %d\n',i);
        end;
    end;     


    % ----------------------- Get the spikes -----------------------

    % Return all the spike times and the unitinfo.     
    spktrainset =[ allSpikeTimes88and91{ thisunitinfo.wholeListInds }];

    % In the orginal classifier we preproccess the spikes in several ways.
    % We replicate that here (not great programming). 

    % Cut down to within 20-100ms
    spktrainset = cellfun(@(x)  x(x>20 & x<=100),spktrainset,'UniformOutput',false);

    % Remove any conditions with no spikes in any sweep.
    spktrainset_chk = cellfun(@(x) ~isempty(x), spktrainset);
    conditionswithspikes = find(sum(spktrainset_chk,1));
    spktrainset = spktrainset(:,conditionswithspikes);
    amfreqset = thisunitinfo.allamfreqs(conditionswithspikes)';

    % And modify unit info to reflect that.
    % This means we have decided which are the "usedamfreqs" afresh (which
    % is good). 
    thisunitinfo.usedamfreqs = thisunitinfo.allamfreqs(conditionswithspikes)';
    thisunitinfo.conditionswithspikes = conditionswithspikes;
        
    spktrainsets{ui} = spktrainset;
    unitinfoout(ui) =  thisunitinfo;
end;




