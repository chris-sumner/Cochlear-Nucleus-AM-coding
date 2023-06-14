%% ------------------------------------------------------------------------
%                    Calculating VS
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------
%    set up paths, load raw data.  

% Add the paths required.
addpath functions;

% Load up spike times. 
% Split across 4 files to fit on github.
allSpikeTimes88and91 = loadSpikeTimes;

% Loads unit info and other supporting data.
load('rawdata\unitList');

%% -----------------------------------------------------------------------
%         This calulates VS test for every condition, every unit. 

win100 = [20 100];          % N.B. Any tones >100ms are only analysed to 100ms. 
win50  = [20 50];
numPHbins = 52;
phbin = 0:1/numPHbins:1; % see Malone et al. 2007

VS = cell(size(allSpikeTimes88and91));

completedvector = cell(size(allSpikeTimes88and91));

for ii = 1:length(allSpikeTimes88and91)
    theseSpikeTimes = allSpikeTimes88and91{ii};
    duration = wholeList171212(ii,14);

    if duration>=100
        win2use = win100;
    elseif duration == 50
        win2use = win50;
    else
        win2use=win100;
    end;
    
    if ~isempty([theseSpikeTimes{:}]) 
        
        % Select out the right bit of the spike train.
        spktrainset = cellfun(@(x)  x(x>=win2use(1) & x<=win2use(2)),theseSpikeTimes,'UniformOutput',false);

        % Some info about the spike train.     
        nsweeps = wholeList171212(ii,19);
        thisModFreq = wholeList171212(ii,13);
        per = 1000 / thisModFreq;
        
        % ------------- vector strength ----------------
                
        % calculate spike-phase as a fraction of the period.
        st = [spktrainset{:}];
        p_st = (st - floor(st/per)*per)/per; % spike times in periods
        l_st = sum(~isnan(p_st));            % number of spikes without nans.

        VS{ii}.VSvalue = sqrt(sum(cos(2*pi*p_st))^2+sum(sin(2*pi*p_st))^2)/l_st;
        VS{ii}.Rayleigh = 2*l_st*(VS{ii}.VSvalue^2);
        VS{ii}.PH = hist(p_st,phbin);
        VS{ii}.PHbins = phbin;
        VS{ii}.spikes_per_sweep = length(st)/nsweeps;
        
        completedvector{ii}.value = 1;
        fprintf('.');
    end;
        
end

% save the calculated statistics. 
save datafiles\VS VS;



%% ------------------------------------------------------------------------
%              Reshape this into individual sets for each MTF.

clear; 

% Add the paths required.
addpath functions;

% Loads unit info and other supporting data.
load('rawdata\unitList');

% Load up the CI and VSs . 
load datafiles\VS VS;

% There is one decision to be made here which is the range of
% frequencies submitted to the classifier - if not all of them. This is
% because this affects the classifier results at other frequencies.

% indices of lists to keep: modFreq<=2000 - this is most of them anyway.
keepConds = wholeList171212(:,find(strcmp(EXPLOGLIST,'modFreq')))<=2000; 
% Shorten the data to exclude MTFs we will not analyse. 
wholeList171212 = wholeList171212(keepConds,:);
allAMtypes88and91 = allAMtypes88and91(keepConds);
% reduce the size of VS accordingly.
VS = VS(keepConds);

% dividing the data up into individual datasets in units.
[newdatainds dividingdata] = divideData(wholeList171212,EXPLOGLIST,allAMtypes88and91);

% Which data to process. 
datalist = 1:length(newdatainds);
minspkpersweep = 2;     % Minimum spikes per sweep: 2 = 25 spikes per sec for an 80ms window.
                        % Below this VS is unreliable.

% Loop to restructure the data into MTF dataset groups.
for thisdatastartind = datalist
    
    % The indexes in the wholelist format where the data should exist.
    thisDataWholeListInds = [];
    if thisdatastartind<length(newdatainds)
        endind = newdatainds(thisdatastartind+1)-1;
    else
        endind = size((newdatainds),1);
    end
    thisDataWholeListInds = [newdatainds(thisdatastartind):endind];

    % Which of these conditions have enough spikes in them.
    conditionswithspikes = [];
    if ~isempty( thisDataWholeListInds )
        % Take this from the VS analysis
        vstmp = VS(thisDataWholeListInds);
        inds = find(cellfun(@(x) ~isempty(x),vstmp));
        spikespercondition = zeros(size(vstmp));
        spikespercondition(inds) = cellfun(@(x) x.spikes_per_sweep, vstmp(inds) );
        conditionswithspikes = find(spikespercondition>=minspkpersweep);
    end;        

    if ~isempty(thisDataWholeListInds) && ~isempty(conditionswithspikes)
        % This stores the actual VS analysis
        unitVSoutputs(thisdatastartind).VSstats = ...
            [VS{ thisDataWholeListInds( conditionswithspikes )   }];
    else
        unitVSoutputs(thisdatastartind).VSstats = [];
    end;
                
    % make the unitinfo.
    unitVSoutputs(thisdatastartind).unitinfo = ...
        constructUnitInfo(wholeList171212,EXPLOGLIST,typeNames,thisdatastartind, ...
         thisDataWholeListInds,newdatainds,conditionswithspikes);
        
    fprintf('.');     
end;

% Add this to the file.
save datafiles\VS unitVSoutputs -append;