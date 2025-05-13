%% ------------------------------------------------------------------------
%                    Calculating VS, mode-locking, CI .... etc
% -------------------------------------------------------------------------

% Takes more than two days to run all analyses execpt SAC peak-picking,  
% which takes weeks. 

% ------------------------------------------------------------------------
%    set up paths, load raw data.  

clear;

% Add the paths required.
addpath functions;
addpath functions\SACcode;  % You need a path to the correlation index code. 
                            % Code is from Joris.    

% Load up spike times. 
% Split across 4 files to fit on github.
allSpikeTimes88and91 = loadSpikeTimes;

% Loads unit info and other supporting data.
load('rawdata\unitList');

% Options - which statistics to recalculate.
calcVS = true;     % Calculate vector strength statistics.
calcCI = true;     % Calculate correlation index statistics.
calcML = true;     % Calculate mode-locking statistics.
calcGC = true;     % Calculate reliability and envelope fluctuation statistics from Gai and Carney.
calcSACpeaks = false; % Calculate the SAC peaks from spike trains.  
                   % Takes weeks to run, and occasionally crashes. Output provided on GitHub.  
                   % Only set to true if you are brave. 

%% -----------------------------------------------------------------------
%         This calculates statistics for every condition, every unit. 

win100 = [20 100];          % N.B. Any tones >100ms are only analysed to 100ms. 
win50  = [20 50];
numPHbins = 52;
phbin = 0:1/numPHbins:1;    % see Malone et al. 2007
coinc_win = .05;            % Coincidence window for Correlation Index.
Ndata = length(allSpikeTimes88and91);
SACfigdir = 'SACfigs';
startii = 1;

% Variables for the peak analysis.
global CRITERION PLOTOPT;
CRITERION =.05;
PLOTOPT = 0;
SACpeaks_minspkpersweep = 4;    % We arbitrarily don't run the SAC peak picking
                                % on anything less than 4 spikes per
                                % sweep.

% Initialise outputs for which ever are being calculated.
if calcVS
    VS = cell(size(allSpikeTimes88and91));
end;
if calcML
    ML = cell(size(allSpikeTimes88and91));
end;
if calcCI 
    CI = cell(size(allSpikeTimes88and91));
end;
if calcGC 
    GC = cell(size(allSpikeTimes88and91));
end;
if calcSACpeaks
    SACpeaks = cell(size(allSpikeTimes88and91));
    
    % Pick out where it last crashed.
    load('datafiles\SACpeaks.mat');
    startii = max(find(cellfun(@(x) ~isempty(x), SACpeaks))) +1;
    
    % If not recalculating load up the correlation index, which provides
    % the SAC for the main analysis.
    if ~calcCI
        load('datafiles\SpikeStats.mat','CI');
    end;
end;

completedvector = cell(size(allSpikeTimes88and91));

for ii = startii:Ndata
    theseSpikeTimes = allSpikeTimes88and91{ii};
    duration = wholeList171212(ii,14); 

    fprintf('Data %d /%d:',ii,Ndata);
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
        if calcVS
            fprintf('VS,');
    
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
        end;

        % ----------- mode-locking -------------------
        if calcML
            fprintf('ML,');
            % Call getModelockCriteria to perform the full set of tests.
            ML{ii} =getModelockCriteria(theseSpikeTimes, thisModFreq, win2use, 1); 
            ML{ii}.modFreq = thisModFreq;
        end;

        % ------------- correlation index ----------------
        
        if calcCI
            % Calcuate the correlation index. 
            fprintf('CI');
            [h, bc] = SPTCORR(spktrainset,'nodiag',20,coinc_win,duration,'LouageNorm');
            CI{ii}.unsmoothedSAC = h;
            CI{ii}.CIvalue_0lag = h(bc==0);    % CI is nominally the  zero lag value.
            CI{ii}.CIvalue_max = max(h);       % But that might not be the maximum value.
            CI{ii}.lags = bc;
            
            % Calculate significance p-value associated with the CI,
            % using the Joris method. 
            if sum(~cellfun('isempty', spktrainset))>0
                % remove any empty cells
                spktrainset = spktrainset(~cellfun('isempty', spktrainset));
                pval = sacpeaksign(spktrainset,coinc_win,500,duration);  % Very slow.
                CI{ii}.p = pval;
            else
                CI{ii}.p =  NaN;
            end
        end;


        % ------------- SAC peak statistics ----------------
        
        % An elaborate analysis of the SAC which looks for peaks additional
        % to the expected ones fom phase-locking. VERY slow. 
        % Cannot emphasize enough that running this is a big deal. 
        if calcSACpeaks
            % Calculate the correlation index. 
            fprintf('SAC peaks...');
            [SACpeaks{ii}, details] = getSACpeaks(spktrainset,CI{ii},win2use, ...
                thisModFreq, SACpeaks_minspkpersweep,SACfigdir,ii);
            if isfield(SACpeaks{ii},'amfreq')
                 SACpeaks{ii}.details = details;
            end;

            % save every 1000th dataset. Or set smaller for paranoia.
            if rem(ii,1000) == 0
                save datafiles\SACpeaks SACpeaks;
            end;
        end;

        % ------------- Gai & Carney statistics ----------------
            
        % Calculates: "reliability" (a PSTH correlation measure) and  
        % "envelope fluctuation".
        if calcGC
            GC{ii} = getGaiCarneyStats(spktrainset,win2use);
        end;

        fprintf('.\n');
    end;
    
end

% Save the calculated statistics. 
save datafiles\SpikeStats_tmp;
% Run and save separately, and save incrementally as they take so long (set all other statistics to false.
%save datafiles\SACpeaks SACpeaks;

