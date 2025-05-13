function outdata=getModelockCriteria(st, fm, win, doNHPP, ntrains, nsweepsPerTrain, plFlag)
% getModelockCriteria Computes a set of criteria used to test for
%                     mode-locking. Laudanski et al. (2010). 
%
% outdata = getModelockCriteria(st, fm, win, doNHPP, ntrains, nsweepsPerTrain, plFlag)
%   st: cell array with each element being a vector of spike times of a response
%       to a single presentation of the stimulus. Time in milliseconds post stimulus onset.
%   fm: Moduation frequency in Hz.
%   win: vector of length 2 with the mimumum and maximum spike time to include in 
%       calculation.
%   doNHPP: Flag to simulate a non-homogeneous Poisson process (1ms
%           dead-time + minimum 1ISI). 0(no)/1(yes). default = 1 (yes).
%   ntrains: 
%   nsweepsPerTrain:
%   plFlag: Whether to process any data handed in (1) or to only 
%           process data where phase-locking to the envelope frequency is
%           significant (0; default), determined by a Rayleigh statistic > 13.8.  
%
% Function originally written by Jonathan Laudanski (2009), modified by
% Chris Scholes (2012/13), updated and documented by Chris Sumner (2024).

%  ----------- Set up parameters and constants -------------- 

% Absolute refractory time of 1.0 ms for NHPP surrogate spike trains.
refrac = 1.0;       

% Length of 1 period in ms. 
period = 1000/fm;

% Flag for wwhether to run the non-homogenous Poisson process simulation.
% 0 means only a homogeneous version is run. 
if nargin<4
    doNHPP = 1;
end

% Number of times to recompute the Z-score.
% Its a noisy measure if there are small number of spikes.
% Compute multiple times to check reliability.
numZiterations = 10;

% Bootstrap parameters for Z-score estimation
if nargin<6
    ntrains = 100;        % Essentially the number of times we 
                          % resample to calculate the internal variability.
    %nsweepsPerTrain = 4; % should be less than 50% of nsweeps
    nsweepsPerTrain = 9; % The number of stimulus sweeps used at each resample.
                         % should be less than 50% of nsweeps
                         % This works for most of the data.
                         % If 50% rule is exceeded, Z is not computed.
end

% To only calculate stats for phase-locking units
if nargin<7
    plFlag=0; % Ray>13.8 for stats to work
else
    plFlag=1; % stats done on everything
end

% Calculate the duration of the analysis window
winDur = win(2) - win(1);

% ------------- Preprocessing & Initialisations --------------

% Sets up a default structure for each PH test.
tmpphtest = struct('st_PHsur',nan,'ISI_PHsur',nan,'ph_PHsur',nan,'ISIcell_PHsur',{},'Z',nan);
NHPP_fixeddt = tmpphtest;   % Fixed dead time
NHPP_1percdt = tmpphtest;   % Dead time which depends on the ISIs in the spike train - set to 1st percentile of ISIs 
                            % This is close to the minimum ISI but not quite as strict.

% Convert spike trains into a 2D cell of number of sweeps x max number of
% spikes. Each row is a sweep. Empty spike times are nans.
% C Scholes did this make it more convenient to manipulate. 
nSweeps = length(st);                               % Number of "sweep" - i.e. # of repetitions of stimulus
stInWin = sTimesCellTo2D(st,nSweeps);

% Remove spikes that are not in the time window
%nSweeps = size(st,2); % Old way. DELETE LATER
%stInWin = st;
stInWin(stInWin<win(1)|stInWin>win(2))=nan;
% put spiketimes in to 1D array
stInWin1D = stInWin(~isnan(stInWin));

% Calculate various statistics used in the analysis
%nSweeps = size(st,2);                               % Number of "sweep" - i.e. # of repetitions of stimulus
sweepsWithSpike = find(nansum(stInWin)~=0);         % Indexes to sweeps that contain at least one spike.
numSweepsWithSpike = length(sweepsWithSpike);       % Number of sweeps that contain at least one spike.
spikesPerSweep = nansum(stInWin./stInWin);          % Number of spikes in each sweep. 
sweepsWithOverOneSpike = find(spikesPerSweep>1);    % Index of sweeps with more than one spike (ISIs are nans if not)
numSweepsWithOverOneSpike = length(sweepsWithOverOneSpike);         % Number of sweeps with more than one spike
periodsPerWindow = winDur/period;                                   % Number of modulation periods within the analysis window.
spikesperperiod = length(stInWin1D)/(periodsPerWindow*nSweeps);     % Average number of spikes in a period.
spikespersweep = length(stInWin1D)/nSweeps;                         % Average number of spike in a sweep.

% ------------ ISI (interspike interval) calculations --------------

% Calculate ISIs (interspike intervals) of spikes times
ISI = diff(stInWin);
% put ISIs in to 1D array
ISI1D = ISI(~isnan(ISI));

% Constructing a cell of the ISIs in each trial
ISIcell = cell(1,nSweeps);
for ii=1:nSweeps
    if ii<=length(ISI)
        thisISI = ISI(:,(ii))';
        % Any sweeps with only one spike end up empty. 
        ISIcell{ii}=thisISI(~isnan(thisISI));
    end;
end

% ------------ Phase (period histogram) calculations --------------

% Convert all spike times to a vector phase within the period 
ph = (stInWin1D - floor(stInWin1D/period)*period)/period;
% Compute vector strength.
VS = sqrt(sum(cos(2*pi*ph))^2+sum(sin(2*pi*ph))^2)/length(ph);
% Compute Rayleigh statistic. 
R = 2*length(stInWin1D)*VS^2;      

% ------------ Miscellaneous calculations --------------

% The entrainment index (as per Kalluri)
EI = sum(ISI1D<fm*1.5)/(period*nSweeps);
% This is not intrinsic to the analysis. 


% ------------------------ Phase shuffling ----------------------------

% Check whether the anlysis should be done
% Gnerally if there is not significant phase locking in the first place,
% there is little point looking at mode-locking. 
if plFlag
    goAhead = ~isempty(stInWin1D);   % Requires there to be some spikes
else
                                     % Requires phaselocking be significant
                                     % and there to be at least one
                                     % interval.
    goAhead = R>=13.8 && ~isempty(stInWin1D) && ~isempty(ISI1D);
end

% Do it if we should.
if goAhead

    % Generate a surrogate spike train which has the same period
    % distribution but the spike phases are shuffled both between periods
    % both within and across sweeps.
    [st_PHsur, ISI_PHsur, ph_PHsur, ISIcell_PHsur] = PHsurrogate(stInWin1D, period, spikesPerSweep,ISI1D, 0);

    % The Z-score is a measure of the difference in ISI distribution before
    % and after shuffling, normalised to the variability of the unshuffled 
    % distribution.
    % ztest is a wrapper function for getZscore_groupsweeps.
    [Z, zDist, zTest, zTerms] = ztest(ISIcell,ISIcell_PHsur,ntrains,nsweepsPerTrain,numZiterations);
    spikesperperiod_PHsur = length(st_PHsur)/(period*nSweeps);
    % If Z = nan a Z score could not be computed. Normally not enough
    % spikes.
    
    % A kolmogorov-smirnov test is used to see if the interspike intervals
    % distributions are different as a result of shuffling. 
    % In practice this is not a sensitive measure - the difference is 
    % significant in most neurons. We did not know that in Laudanski et al.
    % (2010). 
    highPassCrit = 1;                   % Do not test on internals below 1ms. 
                                        % This makes it a little less
                                        % senstive because all neurons have
                                        % a refractory period. 
    allISIs = [ISIcell{:}];                         % Put all the data ISIs together.
    allISISurs = [ISIcell_PHsur{:}];                % Put all the shuffed ISIs together.
    allISIs = allISIs(allISIs>=highPassCrit);               % Remove the really short intervals.
    allISISurs = allISISurs(allISISurs>=highPassCrit);      % Remove the really short intervals.
    
    if ~isempty(allISIs) && ~isempty(allISISurs)
        [ksRes,ksP,ksStat] = kstest2(allISIs,allISISurs,0.01);  % Perform the test at alpha = 0.01.
    else
        ksRes = nan;
        ksP = nan;
        ksStat = nan;
    end;
    
    % Entrainment index measures
    % the reliablity with which a neuron fires in each phase ie. 1 = fires
    % in every period, 0.5 = misses every second period.
    EI_ISI_PHsur = sum(ISI_PHsur<fm*1.5)/(period*nSweeps);
    
    if doNHPP
        % ********* Simulate a non-homogeneous Poission train **********

        % This provides an additional comparison for spike trains, against 
        % models which are not pure Poisson but do not exceed 
        % refractory behaviour.

        if floor(diff(win)/period)~=0 % if fmod is too low then get outta here
            
            % for input to NHPP function - should not need to recompute this 
            %ip = floor(stInWin1D'/period);        % Index of the period for each spike (1st period, 2nd period, ...).
            %ph  = (stInWin1D' - ip*period)/period;   % Phase of each spike (ph1 ph2 ph3...).
            
            % With a fixed refactory period. 
            NHPP_fixeddt = NHPPsurrogate(ph, period, ISI1D, nSweeps, refrac, win);
            NHPP_fixeddt.Z = ztest(ISIcell,NHPP_fixeddt.ISIcell,ntrains,nsweepsPerTrain);
            NHPP_fixeddt.EI = sum(NHPP_fixeddt.ISI<(fm*1.5))/(period*nSweeps);
            NHPP_fixeddt.spikesperperiod = length(NHPP_fixeddt.st)/(period*nSweeps);
            
        end
    end
else
    st_PHsur = NaN;
    ISI_PHsur = NaN;
    ph_PHsur = NaN;
    ISIcell_PHsur = {};
    Z = NaN;
    zTerms = NaN;
    zDist = NaN*ones(1,numZiterations);
    EI_ISI_PHsur = NaN;
    spikesperperiod_PHsur = NaN;
    zTest = nan;
end;

outdata.stats.names={'VS','VS_Isur','R','R_Isur','Z','ksRes','ksP','ksStat'};
outdata.st = stInWin1D;
outdata.zDist = zDist; 
outdata.zTerms = zTerms;
outdata.ISI = ISI1D;
outdata.ph = ph;
outdata.ModeLockingTest_Z = zTest;

% PH surrogate info
outdata.PHsur.st = st_PHsur;
outdata.PHsur.ISI = ISI_PHsur;
outdata.PHsur.ph = ph_PHsur;
outdata.PHsur.ISIcell = ISIcell_PHsur;

% --------------- Interspike Interval Shuffling ----------------

% In Laudanski et al. (2010), we considered that in "mode-locking" neurons, 
% ISI structure depends on factors other than stimulus phase, so this 
% shuffling is likely to deystroy the phase distribution. We (fairly 
% arbitrarily) suggested that a non-significant Rayiegh statistic following
% shuffling is one of two criteria for mode-locking. 

% If the Z-score calculation appears successful (i.e. did not output a nan), 
% perform the ISI shuffling.
if ~isnan(Z)
    % this includes the Compulsory Stop msg
    [st_Isur, ISI_Isur, PH_Isur, VS_Isur, Rs] = serialISIshuffle(stInWin1D,ISI1D,period,0);
    
    outdata.stats.vals=[VS, VS_Isur, R, Rs, Z, ksRes,ksP,ksStat];
    outdata.Isur.st = st_Isur;
    outdata.Isur.ISI = ISI_Isur;
    outdata.Isur.ph = PH_Isur;
    outdata.ModeLockingTest_ISIshuf = double(Rs<13.8);
else
    outdata.stats.vals=[VS, NaN, R, NaN, Z, NaN, NaN, NaN];
    outdata.Isur.st = NaN;
    outdata.Isur.ISI = NaN;
    outdata.Isur.ph = NaN;
    outdata.ModeLockingTest_ISIshuf = nan;
end;

%Add the extra tests to the output.
outdata.NHPP.fixedDt = NHPP_fixeddt;

% spikecount info
outdata.spikecounts.total = length(stInWin1D);
outdata.spikecounts.numSweeps = nSweeps;
outdata.spikecounts.period = period;
outdata.spikecounts.spikesperperiod = spikesperperiod;
outdata.spikecounts.spikesperperiod_PHsur = spikesperperiod_PHsur;

% EI info
outdata.EI.val = EI;
outdata.EI.ISI_PHsur = EI_ISI_PHsur;




% -------------------------------------------------------------------------
%                           SUPPORTING FUNCTIONS
% * -----------------------------------------------------------------------

% -------- ztest - top level function for computing Z-scores --------------
function [Z, wholeZdist,zTest, zTerms] = ztest(ISIcell, ISIcell_PHsur,ntrains,nsweepsPerTrain,zIter)
% ztest: Function to generate a Z-score.

% Initialise.
Zterms = nan;

% ----------- Call the main function to compute Z-score -----------

% Check equal # of interval in ISI and surrogate.
if length(ISIcell)>1 && length(ISIcell_PHsur)==length(ISIcell)

    maxint = max([ISIcell{:} ISIcell_PHsur{:}]);        % Get longest ISI of all. 

    % If iterations is no specified (zIter)
    if nargin<5
        [Z, zTest, zTerms] = getZscore_groupsweeps(ISIcell,ISIcell_PHsur,maxint(1),ntrains,nsweepsPerTrain);
        wholeZdist = NaN;
    else
    % If iterations are specified we bootstrap on the Z calculation.    
        wholeZdist = zeros(1,zIter);  
        wholeZdist(:) = NaN;
        for ii = 1:zIter
            [dTmp, testTmp, tTmp] =  getZscore_groupsweeps(ISIcell,ISIcell_PHsur,maxint(1),ntrains,nsweepsPerTrain);
            ZTests(ii) = testTmp;
            wholeZdist(ii) = dTmp;
            allZTerms(ii) = tTmp;
        end
        Z = nanmean(wholeZdist);
        % Z terms are only carried through as the mean.
        zTerms.mean_RMSe_comp = nanmean([allZTerms(:).mean_RMSe_comp]);
        zTerms.mean_RMSe_data = nanmean([allZTerms(:).mean_RMSe_data]);
        zTerms.std_RMSe_data = nanmean([allZTerms(:).std_RMSe_data]);
        zTest = round(mean(ZTests));
    end
else
    Z = NaN;wholeZdist=NaN;zTerms = nan;zTest = nan;
    disp([length(ISIcell) length(ISIcell_PHsur)])
end;


% -------- minISI - compute min ISI at a percentile --------------
function [minisi, maxisi] = minISI(ISI2test,perc)
mISI = min(ISI2test);
MISI = max(ISI2test);
xISI = [mISI:(MISI-mISI)/100:MISI];
Nisi= histc(ISI2test,xISI);
idxmin = find((Nisi/sum(Nisi))>=perc,1,'first');
idxmax = find((Nisi/sum(Nisi))<=(1-perc),1,'last');

minisi = xISI(idxmin);
maxisi = xISI(idxmax);
