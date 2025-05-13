function output = VSControlModel(pars,classifier)

% -------------------------------------------------------------------------

pars.srate = 10020;                       % Odd sample rate prevents artifacts in period histogram
phbins = 0:1/32:(1-1/32);
win2use = [0 80];


calcML= true;

GC = [];
CI = [];
ML = [];
SACPeaks = [];

% This make a simple set of Poisson spike trains. Just like the expt. 
% This uses von Mises distribution as the basis for the time varying
% proability. Based on Kessler et al. 2021. 
for fi = 1:size(pars.amfreqs,1)    
    fprintf('AM frequency: %d\n',pars.amfreqs(fi));
    
    spiketimes(1:pars.nsweeps,fi) = simVS(pars.amfreqs(fi), ...
        pars.phaselocking(fi,1),pars.phaselocking(fi,2), ...
        pars.dcrate(fi),pars.srate,pars.dur_s,pars.nsweeps);
    
    % Implement a 0.8ms deadtime.
    %dt = diff(spiketimes{ri,fi});
    %spiketimes{ri,fi} = spiketimes{ri,fi}(find(dt>0.8e-3)+1);

    % Compute relative phases. N.B. spike times in milliseconds. 
    spikephase{fi} = rem([spiketimes{1:pars.nsweeps,fi}],1e3/pars.amfreqs(fi))/ ...
        (1e3/pars.amfreqs(fi));
    % Compute the period histogram.
    ph{fi} = hist(spikephase{fi},phbins);
    % Mean polar coords of phase. 
    polar(fi) = mean(cos( 2*pi*spikephase{fi} ) + i*sin(2*pi*spikephase{fi}));
    
    % Compute the Gai and Carney statistics
    if isfield(pars,'calcGC') && pars.calcGC == true
        fprintf('\tGai & Carney statistics\n');
        GC{fi} = getGaiCarneyStats(spiketimes(:,fi),win2use);
    end;

    % Correlation index and SAC peaks
    if isfield(pars,'calcSACpeaks') && pars.calcSACpeaks == true
        fprintf('\tSAC peaks... ');
        [CI{fi}, SACPeaks{fi}] = calcCIandSACpeaks( spiketimes(:,fi), win2use, pars.amfreqs(fi),true,true);
    elseif isfield(pars,'calcCI') && pars.calcCI == true
        fprintf('\tCI\n');
        [CI{fi}, SACPeaks{fi}] = calcCIandSACpeaks( spiketimes(:,fi), win2use, pars.amfreqs(fi),true,false);
    end;        
    
    % mode-locking 
    if isfield(pars,'calcML') && pars.calcML == true
        fprintf('\tML,');
        % Call getModelockCriteria to perform the full set of tests.
        %try 
        ML{fi} =getModelockCriteria(spiketimes(:,fi), pars.amfreqs(fi), win2use, 1); 
        %catch end;
        ML{fi}.modFreq = pars.amfreqs(fi);
    end;
           
end;

%spiketrainset{1} = spiketimes;
stimset = ones(size(spiketimes,1),1)*pars.amfreqs';

% WR classifier.
if nargin<2
    classifier.tau =1;
    classifier.iterations = 1000;
    classifier.store_distances = 0;
end;

% Run the model.
fprintf('WR classifier...\n');
wr_test1 = WR_Classifier(spiketimes,stimset,classifier);

% Assign outputs.
output.wr = wr_test1;
output.GC = GC;
output.ML = ML;
output.CI = CI;
output.SACPeaks = SACPeaks;
output.ph = ph; 
output.vs(:,1) = abs(polar);
output.vs(:,2) = angle(polar);


function [CI, SACpeaks] = calcCIandSACpeaks( spktrainset, win2use , thisModFreq,calcCI,calcSACpeaks)
% Code taken from SpikeStatsCalc. 

% Variables for the peak analysis.
global CRITERION PLOTOPT;
CRITERION =.05;
PLOTOPT = 0;
SACpeaks_minspkpersweep = 4;

coinc_win = .05;            % Coincidence window for Correlation Index.
duration = diff(win2use);             

SACfigdir = [];
  

% ------------- correlation index ----------------

if calcCI
    % Calcuate the correlation index. 
    fprintf('CI... ');
    [h, bc] = SPTCORR(spktrainset,'nodiag',20,coinc_win,duration,'LouageNorm');
    CI.unsmoothedSAC = h;
    CI.CIvalue_0lag = h(bc==0);    % CI is nominally the  zero lag value.
    CI.CIvalue_max = max(h);       % But that might not be the maximum value.
    CI.lags = bc;

    % Calculate significance p-value associated with the CI,
    % using the Joris method. 
    if sum(~cellfun('isempty', spktrainset))>0
        % remove any empty cells
        spktrainset = spktrainset(~cellfun('isempty', spktrainset));
        pval = sacpeaksign(spktrainset,coinc_win,500,duration);  % Very slow.
        CI.p = pval;
    else
        CI.p =  NaN;
    end
end;


% ------------- SAC peak statistics ----------------

% An elaborate analysis of the SAC which looks for peaks additional
% to the expected ones fom phase-locking. VERY slow. 
if calcSACpeaks
    % Calculate the correlation index. 
    fprintf('SAC peaks...');
    [SACpeaks, details] = getSACpeaks(spktrainset,CI,win2use, ...
        thisModFreq, SACpeaks_minspkpersweep,SACfigdir,1);
    if isfield(SACpeaks,'amfreq')
         SACpeaks.details = details;
    end;

end;




