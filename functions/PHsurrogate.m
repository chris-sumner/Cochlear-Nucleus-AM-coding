function [st_sur, ISI_sur, ph_sur, ISIcell_PHsur] = PHsurrogate(st, per, spikePerTrial, ISI, toplot)
% PHsurrogate Produce surrogate spike trains preserving phase distribution 
%
% This function produces a new set of spiketrains which have the same spikes 
% but shuffled between periods, such that the number of spikes in each
% individual sweep is preserved. The phase distribution is guaranteed to be
% identical. The interspike interval distribution would be almost unchanged 
% if the spike trains were pure homogeneous Poisson in behaviour, the ISI will 
% be any altered if there is any inhomogeneity (e.g. serial dependence of 
% spiking) in the original spike trains. 
%
% [st_sur, ISI_sur, ph_sur, ISIcell_PHsur] = PHsurrogate(st, ISI, per, spikePerTrial, toplot)
%      st: vector of spike times in ms re: stimulus onset. 
%      per: period of modulation in ms. 
%      spikePerTrial: a vector of the number of spikes in each trial.
%      ISI: vector of ISIs in ms. Only needed for plotting.
%      toplot: Makes a plot of the process if true
%
% Original function writted by Jonathan Laudansiki (~2009), modified by
% Chris Scholes (2012) and Chris Sumner (2024).
%

numSweeps = length(spikePerTrial);

% ----------- initial phase calculation ------------

% This separates out the phase of the spike and which period of the 
% sweep they appear in. 
ip = floor(st'/per);        % Index of the period for each spike (1st per, 2nd per, ...).
ph  = (st' - ip*per)/per;   % Phase of each spike (ph1 ph2 ph3...).
numSpike = length(ph);

% ----------- shuffle all the spike phases --------------

% This is a completely random ordering of the phases
ph_sur = ph(randperm(numSpike));    % random perm of integer indexes

% ------ Creating surrogate indexes of periods  ---------

% This assigns the spikes to particular sweeps (stimulus presentations)
% and periods. 

ALLip = ip;                                 % A pool of all period index in the data.
NBipPOOL = length(ALLip);   
idx_repe = [0,cumsum(spikePerTrial)];

st_sur = [];ISI_sur = [];
for ii=1:numSweeps
    % For each sweep....

    % pick the correct number of shuffled phases for this sweep
    TEMPph_sur = ph_sur(idx_repe(ii)+1:idx_repe(ii+1)); 

    % Randomise which period a spike occurred in.
    % Randomly pick ip NBspt in the pool
    numSpikeInSweep = spikePerTrial(ii);            % Number of Spk Per Trial (i.e. NBspt).
    idxPICKED = randperm(NBipPOOL);                 % random indexes from 1:number of ips       
    idxPICKED = idxPICKED(1:numSpikeInSweep);       % pick the first numSpikeInSweep inds
    ip_temp = ALLip(idxPICKED);                     % get the ips       
    ip_temp = sort(ip_temp);                        % sort the ips

    % Reconstruct a spike train with absolute times.     
    TEMPst_sur = sort((TEMPph_sur+ip_temp)*per);    % Sort the spkt of a single surrogate trial.
    % The phases, AND which period the spike train appears in 
    % are both randomised. I am not convinvced you need both but this is
    % what Jonathan did (CS 11/7/24). 

    % Remove from the pool the ip selected
    idx2removeIP = ones(size(ALLip));
    idx2removeIP(idxPICKED) = 0.0;
    idx2removeIP = logical(idx2removeIP);
    ALLip = ALLip(idx2removeIP);
    NBipPOOL = length(ALLip);
    
    % OUTPUTS
    st_sur = [st_sur TEMPst_sur];               % Accummulated vector of spike times from each sweep.
    ISI_sur=[ISI_sur diff(TEMPst_sur)];         % Accumulated vector of ISIs
    ISIcell_PHsur{ii}=diff(TEMPst_sur);         % Cell array version of the ISIs - WHY HAVE BOTH?!
end    

% ------------------- plot, if toplot is true -----------------

if(toplot)
    subplot(421)
    plotISIscatt(st,[0 25]);
    subplot(422)
    plotISIscatt(st_sur,[0 25]);
    subplot(423)
    hist(ISI,[0.5:0.5:25]);xlim([0 25]);
    subplot(424)
    hist(ISI_sur,[0.5:0.5:25]);xlim([0 25]);
    subplot(425)
    hist(ph,[0:0.05:1]);xlim([0 1]);
    subplot(426)
    hist(ph_sur,[0:0.05:1]);xlim([0 1]);
end