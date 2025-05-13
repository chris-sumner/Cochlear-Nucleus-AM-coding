function [output] = getGaiCarneyStats(spktrainset,win2use)
% getGaiCarneyStats calculates PSTH reliability and envelope fluctuations
% 
% Replicates analysis measures used in Gai & Carney (2006): 
%  Reliability is the correlation between smoothed PSTHs made from odd and 
%  even presentations of the stimulus.
%  Envelope fluctuation is the normaised sum of the absolute differential of
%  the smoothed PSTH.
%
%  [output] = getGaiCarneyStats(spktrainset,win2use)
%       spktrainset: cell array of spike times from each sweep.
%       win2use: vector of two values - minumum and maximum times to analyse.    
%                N.B. This determiines the PSTH bins for the analysis but
%                     spikes are not removed outside this window.
%       output: Structure of:
%               temporalWindows: vector of the window sizes used.
%               binwidths: bin wdiths of resulting PSTHs (acutally the
%                       same)
%               splitpsths: cell array of psths for each temporal window.
%               reliability_R: Peason correlation coefficients between
%                       PSTHs for each temporal window.
%               envFluct: Envelope fluctuation statistics - one for each 
%                       temporal window.
%
% Written by Chris Sumner (2024). 

% -------------- Reliablity and Envelope fluctuation --------------

% From Gai & Carney:
% PSTHs were convolved with a Gaussian smoothing function before
% the correlation was computed. The SD of the Gaussian smoothing
% function was referred to as the temporal analysis window, which was
% varied over a large range (0.2, 0.6, 1.6, 4.5, and 12.8 ms) to examine
% both slowly and rapidly changing temporal information. The Gaussian
% function was truncated to 3 SD in length (as in Martin et al. 2004). For
% each temporal analysis window, the PSTH was convolved with the
% Gaussian smoothing function, and the bin width of the resulting PSTH
% was chosen to match the temporal analysis window.
% The correlation of the PSTHs made from the odd and even trials.

% Analysis windows are much smaller than G&C because windows larger than 
% the fmods make no sense (and the s.d. is much smaller than the 
% effective period of the Gaussian).  
analWins = [0.08, 0.16, 0.24, 0.48, 0.64, 1.1, 1.6, 2.4, 3.2];
nTrains = length(spktrainset);
initialBinWidth = .25e-4;
maxT = max([spktrainset{:}]);

% Split the trials.
% N.B. This is not sensitive to which sweep spikes came from.
train1 = [spktrainset{1:2:end}];
train2 = [spktrainset{2:2:end}];

% Initialise the output structure for reliability.
output = struct();
output.temporalWindows = analWins;  
output.binwidths = analWins;      % These are the same but to avoid doubt.

% Make PSTHs at an arbitrarily small bin width. 
for ii = 1:length(analWins)

    % ---------------- Reliability -----------------

    % Make the bins and inialise arrays for this analysis window.
    bins = [win2use(1):analWins(ii):win2use(2)];
    psth1 = zeros(size(bins));
    psth2 = zeros(size(bins));
    winMax =  analWins(ii)*3;
    
    % Note about bin widths and sums. 
    % As the bin width is effectively 1sd, each spike will contribute
    % to 7 bins. The PDF across this range sums to ~1, so overall
    % each spike contributes 1 to the PSTH.
    
    % Loop through each bin.
    % This operates on raw spike times to make the PSTH in one go. 
    % Avoids the issue of making a PSTH, filtering and resampling,
    % which I was not clear about from the methods in Gai & Carney. 
    for bini = 1:length(bins)
        tmpt = normpdf( train1(  abs(train1-bins(bini))<winMax ) - bins(bini), ...
                        0, analWins(ii) );
        psth1(bini) = sum(tmpt);
    
        tmpt = normpdf( train2(  abs(train2-bins(bini))<winMax  ) - bins(bini), ...
                        0, analWins(ii) );
        psth2(bini) = sum(tmpt);        
    end;
    % Reliability is the correlation between the two psths
    Rtmp = corrcoef(psth1,psth2);
 
    % ------------------ Envelope fluctuation ------------------

    % Discharge times were smoothed by the Gaussian smoothing function
    % described above and combined across all stimulus presentations
    % to form the PSTH. The discharge count in the previous bin was
    % subtracted from the count in each bin of the PSTH. The absolute
    % values of all the subtracted numbers were summed. The fluctuation of
    % the PSTH was defined as the sum divided by the number of bins in the
    % PSTH for a given temporal analysis window.  
    
    % Since the filtering is linear we can add the two psths.
    psth0 = psth1 + psth2;
    envfluctuation = sum(abs(diff(psth0))) / length(psth0);

    % Make outputs
    output.splitpsths{ii} = [psth1; psth2];
    output.reliability_R(ii) = Rtmp(2);
    output.envFluct(ii) = envfluctuation;

end;


