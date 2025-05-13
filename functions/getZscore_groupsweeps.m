function [Z, ZTest, Zterms]= getZscore_groupsweeps(ISI,shuf_ISI,ISIlim,ntrains,nsweepspertrain)
% getZscore_groupsweeps compute a Z-score from multiple stimulus
%                       presentations.
%
% [Z, Ztest, Zterms]= getZscore_groupsweeps(ISI,shuf_ISI,ISIlim,ntrains,nsweepspertrain)
%      ISI: cell array of ISIs computed from each sweep of a stimulus.
%      shuf_ISI: cell array of ISIs, arranged into sweeps, but after phase
%                shuffling.
%      ISIlim: maximum interval to calcuate for (usualy just the real max). 
%      ntrains: the number of times to resample to compute the Z-score.
%               Larger number is more reliable but slower. Try 50-100 & 
%               check reliability. 
%      nsweepspertrain: the number of sweeps to use to make each spike
%               train. A larger number reduces the variance, increases the 
%               Z-score. Must be less than 1/2 the total number of sweeps.
%      Z: the final calculated Z-score.
%      Zterms: the invidual terms that make up the Z-score. Mean-squared
%              errors for comparing ISI histograms for the original 
%              spike-train, shuffled vs. original and the standard 
%              deviation of the original spike train. 
%     Outputs:
%           Z: The Z-score computed.
%           Zest: Whether Z-score > 2 (coded as 1,0 or NaN if not
%           computed).
%
% Computes a Z-score as per Laudanski et al (2010), but is generalised to
% be more usuable with short stimuli. 
%
% In Laudanski et all we suggested that a Z-score of 2 (2 s.d) was one of
% the criteria required for "mode-locking".
% 
% Modified from original code by Jonathan Laudanski (~2009), by Chris Scholes 
% (~2012) and Chris Sumner (2024).

% Initialise return parameters
Z=nan;
Zterms.mean_RMSe_comp = nan;
Zterms.mean_RMSe_data = nan;
Zterms.std_RMSe_data = nan;
ZTest = nan;

% Check that there are enough sweeps in the data to perform the 
% calculation.
if length(ISI)<nsweepspertrain*2
    return;
end;

% Randomly sample sets of ISIs ntrains times. 
% Each time pick only nsweepspertrain of the available sweeps.
for ii=1:ntrains
    inds = randperm(length(ISI));                                   % A random indexing of all the sweeps
    ISIresample{ii} = [ISI{inds(1:nsweepspertrain)}];               % A list of the sampling of ISIs from raw data. 
    shuf_ISIresample{ii} = [shuf_ISI{inds(1:nsweepspertrain)}];     % A list of the sampling of ISIs from phase shuffled data. 
end;
ISI = ISIresample;
shuf_ISI = shuf_ISIresample;

% Check they have the same number of trains (CS: Don't think this can happen). 
if(length(ISI)~=length(shuf_ISI))
    disp('Error: ISI and surrogate do not have same repe number!');
    disp([length(ISI) length(shuf_ISI)]);
    return;
end

% Make bins for histograms.
dISI = ISIlim/60;
x = 0:dISI:ISIlim;

% Inialise storage for histograms
Ndata = zeros(length(ISI),length(x));   % From data.
Nshuf = zeros(length(ISI),length(x));   % From period shuffled surrogate data.

% Compute histograms for each sample.
for ii=1:length(ISI)
    % If any of the histograms are empty then there are no enough 
    % spike to compute a Z-score - so get out. 
    if(isempty(histc(ISI{ii},x)) || isempty(histc(shuf_ISI{ii},x)))
        return;
    else
        Ndata(ii,:) = histc(ISI{ii},x);
        Nshuf(ii,:) = histc(shuf_ISI{ii},x);
    end
end

% We compare the difference (MSE - mean squared error) between a sample 
% taken from ISI and another sample taken from ISI against the
% difference (MSE) between a sample taken from ISI and a sample taken 
% from ISI shuf to give the Z score

% we do this for all possible combinatoins of samples.
% nchoosek gives all possible combinations of inds without comparing like
% for like (ie. ind 1 with ind 1)
comb = nchoosek(1:ntrains,2);    % List of all combinations in between 50 pres.
for ii = 1:nchoosek(ntrains,2)
    n1 = comb(ii,1);        % idx1 of the repetitions
    n2 = comb(ii,2);        % idx2 of the repetitions
    % ***   Testing diff between ISI p.d.f for presentation n1 and n2 ***
    MSe_data(ii)= mean((Ndata(n1,:)-Ndata(n2,:)).^2)...
        /sum(Ndata(n1,:))^2;
    % *** idem but against shuf_ISI p.d.f for presentation n2 ***
    MSe_comp(ii)= mean((Ndata(n1,1:end)-Nshuf(n2,1:end)).^2)...
        /sum(Ndata(n1,1:end))^2;
end

% Note of approximation made by Jonathan (from Chris 20/11/16):
% The denominators above are there to convert the ISI histograms to PDFs.
% However they are only approximate, as they assume the same number of
% intervals for both histograms. Precise would be to convert each to a PDF 
% before subtraction.   : 
%    MSe_data(ii)= mean( ( (Ndata(n1,:)/sum(Ndata(n1,:)))- ...
%                         (Ndata(n2,:)/sum(Ndata(n2,:)))  ).^2);
% However, I've left it as it was. 

% Z-score : Distance between the mean normalised by the variance.
MSe_comp(isinf(MSe_comp)) = nan;
MSe_data(isinf(MSe_data)) = nan;
mean_RMSe_comp = nanmean(sqrt(MSe_comp));
mean_RMSe_data = nanmean(sqrt(MSe_data));
std_RMSe_data = nanstd(sqrt(MSe_data),0);
Zterms.mean_RMSe_comp = mean_RMSe_comp;
Zterms.mean_RMSe_data = mean_RMSe_data;
Zterms.std_RMSe_data = std_RMSe_data;
Z = (mean_RMSe_comp-mean_RMSe_data)/std_RMSe_data;
ZTest = double(Z>2);

