function [s_st, sISI, sPH, sVS, sR] = serialISIshuffle(st,isi,per,toplot)
% serialISIshuffle - shuffles intervals in a set of spike trains
% 
% Generates a surrogate spike train from the data with similar
% pair conditional probability. If the spike train resemmbes a Poisson
% process, the phase-distribution within the modulation period will
% remain intact. Returns urrogate spk times, ISIs, phases, 
% VS, and Rayleigh stat.
%
% In Laudanski et al. (2010), we considered that in "mode-locking" neurons, 
% ISI structure depends on factors other than stimulus phase, so this 
% shuffling is likely to destroy the phase distribution. We (fairly 
% arbitrarily) suggested that a non-significant Rayiegh statistic following
% shuffling is one of two criteria for mode-locking. 
%
% [s_st, sISI, sPH, sVS, sR] = serialISIshuffle(st,isi, per,toplot)
%      st: Vector of spike times (ms) relative to stimulus start.
%      isi: Vector of isi (in ms). 
%      per: the period of the modulation in milliseconds.
%      toplot: whether (1) or not (0) to plot. 

% Calculate the spike phases. Only needed for display.
ph  = (st - floor(st/per)*per)/per;

% We need interval pairs.
isiPRE = isi(1:end-1);  % Intervals 1->N-1
isiPOST= isi(2:end);    % Intervals 2->N

% Begin by picking an interval at random.
rdn_idx = randperm(length(isi));
sISI(1) = isi(rdn_idx(1)); % choose a random ISI

% ----------------- Interval shuffling -------------------------

% For each interval we pick, we find all the intervals were
% the preceding interval was similar - the same or a bit longer.
% By a bit, we choose maximum dISI to be random between 0.5 and 1.8 ms
% longer than the interval we just picked. 
for ii=2:length(isiPRE)
    % find all ISIs that are greater than the last ISI and less than a random ISI
    % + random number between 0.5 and 1.8
    isi_idx = (isiPRE>=sISI(ii-1) & isiPRE<(sISI(ii-1)+ 1.3*rand+0.5));
    % Store all the ISIs which follow the intervals we just found.     
    tempISI = isiPOST(isi_idx); 
    % Randomly pick of those new intervals. 
    % N.B. This is a resample with replacement as we do not remove this 
    %      interval from isiPOST. 
    if(~isempty(tempISI))
        idx = randperm(length(tempISI));
        sISI(ii) = tempISI(idx(1)); % then pick a random ISI from the list of ISIs
    else
        % Stops if it fails to find a pair of intervals to satisfy.
        disp(['Compulsory stop at:' num2str(ii)])
        break;
    end
end

% -------------- Reconstruct the shuffled spike train. --------------

% We reconstruct the spike time by accumulating intervals. 
% If the interval statistics really reflect the phase locking properties
% then phase locking will not be completely destroyed. 

% Note that this is reconstructed as one long spike train. Much longer than
% the original stimulus.,
s_st = cumsum([st(1) sISI]);

% Compute statistics from it. 
sPH = (s_st - floor(s_st/per)*per)/per;
sVS = sqrt(sum(cos(2*pi*sPH))^2+sum(sin(2*pi*sPH))^2)/length(sPH);
sR = 2*length(s_st)*sVS^2;

% Plot
if(toplot)
    subplot(221)
    hist(ph,50);ylim([0 200])
    title('DATA period-histo')
    subplot(222)
    scatter(isiPRE,isiPOST,'x');xlim([0 50]);ylim([0 50]);
    subplot(223)
    hist(sPH,50);ylim([0 200])
    title('Conditional ISI shuffling period-histo')
    subplot(224)
    scatter(sISI(1:end-1),sISI(2:end),'x')
    xlim([0 50]);ylim([0 50]);
end