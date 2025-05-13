function [peakstats, peakstats_details] = getSACpeaks(spktrainset,CIin,analwindow,amfreq,minspkpersweep,plotdir,unitid)
% getSACpeaks  main function to analyse the peaks in a single SAC. 
%  [peakstats, peakstats_details] = getSACpeaks(spktrainset,CIin,analwindow,amfreq,minspkpersweep,plotdir,unitid)

global PLOTOPT;

% If a plot directory is specified then save any plot.
if nargin>5 && ~isempty(plotdir) && PLOTOPT == true
    saveplot = true;
else
    saveplot = false;
end;

% Only plot am freqencies<1000 and up to index 10000
if PLOTOPT && amfreq<1000 && unitid<10000 && CIin.p>0.99
    plotopt = true;
else
    plotopt = false;
end;

% Reassure user something is happening.
fprintf('getSACpeaks: ');
if nargin>6
    fprintf(' (#%d) ',unitid);
end;
fprintf('AM frequency:%dHz. ',amfreq);  

% Analysis window.
tstart = analwindow(1);
tend = analwindow(2);
duration = diff(analwindow);

% Parameters for the peak picking. 
am_period = 1e3/amfreq;
thresholdpvalue = .05;  % DOES THIS DO ANYTHING?
randlagfactor = 2;      % ???
wfr = 1/50;             % Length of the smoothing window (in fraction of a period).
maxlagfactor = 1.5;     % Maximum lag to compute re period. N.B. Do not compute further than you need to.
                        % Lag should be much smaller than simulus duration.
                        % But too short can cause end effect problems.
maxlag_ms= maxlagfactor*1e3/amfreq;

% Initialise outputs.
peakstats_details = struct();
peakstats = struct();

% Check the CI input is not empty.
% Most likelt this would indicate not enough spikes.
if isempty(CIin)
    fprintf('Empty correlation index input. Stopping.\n');
    return;
end;

% Check the number of spikes
allspks = [spktrainset{:}];
numSweeps = size(spktrainset,1);
spikesPerSweep = length(allspks)/numSweeps;
if spikesPerSweep<minspkpersweep
    fprintf('Not enough spikes for analysis. Stopping.\n');
    return;
end;

% This speeds things up a bit. 
if amfreq>1000
    fprintf('AM frequency >1kHz. Stopping.\n');
    return;
end;

% Calculate window parameters for analysis. 
wlen = round(1e3*wfr/(amfreq*(CIin.lags(2)-CIin.lags(1))))+1;   % Smoothing window.
isel = find(abs(CIin.lags)<maxlag_ms*1.1);                      % Analysis window. A little more than 1 period. 
lagaxis = CIin.lags(isel);

% Smooth with a symmetrical n-point moving average. 
% N.B. Also restrict SAC to a little over 1 period. 
if wlen>1
    smoothsac = filtfilt(ones(1,wlen)/wlen,1,CIin.unsmoothedSAC(isel));
else
    smoothsac = CIin.unsmoothedSAC(isel);
end;

% This is the function which actually picks peaks in the smoothed SAC 
[peaksets, trofsets salsets] = findmainpeaks5(smoothsac);

if plotopt
    sacfigh = figure; set(gcf,'Position',[200 200 600 400]);
    % Plot the results
    sacploth = subplot(2,1,1);
    plot(lagaxis,smoothsac,'color',[.5 .5 1]); hold on;
    plot(lagaxis,smoothsac/10,':','color',[.7 .7 1]);
    plot(lagaxis,smoothsac/5,':','color',[.7 .7 1]);
    xlabel('Delay (ms)');
    titlestr = [' SAC for AM freq:' num2str(amfreq) 'Hz,  p(CI@0ms):' num2str(CIin.p) ];
    if nargin>6
        titlestr  = ['Dataset ID:' num2str(unitid) titlestr];
    end;
    ylabel('CI'); title(titlestr); set(gcf,'Name',titlestr);
    line([am_period am_period],[0 max(ylim)]);
    cellfun(@(x,y) plot(lagaxis(x),y,'+'),peaksets,salsets); 
    cellfun(@(x) plot(lagaxis(x),smoothsac(x),'o'),peaksets); 
    %plot(peaksets,salsets,'+');
end;


% -- Bootstrapping to get a null distribution for significance  ----

% The peak-picking works well, but it can pick very small peaks which could
% be noise. Here we run it for spike trains made from shuffled interspike
% intervals. This retains the first-order ISI statistics but destroys the 
% higher order relationships. 

randn = 500;
fprintf('\nComputing bootstrapped distributions of peak saliences:\n');
%clear  rand_spikesets randpeaksets randtrofsets randsaliences

randpeaksets = cell(randn,1);
randtrofsets = cell(randn,1);
randsaliences = cell(randn,1);
allrandsmoothsacs = cell(randn,1);

% This can be a parfor loop. 
for ri = 1:randn
    % Define variables
    rand_spikesets = cell(length(spktrainset),1);
    isis = cell(length(spktrainset),1);
     
    % --- Randomise intervals across all trials ---
    
    % Accumulate all intervals. 
    for sweepi = 1:length(spktrainset)                      
        % This is another way - maintain the ISIs. 
        isis{sweepi} =  diff(spktrainset{sweepi});
    end;

    % Convert to array and then shuffle.
    isiarray = [isis{:}];
    surrograteisis = isiarray(randperm(length(isiarray)));

    % Now convert back to spike times - same # spikes in each sweep. 
    counter = 1;
    for sweepi = 1:length(spktrainset)                      
        inc = length(spktrainset{sweepi})-1;
        if length(spktrainset{sweepi})==1
            rand_spikesets{sweepi} =  [spktrainset{sweepi}(1)];
        elseif length(spktrainset{sweepi})>1
            rand_spikesets{sweepi} =  cumsum([spktrainset{sweepi}(1) surrograteisis(counter:counter+inc-1)]);
            counter = counter+inc;
        end;
    end;        

    % Compute the SAC.
    % N.B. maxlag + 10% solves end effects problems. 
    [randh2plot randbc2plot] = SPTCORR(rand_spikesets,'nodiag',1.1*maxlag_ms,.05,duration,'LouageNorm');

    % Recompute smoothing window and analysis window for peakpeaking.
    % It is likely the SAC is shorter than the one handed in. 
    if ri == 1
        wlen = round(1e3*wfr/(amfreq*(randbc2plot(2)-randbc2plot(1))))+1;   % Smoothing window.
        isel = find(abs(randbc2plot)<maxlag_ms*1.1);                            % Analysis window.
    end;

    % Smooth the SAC function.
    if wlen>1
        randsmoothsac = filtfilt(ones(1,wlen)/wlen,1,randh2plot);
    else
        randsmoothsac = randh2plot;
    end;
    % Store all the smoothed sacs.
    allrandsmoothsacs{ri} = randsmoothsac;
    
    if rem(ri,50)==0 fprintf('.'); end;
end;    
 
% findmainpeaks5 will not work within the parfor loop. However taking it out does not slow things down significantly.  
for ri = 1:randn
    % Find the peaks and troughs in this random set. 
    [randpeaksets{ri}, randtrofsets{ri}, randsaliences{ri}] = findmainpeaks5(allrandsmoothsacs{ri});
    
    % Plot example random set.
    if ri<10 & plotopt
        %figure(sacfigh);
        subplot(sacploth);
       % plot(randbc2plot(isel),allrandsmoothsacs{ri}(isel),'color',[1 .7 .7]);
         plot(randbc2plot,allrandsmoothsacs{ri},'color',[1 .7 .7]);
    end;
end;
% Store one example.
randsmoothsaceg = allrandsmoothsacs{1};



% Statistics of the saliences (peaks and troughs) found in the random
% data. 
salvals = cell(120,1);
saltmp = [randsaliences{:}];
for i=1:length(saltmp)
    sallen =  length(saltmp{i});    
    
    %  --------- One way of computing the salience statistics --------
    % This gives you the full distribution of the saliences of all the
    % peaks given a number of picked peaks. It does not take account of
    % the fact that peaks will get smaller the more you pick. So it is 
    % conservative.

    %salvals{sallen}  = [ salvals{sallen} saltmp{i} ]; 
    % salval{1} is the saliences of all the peaks and troughs in the random data
    % if you only pick one peak.      
    % salval{2} is the saliences of all the peaks and troughs in the
    % random data when you pick two peaks - so by definition the largest
    % and next largest. 

    %  --------- Another way of computing the salience statistics --------
    % This is more generous - only stores the least salient peak
    % This is the basis for computing the probabity of finding a peak of
    % this size conditioned on the number of peaks you have already picked.
    salvals{sallen}  = [ salvals{sallen} min(saltmp{i}) ];    
end;

% Histogram of the saliences found for any given number of peaks.
salbins = [0:.01:10]; clear salhist salpval;
for i=1:length(salvals)
    if ~isempty(salvals{i})
        % Histogram of salience values. 
        salhist(i,:) = hist(salvals{i},salbins);
        % P values for a given number of peaks and salience value. 
        salpval(i,:) = 1- cumsum(salhist(i,:))/sum(salhist(i,:));
        
        % Some additional statistics about the saliences. 
        salmean(i) = nanmean(salvals{i});
        salmedian(i) = nanmedian(salvals{i});
        salsd(i) = nanstd(salvals{i});    
        salquantiles(i,:) = quantile(salvals{i},[.025 .25 .50 .75 .975]); % a useful summary of x            
    else
        salhist(i,:) = zeros(1,length(salbins)) ;
        salpval(i,:) = zeros(1,length(salbins)) ;

        salmean(i) = nan;
        salmedian(i) = nan;
        salsd(i) = nan;    
        salquantiles(i,:) = [nan nan nan nan nan];        
    end;
end;    

if plotopt
    % Make a figure of the salience distributions.
    salploth = subplot(2,1,2);
    colorspecs = get(gca,'ColorOrder');
    plot(salbins,salhist);
end;

% Calculate the saliences of the individual peaks and the p value of
% each.
peakpvalues = [];
saliences = cell(length(peaksets),1);
for pi = length(peaksets):-1:1
    % Compute the salience for each peak-trough pair.
    % N.B. this is only one way - need to compute the MINIMUM salience.
    setlen = sum(~isnan(peaksets{pi}));    
    saliences{setlen} = salsets{pi};

    % % From this salience compute the probability it would occur in
    % the randomised peak statistics.        
    for pii = 1:setlen
        pind = max(find(salbins<saliences{setlen}(pii)));
        if ~isempty(pind)
            peakpvalues{setlen}(pii) = salpval(setlen,pind);
        end;
    end;

    if plotopt
        subplot(salploth);
        hold on;
        colorspecs = get(gca,'ColorOrder');
        plot(saliences{setlen},pi*15*ones(1,length(saliences{setlen})), ...
            'o','color',colorspecs(rem(setlen-1,7)+1,:),'markerface',colorspecs(rem(setlen-1,7)+1,:));
        maxx =  max(saliences{setlen});
        if ~isempty(maxx)
            text(maxx*1.1,pi*50,sprintf('%d:%.2g',setlen,max(peakpvalues{setlen})),'fontsize',6);
        end;
    end;

end;

% More reassurance.
fprintf('\tdone.\n');

% Put all this information into a structure. 
peakstats.wlen =wlen;
peakstats.amfreq = amfreq;
peakstats.am_period = am_period;
peakstats.peaksets =peaksets;
peakstats.trofsets = trofsets;
peakstats.salsets =salsets;
peakstats.maxlag_ms = maxlag_ms;
peakstats.smoothsac =smoothsac;
peakstats.randsmoothsaceg =randsmoothsaceg;
peakstats.peakpvalues = peakpvalues;
peakstats.lagaxis = lagaxis;
peakstats.spikespersweep = spikesPerSweep;

% Statistics about the saliences of the random data.
peakstats.randsals.mean = salmean;
peakstats.randsals.median = salmedian;
peakstats.randsals.sd = salsd;
peakstats.randsals.quantiles = salquantiles;
peakstats.randsals.quantile_values = [.025 .25 .50 .75 .975];

% Find the largest set of peaks that satisfies p<.05
peakstats.setofsigpeaks = [];
peakstats.setofsigpeaklags = [];
if ~isempty(peakpvalues)
    critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
    ptestinds = find([critpvaleachset{:}]<.05);
    if ~isempty(ptestinds)
        % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
        peakstats.setofsigpeaks = peaksets{max(ptestinds)}(~isnan(peaksets{max(ptestinds)}));           
    end;
end;
if ~isempty(peakstats.setofsigpeaks)
   % peakstats.setofsigpeaklags = randbc2plot( isel( peakstats.setofsigpeaks ) );
    peakstats.setofsigpeaklags = randbc2plot( peakstats.setofsigpeaks );
end;

% Find the largest set of peaks that satisfies p<.01
peakstats.setofsigpeaks_p01 = [];
peakstats.setofsigpeaklags_p01 = [];
if ~isempty(peakpvalues)
    critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
    ptestinds = find([critpvaleachset{:}]<.01);
    if ~isempty(ptestinds)
        % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
        peakstats.setofsigpeaks_p01 = peaksets{max(ptestinds)}(~isnan(peaksets{max(ptestinds)}));           
    end;
end;
if ~isempty(peakstats.setofsigpeaks_p01)
    %peakstats.setofsigpeaklags_p01 = randbc2plot( isel( peakstats.setofsigpeaks_p01 ) );
    peakstats.setofsigpeaklags_p01 = randbc2plot( peakstats.setofsigpeaks_p01 );
end;

% Find the largest set of peaks that satisfies p<.001
peakstats.setofsigpeaks_p001 = [];
peakstats.setofsigpeaklags_p001 = [];
if ~isempty(peakpvalues)
    critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
    ptestinds = find([critpvaleachset{:}]<.001);
    if ~isempty(ptestinds)
        % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
        peakstats.setofsigpeaks_p001 = peaksets{max(ptestinds)}(~isnan(peaksets{max(ptestinds)}));           
    end;
end;
if ~isempty(peakstats.setofsigpeaks_p001)
    peakstats.setofsigpeaklags_p001 = randbc2plot( peakstats.setofsigpeaks_p001 );
end;

% Find the largest set of peaks that satisfies p<.001
peakstats.setofsigpeaks_p0001 = [];
peakstats.setofsigpeaklags_p0001 = [];
if ~isempty(peakpvalues)
    critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
    ptestinds = find([critpvaleachset{:}]<.0001);
    if ~isempty(ptestinds)
        % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
        peakstats.setofsigpeaks_p0001 = peaksets{max(ptestinds)}(~isnan(peaksets{max(ptestinds)}));           
    end;
end;
if ~isempty(peakstats.setofsigpeaks_p0001)
    peakstats.setofsigpeaklags_p0001 = randbc2plot( peakstats.setofsigpeaks_p0001 );
end;

% Additional data structure for storing extra infomration about the
% bootstrap.
peakstats_detail.salhist = salhist;
peakstats_detail.salpval = salpval;
peakstats_detail.randpeaksets = randpeaksets;
peakstats_detail.randtrofsets = randtrofsets;
peakstats_detail.randsaliences = randsaliences;
peakstats_detail.allrandsmoothsacs = allrandsmoothsacs;

if plotopt
    subplot(sacploth);
    plot(lagaxis(peakstats.setofsigpeaks),smoothsac(peakstats.setofsigpeaks),'*r'); 
    plot(lagaxis(peakstats.setofsigpeaks_p01),smoothsac(peakstats.setofsigpeaks_p01),'*m'); 
    plot(lagaxis(peakstats.setofsigpeaks_p001),smoothsac(peakstats.setofsigpeaks_p001),'*g'); 

    subplot(salploth);
    xlabel('Salience (peak-trough difference)');
    ylabel('# peaks'); 
    titlestr = ['Salience histograms AM freq:' num2str(amfreq) 'Hz'];
    if nargin>6
        titlestr  = ['Dataset ID:' num2str(unitid) titlestr];
    end;
    title(titlestr);
    set(gcf,'Name',titlestr);
    if ~isempty(saliences)
        xlim([0 1.2*max([saliences{:}])]);
    end;

    if saveplot 
        if nargin<6
           fprint('Cannot save: dataset number not specified.');
        else
            figname = sprintf('N%d_mod%dHz.tiff',unitid,amfreq); 
            %try 
            print('-dtiff','-r200',[plotdir '\' figname]);     
            %catch; end;
        end;
      else
        fprintf('Press a key to continue...');
        pause;
    end;
    close(sacfigh);        
end;



