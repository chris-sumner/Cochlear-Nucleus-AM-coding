function [unitpeakstats  thisdatainds conditionswithspikes] = findSACpeaksPar(thisdatastartind,unitid,maxlagfactor,wfr,analwindow,figpath,newdatainds, minspkpersweep, wholeList171212, allSpikeTimes88and91, EXPLOGLIST)
% findSACpeaks  uses a peak picking algorithm to find peaks in the SAC. 
% N.B. several global variables used

thresholdpvalue = .05;
tstart = analwindow(1);
tend = min(analwindow(2), wholeList171212(thisdatastartind,find(strcmp(EXPLOGLIST,'amToneDur')))) ;

% ------------ Get the spike conditions for analysis ------------

    % Find the rows this piece of data occupy.
    if thisdatastartind<length(newdatainds)
        endind = newdatainds(thisdatastartind+1)-1;
    else
        endind = size((newdatainds),1);
    end
    thisdatainds = [newdatainds(thisdatastartind):endind];
    
    if isempty( thisdatainds )
        unitpeakstats = [];
        thisdatainds = [];
        conditionswithspikes = [];
        return;
    end;        
    
    % Get all the spike trains.
    spktrainset = [ allSpikeTimes88and91{thisdatainds} ];
    % Select out the analysis window only (and subtract start time from spike times)
    spktrainset = cellfun(@(x)  x(x>tstart & x<=tend)-tstart,spktrainset,'UniformOutput',false);
    
    allamfreqs = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'modFreq')));
    
    % Remove any conditions with a mean of a criterions spikes per sweep.
    spktrainset_chk = cellfun(@(x) (~isempty(x) && length(x)>minspkpersweep), spktrainset);
    spikespersweep = mean(cellfun(@(x) length(x), spktrainset));
    conditionswithspikes = find(sum(spktrainset_chk,1));
    
    % If there is no data return empty  
    if isempty(conditionswithspikes)
        unitpeakstats = [];
        thisdatainds = [];
        conditionswithspikes = [];
        return;
    end;
    
    spktrainset = spktrainset(:,conditionswithspikes);
    spikespersweep = mean(cellfun(@(x) length(x), spktrainset));
    amfreqset = ones(size(spktrainset,1),1)*(allamfreqs(conditionswithspikes)');
    amfreqs = amfreqset(1,:)';


% ------------- Otherwise contiue to do the peak picking --------

% Cutdown to 20-100ms.
duration = diff(analwindow);
randlagfactor = 2;

fprintf('findSACpeaks: analysing unit %d.\n',unitid);
fprintf('AM frequencies:\n');  fprintf('\t%dHz\n ',amfreqs);

% Make a directory for that unit if it does not exist. 
if nargin>6
    unitfigpath = [figpath '\' num2str(unitid) '_N' num2str(thisdatastartind) ];
    lookfordir = dir(unitfigpath);
    if isempty(lookfordir)
        mkdir(unitfigpath);
    end;
end;
   

for amind = 1:length(amfreqs)
    
    fprintf('AM frequency: %dHz\n',amfreqs(amind));

    sacfigh = figure; set(gcf,'Position',[200 200 600 400]);
    maxlag_ms= maxlagfactor*1e3/amfreqs(amind);
    
    spikesets_1am = spktrainset(:,amind);
    % N.B. maxlag + 10% solves a end effects problem. 
    [h2plot bc2plot] = SPTCORR(spikesets_1am,'nodiag',1.1*maxlag_ms,.05,duration,'LouageNorm');

    wlen = round(1e3*wfr/(amfreqs(amind)*(bc2plot(2)-bc2plot(1))))+1;   % Smoothing window.
    isel = find(abs(bc2plot)<maxlag_ms);                                 % Analysis window.
    
    % Smooth with a symmetrical n-point moving average. 
    if wlen>1
        smoothsac = filtfilt(ones(1,wlen)/wlen,1,h2plot);
    else
        smoothsac = h2plot;
    end;
    
    % Plot the results
    sacploth = subplot(2,1,1);
    plot(bc2plot(isel),smoothsac(isel),'color',[.5 .5 1]); hold on;
    xlabel('Delay (ms)');
    titlestr = ['UNIT:' num2str(unitid) ' SAC for AM freq:' num2str(amfreqs(amind)) 'Hz'];
    ylabel('CI'); title(titlestr); set(gcf,'Name',titlestr);
    am_period = 1e3/amfreqs(amind);
    line([am_period am_period],[0 max(ylim)]);
    
    % This works well, but needs a good way to set the criterion. 
    [peaksets, trofsets salsets] = findmainpeaks5(smoothsac(isel));

    if amfreqs(amind) == 850
        fprintf('!');
    end;
        
    % -- BOOTSTRAPPING TO GET RANDOM RESULTS ----
    randn = 500;
    fprintf('Computing bootstrapped distribtutions of peak saliences:\n');
    %clear  rand_spikesets randpeaksets randtrofsets randsaliences
    
    randpeaksets = cell(randn,1);
    randtrofsets = cell(randn,1);
    randsaliences = cell(randn,1);
    allrandsmoothsacs = cell(randn,1);
    
    parfor ri = 1:randn
        
        % Define variables
        rand_spikesets = cell(length(spikesets_1am),1);
        isis = cell(length(spikesets_1am),1);
         
        % --- Randomise intervals but across all trials ---
        
        % Accumulate all intervals. 
        for sweepi = 1:length(spikesets_1am)                      
            % This is another way - maintain the ISIs. 
            isis{sweepi} =  diff(spikesets_1am{sweepi});
        end;
        % Convert to array and then shuffle.
        isiarray = [isis{:}];
        surrograteisis = isiarray(randperm(length(isiarray)));
        % Now convert back to spike times - same # spikes in each sweep. 
        counter = 1;
        for sweepi = 1:length(spikesets_1am)                      
            inc = length(spikesets_1am{sweepi})-1;
            if length(spikesets_1am{sweepi})==1
                rand_spikesets{sweepi} =  [spikesets_1am{sweepi}(1)];
            elseif length(spikesets_1am{sweepi})>1
                rand_spikesets{sweepi} =  cumsum([spikesets_1am{sweepi}(1) surrograteisis(counter:counter+inc-1)]);
                counter = counter+inc;
            end;
        end;        
        
        
        % Compute the SAC.
        % N.B. maxlag + 10% solves a end effects problem. 
        % This is the main delay in running - which has to be within the
        % parfor loop. 
        [randh2plot randbc2plot] = SPTCORR(rand_spikesets,'nodiag',1.1*maxlag_ms,.05,duration,'LouageNorm');

        % Smooth the SAC function.
        if wlen>1
            randsmoothsac = filtfilt(ones(1,wlen)/wlen,1,randh2plot);
        else
            randsmoothsac = randh2plot;
        end;
        % Store all the smoothed sacs.
        allrandsmoothsacs{ri} = randsmoothsac;
        
        fprintf('.'); if rem(ri,50)==0 fprintf('\n'); end;
    end;    
    
   % findmainpeaks5 will not work within the parfor loop. However this does not slow things down significantly.  
   for ri = 1:randn
        % Find the peaks and troughs in this random set. 
        [randpeaksets{ri}, randtrofsets{ri}, randsaliences{ri}] = findmainpeaks5(allrandsmoothsacs{ri}(isel));
        
        % Plot example random set.
        if ri<10;
            %figure(sacfigh);
            subplot(sacploth);
            plot(bc2plot(isel),allrandsmoothsacs{ri}(isel),'color',[1 .7 .7]);
        end;
   end;
   % Store one example.
   randsmoothsaceg = allrandsmoothsacs{1};

    
    fprintf('\nDone.\n');

    % Statistics of the saliences (peaks and troughs) found in the random
    % data. 
    
    salvals = cell(120,1);
    saltmp = [randsaliences{:}];
    for i=1:length(saltmp)
        sallen =  length(saltmp{i});    
        
        % One way of computing the peaks and troughs:
        %salvals{sallen}  = [ salvals{sallen} saltmp{i} ];    
        % salval{1} is the saliences of all the peaks and troughs in the random data
        % if you only pick one peak.      
        % salval{2} is the saliences of all the peaks and troughs in the
        % random data when you pick two peaks - so by definition the largest
        % and next largest. 

        salvals{sallen}  = [ salvals{sallen} min(saltmp{i}) ];    
        % salval{2} is the MINIMUM salience of all the peaks and troughs in the
        % random data when you pick two peaks - so by definition the second largest
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

    % Make a figure of the salience distributions.
    %salfigh = figure; set(gcf,'Position',[200 200 600 400]);
    salploth = subplot(2,1,2);
    colorspecs = get(gca,'ColorOrder');
    plot(salbins,salhist);
    
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
        
        % Plot these on the salience figure. 
        %figure(salfigh);    
        subplot(salploth);
        hold on;
        plot(saliences{setlen},pi*15*ones(1,length(saliences{setlen})), ...
            'o','color',colorspecs(rem(setlen-1,7)+1,:),'markerface',colorspecs(rem(setlen-1,7)+1,:));
        maxx =  max(saliences{setlen});
        if ~isempty(maxx)
            text(maxx*1.1,pi*50,sprintf('%d:%.2g',setlen,max(peakpvalues{setlen})),'fontsize',6);
        end;
        
        % Plot the peaks and troughs, shading the grey so that the darkest
        % are more 'primary'
        subplot(sacploth);
        hold on;
        %         plot((bc2plot(isel(peaksets{pi}))),smoothsac(isel(peaksets{pi})), ...
        %              'o','color',colorspecs(rem(setlen-1,7)+1,:),'markerface',colorspecs(rem(setlen-1,7)+1,:));
        %                hold on;
        inds = peaksets{pi}(~isnan(peaksets{pi}));
        plot((bc2plot(isel(inds))),smoothsac(isel(inds)), ...
             'o','color',colorspecs(rem(setlen-1,7)+1,:),'markerface',colorspecs(rem(setlen-1,7)+1,:));

    end;
    
    % Finalise salience plot.
    %figure(salfigh);
    subplot(salploth);
    xlabel('Salience (peak-trough difference)');
    ylabel('# peaks'); 
    titlestr = ['Salience histograms AM freq:' num2str(amfreqs(amind)) 'Hz'];
    title(titlestr);
    set(gcf,'Name',titlestr);
    if ~isempty(saliences)
        xlim([0 1.2*max([saliences{:}])]);
    end;
    
    % Put all this information into a structure. 
    unitpeakstats(amind).wlen =wlen;
    unitpeakstats(amind).amfreq = amfreqs(amind);
    unitpeakstats(amind).am_period = am_period;
    unitpeakstats(amind).peaksets =peaksets;
    unitpeakstats(amind).trofsets = trofsets;
    unitpeakstats(amind).salsets =salsets;
    unitpeakstats(amind).maxlag_ms = maxlag_ms;
    unitpeakstats(amind).smoothsac =smoothsac(isel);
    unitpeakstats(amind).randsmoothsaceg =randsmoothsaceg;
    unitpeakstats(amind).peakpvalues = peakpvalues;
    unitpeakstats(amind).lagaxis = bc2plot(isel);
    unitpeakstats(amind).spikespersweep = spikespersweep;

    % Statistics abou the saliences of the random data.
    unitpeakstats(amind).randsals.mean = salmean;
    unitpeakstats(amind).randsals.median = salmedian;
    unitpeakstats(amind).randsals.sd = salsd;
    unitpeakstats(amind).randsals.quantiles = salquantiles;
    unitpeakstats(amind).randsals.quantile_values = [.025 .25 .50 .75 .975];
    
    % Find the largest set of peaks that satisfies p<.05
    unitpeakstats(amind).setofsigpeaks = [];
    unitpeakstats(amind).setofsigpeaklags = [];
    if ~isempty(peakpvalues)
        critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
        ptestinds = find([critpvaleachset{:}]<.05);
        if ~isempty(ptestinds)
            % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
            unitpeakstats(amind).setofsigpeaks = peaksets{max(ptestinds)}(~isnan(peaksets{max(ptestinds)}));           
        end;
    end;
    if ~isempty(unitpeakstats(amind).setofsigpeaks)
        unitpeakstats(amind).setofsigpeaklags = bc2plot( isel( unitpeakstats(amind).setofsigpeaks ) );
    end;

    % Find the largest set of peaks that satisfies p<.01
    unitpeakstats(amind).setofsigpeaks_p01 = [];
    unitpeakstats(amind).setofsigpeaklags_p01 = [];
    if ~isempty(peakpvalues)
        critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
        ptestinds = find([critpvaleachset{:}]<.01);
        if ~isempty(ptestinds)
            % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
            unitpeakstats(amind).setofsigpeaks_p01 = peaksets{max(ptestinds)}(~isnan(peaksets{max(ptestinds)}));           
        end;
    end;
    if ~isempty(unitpeakstats(amind).setofsigpeaks_p01)
        unitpeakstats(amind).setofsigpeaklags_p01 = bc2plot( isel( unitpeakstats(amind).setofsigpeaks_p01 ) );
    end;
    
    % Find the largest set of peaks that satisfies p<.001
    unitpeakstats(amind).setofsigpeaks_p001 = [];
    unitpeakstats(amind).setofsigpeaklags_p001 = [];
    if ~isempty(peakpvalues)
        critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
        ptestinds = find([critpvaleachset{:}]<.001);
        if ~isempty(ptestinds)
            % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
            unitpeakstats(amind).setofsigpeaks_p001 = peaksets{max(ptestinds)}(~isnan(peaksets{max(ptestinds)}));           
        end;
    end;
    if ~isempty(unitpeakstats(amind).setofsigpeaks_p001)
        unitpeakstats(amind).setofsigpeaklags_p001 = bc2plot( isel( unitpeakstats(amind).setofsigpeaks_p001 ) );
    end;
    
    % Find the largest set of peaks that satisfies p<.001
    unitpeakstats(amind).setofsigpeaks_p0001 = [];
    unitpeakstats(amind).setofsigpeaklags_p0001 = [];
    if ~isempty(peakpvalues)
        critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
        ptestinds = find([critpvaleachset{:}]<.0001);
        if ~isempty(ptestinds)
            % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
            unitpeakstats(amind).setofsigpeaks_p0001 = peaksets{max(ptestinds)}(~isnan(peaksets{max(ptestinds)}));           
        end;
    end;
    if ~isempty(unitpeakstats(amind).setofsigpeaks_p0001)
        unitpeakstats(amind).setofsigpeaklags_p0001 = bc2plot( isel( unitpeakstats(amind).setofsigpeaks_p0001 ) );
    end;
    
    % Additional data structure for storing extra infomration about the
    % bootstrap.
    unitpeakstats_detail(amind).salhist = salhist;
    unitpeakstats_detail(amind).salpval = salpval;
    unitpeakstats_detail(amind).randpeaksets = randpeaksets;
    unitpeakstats_detail(amind).randtrofsets = randtrofsets;
    unitpeakstats_detail(amind).randsaliences = randsaliences;
    unitpeakstats_detail(amind).allrandsmoothsacs = allrandsmoothsacs;
    
    % PLot the figure if a path is provided. 
    if nargin>6
        figname = sprintf('U%d_mod%dHz.tiff',unitid,amfreqs(amind)); 
        try 
        print('-dtiff','-r200',[unitfigpath '\' figname]);     
        catch; end;
        close(sacfigh);        
    end;
    
    pause;
end;

% Save extra info about the dataset here. 
dataname = sprintf('U%d_N%d',unitid,thisdatastartind); 
save([unitfigpath '\' dataname],'unitpeakstats','unitpeakstats_detail','thisdatainds','conditionswithspikes');

