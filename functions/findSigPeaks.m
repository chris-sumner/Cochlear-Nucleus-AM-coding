function unitpeakstats = findSigPeaks(unitpeakstatsin,pval,minsal,maxperiods)

unitpeakstats = unitpeakstatsin;

% Setup the output strings. 
switch(pval)
    case .05
        peaksfield = ['setofsigpeaks'];
        lagsfield = ['setofsigpeaklags'];
    case .01
        peaksfield = ['setofsigpeaks_p01'];
        lagsfield = ['setofsigpeaklags_p01'];
    case .001
        peaksfield = ['setofsigpeaks_p001'];
        lagsfield = ['setofsigpeaklags_p001'];
    case .0001
        peaksfield = ['setofsigpeaks_p0001'];
        lagsfield = ['setofsigpeaklags_p0001'];
    case .00001
        peaksfield = ['setofsigpeaks_p00001'];
        lagsfield = ['setofsigpeaklags_p00001'];
end;

% Loop around all the AM freqencies. 
for amind = 1:length(unitpeakstatsin)

    peakpvalues = unitpeakstats(amind).peakpvalues;
    % N.B. Due to poor coding peakpvalues can be smaller than peaksets &
    % salsets
    peaksets =  unitpeakstats(amind).peaksets(1:length(peakpvalues));
    salsets  =  unitpeakstats(amind).salsets(1:length(peakpvalues));
    
    % Insitialise to empty.
    unitpeakstats(amind).(peaksfield) = [];
    unitpeakstats(amind).(lagsfield) = [];


    % Find the largest set of peaks that satisfies p
    if ~isempty(peakpvalues)
        % The largest p value for each set.
        critpvaleachset = cellfun(@(x) max(x),peakpvalues,'UniformOutput',false);
        % The minimum salience for each set.
        minsaleachset = cellfun(@(x) min(x),salsets);        
        unitpeakstats(amind).minsals = minsaleachset;
        
        % Find which sets meet that criteria.
        testinds = find([critpvaleachset{:}]<pval & minsaleachset>minsal);
        
        if ~isempty(testinds)
            % N.B. Check to remove nans as the current algortihm will introduce some unpaired troughs. 
            unitpeakstats(amind).(peaksfield) = peaksets{max(testinds)}(~isnan(peaksets{max(testinds)}));           
        end;
    end;

    % Find that set of lags. 
    if ~isempty(unitpeakstats(amind).(peaksfield))
        unitpeakstats(amind).(lagsfield) = unitpeakstats(amind).lagaxis(  unitpeakstats(amind).(peaksfield) ) ;

        % Remove any lags outside of max periods. 
        okperiods = abs(unitpeakstats(amind).(lagsfield)) <= maxperiods*1e3/unitpeakstats(amind).amfreq;
        unitpeakstats(amind).(lagsfield) = unitpeakstats(amind).(lagsfield)(okperiods);
        unitpeakstats(amind).(peaksfield) = unitpeakstats(amind).(peaksfield)(okperiods);
    end;
    
end; 
