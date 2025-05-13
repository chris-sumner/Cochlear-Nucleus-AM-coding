
%% ------------------------------------------------------------------------
%              Reshape this into individual sets for each MTF.

% N.B. uncomment last line to save
clear; 

% Option to make the file as small as possible for github. 
makesmall = true;

% Add the paths required.
addpath functions;

% Loads unit info and other supporting data.
load('rawdata\unitList');

% Load up the CI and VSs . 
load datafiles\SpikeStats;

%load datafiles\SACpeaks; % Too big for github
% shrinkSACpeaks splits it into 3.
% Load up as three parts and reconstitute. 
load datafiles\SACpeaksA;
load datafiles\SACpeaksB;
load datafiles\SACpeaksC;
SACpeaks = {SACpeaksA{:} SACpeaksB{:} SACpeaksC{:}}; 

% indices of lists to keep: modFreq<=2000 - this is most of them anyway.
keepConds = wholeList171212(:,find(strcmp(EXPLOGLIST,'modFreq')))<=2000; 
% Shorten the data to exclude MTFs we will not analyse. 
wholeList171212 = wholeList171212(keepConds,:);
allAMtypes88and91 = allAMtypes88and91(keepConds);

% reduce the size of VS, ML, CI & GC accordingly.
VS = VS(keepConds);
ML  = ML(keepConds);
CI = CI(keepConds);
GC = GC(keepConds);
SACpeaks = SACpeaks(keepConds);

% dividing the data up into individual datasets in units.
[newdatainds dividingdata] = divideData(wholeList171212,EXPLOGLIST,allAMtypes88and91);

% Which data to process. 
datalist = 1:length(newdatainds);
minspkpersweep = 2;     % Minimum spikes per sweep: 2 = 25 spikes per sec for an 80ms window.
                        % Below this VS is unreliable. 
                        
% Loop to restructure the data into MTF dataset groups.
for thisdatastartind = datalist
    
    % The indexes in the wholelist format where the data should exist.
    thisDataWholeListInds = [];
    if thisdatastartind<length(newdatainds)
        endind = newdatainds(thisdatastartind+1)-1;
    else
        endind = size((newdatainds),1);
    end
    thisDataWholeListInds = [newdatainds(thisdatastartind):endind];

    % Which of these conditions have enough spikes in them.
    conditionswithspikes = [];
    if ~isempty( thisDataWholeListInds )
        % Take this from the VS analysis
        vstmp = VS(thisDataWholeListInds);
        inds = find(cellfun(@(x) ~isempty(x),vstmp));
        spikespercondition = zeros(size(vstmp));
        spikespercondition(inds) = cellfun(@(x) x.spikes_per_sweep, vstmp(inds) );
        conditionswithspikes = find(spikespercondition>=minspkpersweep);
    end;        

    if ~isempty(thisDataWholeListInds) && ~isempty(conditionswithspikes)
        
        % --------- This stores the actual VS analysis ---------
        
        % If for github, we remove large fields which are not needed.
        VSstats =  [VS{  thisDataWholeListInds( conditionswithspikes )   }];
        if makesmall
            VSstats = arrayfun(@(x) rmfield(x,{'PH','PHbins'}),VSstats);
        end;
        
        unitVSoutputs(thisdatastartind).VSstats = VSstats; %  ...
            %[VS{  thisDataWholeListInds( conditionswithspikes )   }];

        
        % --------- This stores the ML analysis ----------------
        % Make a temporary copy of all the records for this recording
        % remove a few fields.)
        mltmp = ML( thisDataWholeListInds( conditionswithspikes )   );
        
        if ~isempty(mltmp{1}) 
            % NHPP based Z needs extracting.  
            NHPPind = cellfun(@(x) ~isempty(x.NHPP.fixedDt),mltmp);
            mltmp(NHPPind) = cellfun(@(x) setfield(x,'Z_NHPP',x.NHPP.fixedDt.Z),mltmp(NHPPind),'UniformOutput',false);
            mltmp(~NHPPind) = cellfun(@(x) setfield(x,'Z_NHPP',nan),mltmp(~NHPPind),'UniformOutput',false);

            % Some fields we remove. 
            mlfields2go = {'st','zTerms','ISI','ph','PHsur','Isur','NHPP','EI'};
            mltmp = cellfun(@(x) rmfield(x,mlfields2go),mltmp,'uni',false);

            % For github, also remove zDist. 
            if makesmall
                mltmp = cellfun(@(x) rmfield(x,'zDist'),mltmp,'uni',false);
            end;
                
            unitVSoutputs(thisdatastartind).ML = [mltmp{ : }];
        else
            unitVSoutputs(thisdatastartind).ML = [];
        end;
        
        % --------- This stores the CI analysis ----------------        
        % Make a temporary copy of all the records for this recording
        % remove a few fields.
        cifields2go = {'unsmoothedSAC','lags'};
        citmp = CI( thisDataWholeListInds( conditionswithspikes )   );
    
        if ~isempty(citmp{1}) 
            citmp = cellfun(@(x) rmfield(x,cifields2go),citmp,'uni',false);        
            unitVSoutputs(thisdatastartind).CI = [citmp{ : }];
        else
            unitVSoutputs(thisdatastartind).CI = [];
        end;

        % --------- This stores the SAC peaks analysis --------

        % This is involved becuase we do not pick peak at all the frequencies
        % that other analyses are performed on because the analysis is so slow
        % and sensitive to noisy SACs. 

        sacpeakstmp = SACpeaks( thisDataWholeListInds( conditionswithspikes )   );

        % Which fields we keep. 
        sacpeakfields2keep = {'amfreq','am_period','maxlag_ms','smoothsac', 'lagaxis', ...
             'spikespersweep','setofsigpeaks_p0001' 'setofsigpeaklags_p0001'};
        newfields = {'setofsigpeakhgts_p0001',[], ...
                     'setofsigtrofs_p0001',[], ...
                     'setofsigtrofdepths_p0001', [] ...
                     'setofsigsaliences_p0001', []};
        
        % Makes a list of pairs - parameter, value - for making the
        % structure with.
        sacpeakpars = cell(2*length(sacpeakfields2keep),1);
        sacpeakpars(1:2:end) = sacpeakfields2keep;
        sacpeaksdummystruct = struct(sacpeakpars{:},newfields{:});       % A dummy structure.  

        % Initialise to empty structures of the same fields even if empty.
        unitVSoutputs(thisdatastartind).SACpeaks = cellfun(@(x) sacpeaksdummystruct,  sacpeakstmp  );

        % Loop through each am frequency. 
        for ami = 1:length(sacpeakstmp)
            % If there are some fields in the structure. 
            if isfield(sacpeakstmp{ami},'amfreq')
                               
                % Construct a new stucture. 
                fieldidx = find(ismember(fieldnames(sacpeakstmp{ami}),sacpeakfields2keep)); % Find the fields we want to keep.
                structvals = struct2cell(sacpeakstmp{ami});                                 % Extract the values. 
                sacpeakpars2 = sacpeakpars;                                                 % Copy the pars for call to struct
                sacpeakpars2(2:2:end) = structvals(fieldidx);                               % Find in the values. 
                sacpeakstruct = struct(sacpeakpars2{:},newfields{:});                       % Make the structure. 
                
                % Pull out the peak amplitudes, troughs and saliences for p == 0.0001 
                
                % Find out which of the lists of peaks, troughs and
                % salience we need by looking at the length of setofsigpeaks_p0001
                num_peaks = length(sacpeakstruct.setofsigpeaks_p0001);
                                                
                if num_peaks>0
                    set_lengths  = cellfun(@(x) length(x), sacpeakstmp{ami}.peaksets);
                    which_set = min( find( set_lengths == num_peaks ) );

                    % 2. Get those values for peaksets, trofsets, salsets 
                    sacpeakstruct.setofsigpeaks_p0001 =  sacpeakstmp{ami}.peaksets{which_set};
                    % The next line should not be necessary but there are
                    % some errors from SAC peaks - this relcaluates the delays of those peaks. 
                    sacpeakstruct.setofsigpeaklags_p0001 =  sacpeakstruct.lagaxis(sacpeakstruct.setofsigpeaks_p0001);
                    sacpeakstruct.setofsigpeakhgts_p0001 = sacpeakstruct.smoothsac(sacpeakstruct.setofsigpeaks_p0001);
                    sacpeakstruct.setofsigtrofs_p0001 =  sacpeakstmp{ami}.trofsets{which_set};
                    sacpeakstruct.setofsigtrofdepths_p0001 =  sacpeakstruct.smoothsac(sacpeakstruct.setofsigtrofs_p0001);
                    sacpeakstruct.setofsigsaliences_p0001 =  sacpeakstmp{ami}.salsets{which_set};
                end;
                
                % For github we do not store any of the SACs themselves,
                % UNLESS they are one of the examples used in the figure.  
                if makesmall 
                    % Check for the example unit and the modlation
                    % frequency of the plotted SAC.
                    unitid =  unique(wholeList171212(thisDataWholeListInds( conditionswithspikes ),find(strcmp(EXPLOGLIST,'uniqueID'))));
                    ChSexampleChk = (unitid == 88299021) & (sacpeakstruct.amfreq == 125);
                    PLexampleChk = (unitid == 88340053) & (sacpeakstruct.amfreq == 150);
                    % We set all the others to empty. 
                    if ~ChSexampleChk & ~PLexampleChk 
                        sacpeakstruct.smoothsac = [];
                        sacpeakstruct.lagaxis = [];
                    end;
                end;
                
                % Store this sturcture in place of the dummy empty one. 
                unitVSoutputs(thisdatastartind).SACpeaks(ami) = sacpeakstruct; 

            end;
        end;

        % --------- This stores the actual GC analysis ---------
                
        GCtmp = [GC{  thisDataWholeListInds( conditionswithspikes )   }];        

        % For github    
        if makesmall
            GCtmp = arrayfun(@(x) rmfield(x,'splitpsths'),GCtmp);
        end;
      
        unitVSoutputs(thisdatastartind).GC = GCtmp; %...
            %[GC{  thisDataWholeListInds( conditionswithspikes )   }];

    else
        unitVSoutputs(thisdatastartind).VSstats = [];
        unitVSoutputs(thisdatastartind).ML = [];
        unitVSoutputs(thisdatastartind).CI = [];
        unitVSoutputs(thisdatastartind).GC = [];
        unitVSoutputs(thisdatastartind).SACpeaks = [];
    end;
                
    % make the unitinfo.
    unitVSoutputs(thisdatastartind).unitinfo = ...
        constructUnitInfo(wholeList171212,EXPLOGLIST,typeNames,thisdatastartind, ...
         thisDataWholeListInds,newdatainds,conditionswithspikes);
        
    fprintf('.');  
    if rem(thisdatastartind,40) == 0; fprintf('\n'); end;
end;

% Save this to a file.
%save datafiles\SpikeStatsSets_small unitVSoutputs;