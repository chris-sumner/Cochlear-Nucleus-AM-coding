function statsmat = collate_byDataset(unitCIandVS,unitWR,varargin)
% collate_byDataset brings together two separate sets of analysis.

unitWR = getBestWROutputs(unitWR);      % Select best tau for each set. 
numDatasets = min([length(unitCIandVS), length(unitWR)] );

fprintf('Collating data...');
fprintf('\nPercent complete:');

% ----------------------  Load up  data  ------------------------
% This is annoying - in order to confirm the AM type - we need to do the
% correspondence checking all over again. 
% Should have done this at an earlier stage. 

load('rawdata\unitList');
colfn =  @(s) find(strcmp(EXPLOGLIST,s));

% Can only use this list if you first remove any frequencies>2000
keepConds = wholeList171212(:,find(strcmp(EXPLOGLIST,'modFreq')))<=2000; 
% Shorten the data to exclude MTFs we will not analyse. 
wholeList171212 = wholeList171212(keepConds,:);
allAMtypes88and91 = allAMtypes88and91(keepConds);

pcdone = 0;

% Loop through every data set.
for i=1:numDatasets
    
    % Local copy of the unitinfo for this data set for convenience.
    thisunitinfo = unitCIandVS(i).unitinfo;
    

    % --------------- retrieve the AM type ----------------
    datasetChk = 1;
    if ~isempty( thisunitinfo.wholeListInds )
        AMType = unique( allAMtypes88and91( thisunitinfo.wholeListInds ) );

         % Store: give warning if ambiguous.
         if length(AMType)>1
             error('AMType is ambiguous for this dataset: %d\n',i);
             AMType = 'BAD';
         end;

        % Double check the unit info against to the conditions in whole unit list.
        % Extremely unlikely to be wrong but...

        % The Scholes' wholeList for this dataset.
        dataSetScholesList = wholeList171212(thisunitinfo.wholeListInds,:); 

        % Check all these things agains the unitinfo.
        % N.B. Special treatment of carrier frequency in case of nans (ie.
        % noise carrier, cf and cf_thr can all be nan)
        datasetChk = prod( ...
                dataSetScholesList(:,colfn('Unit')) == rem(thisunitinfo.unitid,1000)  & ...
                dataSetScholesList(:,colfn('Expt')) == round(thisunitinfo.unitid/1000)  & ...
                ( dataSetScholesList(:,colfn('cf')) == thisunitinfo.cf | ...
                      isnan(dataSetScholesList(:,colfn('cf'))) & isnan(thisunitinfo.cf) ) & ...
                ( dataSetScholesList(:,colfn('cf_thr')) == thisunitinfo.cf_thr | ...
                      isnan(dataSetScholesList(:,colfn('cf_thr')))  & isnan( thisunitinfo.cf_thr ) ) &...
                ( dataSetScholesList(:,colfn('carrFreq')) == thisunitinfo.carrFreq | ...
                    isnan(dataSetScholesList(:,colfn('carrFreq'))) == isnan( thisunitinfo.carrFreq ) ) & ...
                dataSetScholesList(:,colfn('amToneDur')) == thisunitinfo.amToneDur & ...
                dataSetScholesList(:,colfn('modLevel')) == thisunitinfo.modLevel & ...
                dataSetScholesList(:,colfn('depthMod')) == thisunitinfo.depthMod ...
                );

        % If this fails then the AM type cannot be correct.    
        if ~datasetChk
            error('Failure of correspondence of unitinfo with wholeList171212 for dataset %d\n',i);
        end;
    end;     
     
    
    % --------------------- Basic statistics ----------------------
     
    statsmat.unitid{i} =  unitCIandVS(i).unitinfo.unitid;
    statsmat.modLevel{i} =  unitCIandVS(i).unitinfo.modLevel;
    statsmat.amToneDur{i} =  unitCIandVS(i).unitinfo.amToneDur;
    statsmat.carrFreq{i} =  unitCIandVS(i).unitinfo.carrFreq;
    statsmat.depthMod{i} =  unitCIandVS(i).unitinfo.depthMod;
    statsmat.cf{i} =  unitCIandVS(i).unitinfo.cf;
    statsmat.cf_thr{i} =  unitCIandVS(i).unitinfo.cf_thr;
    statsmat.type{i} =  unitCIandVS(i).unitinfo.TypeName;
    statsmat.typeNum{i} =  unitCIandVS(i).unitinfo.TypeNum;
    statsmat.CV{i} =  unitCIandVS(i).unitinfo.CV;    
    statsmat.allamfreqs{i} = unitCIandVS(i).unitinfo.allamfreqs;
    statsmat.usedamfreqs{i} =  unitCIandVS(i).unitinfo.usedamfreqs;   
    % The frequencies used for the classifier are different. 
    statsmat.usedamfreqs_dprime{i} =  unitWR(i).unitinfo.usedamfreqs;    
    statsmat.amType{i} = AMType;
    statsmat.dataChk{i} = datasetChk;
    % This is a useful index to the raw data which we carry through any
    % reductions in the statsmat size. 
    statsmat.fullDataUnitIndex{i} = i;  
    
    % "rationalised" unit type which groups onset types.. 
    if strcmp(statsmat.type{i},'OnI') || strcmp(statsmat.type{i},'OnC') || strcmp(statsmat.type{i},'OnL')
        statsmat.rationalisedType{i} = 'On';
    % blank out unusual or unclassified.     
    elseif strcmp(statsmat.type{i},'UNC') || strcmp(statsmat.type{i},'Unu') 
        statsmat.rationalisedType{i} = '';
    else    
        statsmat.rationalisedType{i} = statsmat.type{i};
    end;
    
    % "rationalised" signal levels. 
    % Analysis shows that the vast majority of the data are at 30,50
    % and 70dB SPL. So we rationalise to the nearest of those, and
    % round up for the inbetween (40, 60).
    rLvls = [30 50 70];
    statsmat.modLvl_rationalised{i} = max(rLvls(findnearest(statsmat.modLevel{i},rLvls)));
    
    % Rationalised AM frequencies - rounded to the next lowest x50Hz.
    % This may seem odd but most of the dataset are in 50,150, etc. 100Hz steps. 
    statsmat.allamfreqs_rationalised{i} = rationaliseModFreq(statsmat.allamfreqs{i});   
    statsmat.usedamfreqs_rationalised{i} = rationaliseModFreq(statsmat.usedamfreqs{i});
    statsmat.usedamfreqs_dprime_rationalised{i} = rationaliseModFreq(statsmat.usedamfreqs_dprime{i}); 
    
    
    % -------------------- Vector strength ------------------------
    
    statsmat.VS{i} =  arrayfun(@(x) x.VSvalue,unitCIandVS(i).VSstats);
    statsmat.Rayleigh{i} =  arrayfun(@(x) x.Rayleigh,unitCIandVS(i).VSstats);   
        
    % --- Express modulation frequency relative to peak vector strength --
    
    % Find the significant VS and frequencies.
    sigvs_inds = find(statsmat.Rayleigh{i}>13.8);
    sigvs_bmfs = statsmat.usedamfreqs{i}(sigvs_inds);
    sigvs =  statsmat.VS{i}(sigvs_inds);
            
    % Find the BMF.
    [best_vs bmf_ind] = nanmax(sigvs);
    bmf = sigvs_bmfs(bmf_ind);
    bmf_rationalised = rationaliseModFreq(bmf); % 100*ceil(bmf/100);
    if ~isempty(bmf)
        % Calculate modulation frequencies re: BMF for each grouping. 
        statsmat.BMF{i} = bmf;
        statsmat.VS_max{i} = best_vs;
        
        below1kInds = find(statsmat.usedamfreqs{i}<1000 & statsmat.Rayleigh{i}'>13.8);        
        statsmat.VS_meanBelow1k{i} = nanmean(statsmat.VS{i}(below1kInds));
        closeto150Inds = find( (statsmat.usedamfreqs{i}==150 | statsmat.usedamfreqs{i}==200) & statsmat.Rayleigh{i}'>13.8);      
        if ~isempty(closeto150Inds)
            statsmat.VS_150ish{i} = statsmat.VS{i}(closeto150Inds(1));
        else
            statsmat.VS_150ish{i} = nan;
        end;            

        % Calculate the dprimes relative to the best value.
        statsmat.VS_reBestVS{i} = statsmat.VS{i} - best_vs;
        statsmat.VS_pcBestVS{i} = statsmat.VS{i} / best_vs;

        % MTF shape analysis.
        vstmp = statsmat.VS{i};
        vstmp(statsmat.Rayleigh{i}<13.8) = 0.1;  % Replace any non-significant values with .1
        
        % "Corner" is -3dB, cutoff is -dB re: max.  from Sayles. 
        VS_absthr = 0.2;
        tmp = MTFShapeAnal('VS',vstmp,statsmat.usedamfreqs{i},VS_absthr);
        statsmat.VS_LoCornerHz{i} = tmp.MTFLoCorner;
        statsmat.VS_HiCornerHz{i} = tmp.MTFHiCorner;
        statsmat.VS_LoCutOffHz{i} = tmp.MTFLoCutOff;
        statsmat.VS_HiCutOffHz{i} = tmp.MTFHiCutOff;
        statsmat.VS_MTFType{i} = tmp.MTFType;
        statsmat.VS_MTFCornerGoodness{i} = tmp.MTFCornerGoodness;
        statsmat.VS_MTFCutOffGoodness{i} = tmp.MTFCutoffGoodness;
        statsmat.VS_MTFCornerBW{i} = tmp.MTFCornerBW;
        statsmat.VS_MTFCutOffBW{i} = tmp.MTFCutOffBW;
        statsmat.VS_HiFAbsThr{i} = tmp.MTFHiFAbsThr;
        statsmat.VS_LoFAbsThr{i} = tmp.MTFLoFAbsThr;
        statsmat.VS_MTFAbsThrGoodness{i} = tmp.MTFAbsThrGoodness;
        statsmat.VS_MTFAbsThrBW{i} = tmp.MTFAbsThrBW;
        
        % "Used" AM frequencies.
        amfs_rebmf = statsmat.usedamfreqs{i} - bmf;
        statsmat.usedamfreqs_reBMF{i} = amfs_rebmf;                               % store the AMF as relative to BMF.
        % "Used" AM frequencies rationalised.
        % N.B. These are relative not absolute, so are not reationalised in
        % the same way as absolute values. 
        statsmat.usedamfreqs_reBMF_rationalised{i} = 100*round(amfs_rebmf/100);   % store the AMF as relative to BMF.
 
        % "all" AM frequencies.
        amfs_rebmf = statsmat.allamfreqs{i} - bmf;
        statsmat.allamfreqs_reBMF{i} = amfs_rebmf;                               %  store the AMF as relative to BMF.
        % "all" AM frequencies rationalised.
        statsmat.allamfreqs_reBMF_rationalised{i} = 100*round(amfs_rebmf/100);   %  store the AMF as relative to BMF.

        % "used" AM frequencies for the dprime measurements.
        amfs_rebmf = statsmat.usedamfreqs_dprime{i} - bmf;
        statsmat.usedamfreqs_dprime_reBMF{i} = amfs_rebmf;                               %  store the AMF as relative to BMF.
        % "used" AM frequencies for the dprime measurements rationalised.
        statsmat.usedamfreqs_dprime_reBMF_rationalised{i} = 100*round(amfs_rebmf/100);   %  store the AMF as relative to BMF.
        
    else
        statsmat.BMF{i} = nan;
        statsmat.VS_max{i} = nan;
        statsmat.VS_meanBelow1k{i} = nan;
        statsmat.VS_150ish{i} = nan;
        statsmat.VS_reBestVS{i} = nan*statsmat.usedamfreqs{i}; 
        statsmat.VS_pcBestVS{i} = nan*statsmat.usedamfreqs{i}; 
        statsmat.usedamfreqs_reBMF{i} = nan*statsmat.usedamfreqs{i};     % Nans if no significant phase locking
        statsmat.allamfreqs_reBMF{i} = nan*statsmat.allamfreqs{i};     
        statsmat.usedamfreqs_dprime_reBMF{i} = nan*statsmat.usedamfreqs_dprime{i};          
        statsmat.usedamfreqs_reBMF_rationalised{i} = nan*statsmat.usedamfreqs{i};     % Nans if no significant phase locking
        statsmat.allamfreqs_reBMF_rationalised{i} = nan*statsmat.allamfreqs{i};     
        statsmat.usedamfreqs_dprime_reBMF_rationalised{i} = nan*statsmat.usedamfreqs_dprime{i};          
        statsmat.VS_LoCornerHz{i} = nan;
        statsmat.VS_HiCornerHz{i} = nan;
        statsmat.VS_LoCutOffHz{i} = nan;
        statsmat.VS_HiCutOffHz{i} = nan;
        statsmat.VS_MTFType{i} = 'none';
        statsmat.VS_MTFCornerGoodness{i} = nan;
        statsmat.VS_MTFCutOffGoodness{i} = nan;
        statsmat.VS_MTFCornerBW{i} = nan;
        statsmat.VS_MTFCutOffBW{i} = nan;
        statsmat.VS_HiFAbsThr{i} = nan;
        statsmat.VS_LoFAbsThr{i} = nan;
        statsmat.VS_MTFAbsThrGoodness{i} = nan;
        statsmat.VS_MTFAbsThrBW{i} = nan;
    end;
                        
    
   % ------------------ WR classifier results --------------------
   
   % Initialise dprime statistics to nans. etc.
    statsmat.BMFdp{i} = nan;
    statsmat.dprimes_max{i} = nan;
    statsmat.dprimes_reBestDP{i} = nan*statsmat.usedamfreqs_dprime{i};
    statsmat.dprimes_pcBestDP{i} = nan*statsmat.usedamfreqs_dprime{i};
    statsmat.usedamfreqs_reBMFdp{i} = nan*statsmat.usedamfreqs{i};     % Nans if no significant phase locking
    statsmat.allamfreqs_reBMFdp{i} = nan*statsmat.allamfreqs{i};     
    statsmat.usedamfreqs_dprime_reBMFdp{i} = nan*statsmat.usedamfreqs_dprime{i};          
    statsmat.usedamfreqs_reBMFdp_rationalised{i} = nan*statsmat.usedamfreqs{i};     % Nans if no significant phase locking
    statsmat.allamfreqs_reBMFdp_rationalised{i} = nan*statsmat.allamfreqs{i};     
    statsmat.usedamfreqs_dprime_reBMFdp_rationalised{i} = nan*statsmat.usedamfreqs_dprime{i};          
    statsmat.dprime_LoCornerHz{i} = nan;
    statsmat.dprime_HiCornerHz{i} = nan;
    statsmat.dprime_LoCutOffHz{i} = nan;
    statsmat.dprime_HiCutOffHz{i} = nan;
    statsmat.dprime_MTFType{i}    = 'none';
    statsmat.dprime_MTFCornerGoodness{i} = nan;
    statsmat.dprime_MTFCutOffGoodness{i} = nan;
    statsmat.dprime_MTFCornerBW{i} = nan;
    statsmat.dprime_MTFCutOffBW{i} = nan;
    statsmat.dprime_HiFAbsThr{i} = nan;
    statsmat.dprime_LoFAbsThr{i} = nan;
    statsmat.dprime_MTFAbsThrGoodness{i} = nan;
    statsmat.dprime_MTFAbsThrBW{i} = nan;
   
   
   if ~isempty(unitWR(i).bestClassifier)
         statsmat.dprimes{i} = unitWR(i).wroutput.dprime;

         % Find the BMF for the d-prime values.
        [best_dp bmfdp_ind] = nanmax(statsmat.dprimes{i});
        bmfdp = statsmat.usedamfreqs_dprime{i}(bmfdp_ind);
        bmfdp_rationalised = rationaliseModFreq(bmfdp); % 100*ceil(bmfdp/100);
        if ~isempty(bmfdp)
            % Calculate modulation frequencies re: BMF for each grouping. 
            
            % Calculate the dprimes relative to the best value.
            statsmat.dprimes_reBestDP{i} = statsmat.dprimes{i} - best_dp;
            statsmat.dprimes_pcBestDP{i} = statsmat.dprimes{i}/best_dp;
            statsmat.dprimes_max{i} = best_dp;

            % MTF shape analysis.
            dprime_absthr = 1;
            tmp = MTFShapeAnal('DP',statsmat.dprimes{i},statsmat.usedamfreqs_dprime{i},dprime_absthr);
            statsmat.dprime_LoCornerHz{i} = tmp.MTFLoCorner;
            statsmat.dprime_HiCornerHz{i} = tmp.MTFHiCorner;
            statsmat.dprime_LoCutOffHz{i} = tmp.MTFLoCutOff;
            statsmat.dprime_HiCutOffHz{i} = tmp.MTFHiCutOff;
            statsmat.dprime_MTFType{i}     = tmp.MTFType;
            statsmat.dprime_MTFCornerGoodness{i} = tmp.MTFCornerGoodness;
            statsmat.dprime_MTFCutOffGoodness{i} = tmp.MTFCutoffGoodness;
            statsmat.dprime_MTFCornerBW{i} = tmp.MTFCornerBW;
            statsmat.dprime_MTFCutOffBW{i} = tmp.MTFCutOffBW;
            statsmat.dprime_HiFAbsThr{i} = tmp.MTFHiFAbsThr;
            statsmat.dprime_LoFAbsThr{i} = tmp.MTFLoFAbsThr;
            statsmat.dprime_MTFAbsThrGoodness{i} = tmp.MTFAbsThrGoodness;
            statsmat.dprime_MTFAbsThrBW{i} = tmp.MTFAbsThrBW;
                  
            % Store the  BMF-dP
            statsmat.BMFdp{i} = bmfdp;

            % "Used" AM frequencies.
            amfs_rebmfdp = statsmat.usedamfreqs{i} - bmfdp;
            statsmat.usedamfreqs_reBMFdp{i} = amfs_rebmfdp;                               % store the AMF as relative to BMF.
            % "Used" AM frequencies rationalised.
            statsmat.usedamfreqs_reBMFdp_rationalised{i} = 100*round(amfs_rebmfdp/100);   %  store the AMF as relative to BMF.
  
            % "all" AM frequencies.
            amfs_rebmfdp = statsmat.allamfreqs{i} - bmfdp;
            statsmat.allamfreqs_reBMFdp{i} = amfs_rebmfdp;                               %  store the AMF as relative to BMF.
            % "all" AM frequencies rationalised.
            statsmat.allamfreqs_reBMFdp_rationalised{i} = 100*round(amfs_rebmfdp/100);   %  store the AMF as relative to BMF.

            % "used" AM frequencies for the dprime measurements.
            amfs_rebmfdp = statsmat.usedamfreqs_dprime{i} - bmfdp;
            statsmat.usedamfreqs_dprime_reBMFdp{i} = amfs_rebmfdp;                               %  store the AMF as relative to BMF.
            % "used" AM frequencies for the dprime measurements rationalised.
            statsmat.usedamfreqs_dprime_reBMFdp_rationalised{i} = 100*round(amfs_rebmfdp/100);   %  store the AMF as relative to BMF.

        end;
    
           if isempty(statsmat.dprime_MTFCornerGoodness)
                error('Opps');
            end

        
    end;   
    pcdone = round(100*i/numDatasets);
    if rem(i,100)==0
        fprintf('%d,',pcdone);
    end;
    
end;


fprintf('\n%d datasets failed data check\n',sum( ~[statsmat.dataChk{:}] ) );


  


function  tmp = MTFShapeAnal(opt,vals,fs,absthr)
% Analysis the shape of the MTF to extract statistics.

% Normalise the response. 
valnorm = vals/max(vals);

% Find the points that exceed the 70.8% point (-3dB point).
% "Corner frequencies"
thr3dB = 10^(-3/20);

% Low side corner.
abovethr3dB = valnorm>thr3dB;
indL = min(find(abovethr3dB));
if ~isempty(indL) & indL>1 
    tmp.MTFLoCorner = fs(indL);
    % Linearly interpolate to find a cut-off frequency. 
    tmp.MTFLoCorner = fs(indL-1) + ...
        (fs(indL)-fs(indL-1))*(thr3dB-valnorm(indL-1))/(valnorm(indL)-valnorm(indL-1));
else
    tmp.MTFLoCorner = nan;
end;    

% High side corner.
indH = max(find(abovethr3dB));
if ~isempty(indH) && indH<length(fs)
    % Linearly interpolate to find a cut-off frequency. 
    tmp.MTFHiCorner = fs(indH) + ...
        (fs(indH+1)-fs(indH))*(thr3dB-valnorm(indH))/(valnorm(indH+1)-valnorm(indH));
else
    tmp.MTFHiCorner = nan;
end;    


% Decide what kind of function it is.
% And calculate a statistic of "goodness" of the passband. 
if isnan(tmp.MTFHiCorner) & isnan(tmp.MTFLoCorner)
    tmp.MTFCornerGoodness =  sum(valnorm>thr3dB)./length(valnorm);
    if tmp.MTFCornerGoodness>=.75
        tmp.MTFType = 'allpass';
    elseif tmp.MTFCornerGoodness<=.25
        tmp.MTFType = 'allstop';
    else
        tmp.MTFType = 'ambiguous';
    end;
    tmp.MTFCornerBW = nan;
elseif ~isnan(tmp.MTFHiCorner) & ~isnan(tmp.MTFLoCorner)    
    tmp.MTFType = 'bandpass';
    tmp.MTFCornerGoodness  = sum(valnorm(indL:indH)>thr3dB)/(indH-indL+1);
    tmp.MTFCornerBW = tmp.MTFHiCorner - tmp.MTFLoCorner;
elseif ~isnan(tmp.MTFHiCorner) & isnan(tmp.MTFLoCorner)
    tmp.MTFType = 'lowpass';
    tmp.MTFCornerGoodness  = sum(valnorm(1:indH)>thr3dB)/(indH);
    tmp.MTFCornerBW = tmp.MTFHiCorner - fs(1);
elseif isnan(tmp.MTFHiCorner) & ~isnan(tmp.MTFLoCorner)
    tmp.MTFType = 'hipass';
    tmp.MTFCornerGoodness  = sum(valnorm(indL:end)>thr3dB)/(length(valnorm)-indL);
    tmp.MTFCornerBW = fs(end) - tmp.MTFLoCorner;
end;


% Find the points that exceed the 31.6% point (-10dB point).
% "Cutoff frequencies"
thr10dB = 10^(-10/20);

% Low cut-off
abovethr10dB = valnorm>thr10dB;
indL = min(find(abovethr10dB));
if ~isempty(indL) & indL>1
    % Linearly interpolate to find a cut-off frequency. 
    tmp.MTFLoCutOff = fs(indL-1) + ...
        (fs(indL)-fs(indL-1))*(thr10dB-valnorm(indL-1))/(valnorm(indL)-valnorm(indL-1));
else
    tmp.MTFLoCutOff = nan;
end;    

% High side corner.
indH = max(find(abovethr10dB));
if ~isempty(indH) && indH<length(fs)
    % Linearly interpolate to find a cut-off frequency. 
    tmp.MTFHiCutOff = fs(indH) + ...
        (fs(indH+1)-fs(indH))*(thr10dB-valnorm(indH))/(valnorm(indH+1)-valnorm(indH));
else
    tmp.MTFHiCutOff = nan;
end;    
   
%  calculate a statistic of "goodness" of the passband for the cut-off. 
if isnan(tmp.MTFHiCutOff) & isnan(tmp.MTFLoCutOff)
    tmp.MTFCutoffGoodness =  sum(valnorm>thr10dB)./length(valnorm);
    tmp.MTFCutOffBW = nan;
elseif ~isnan(tmp.MTFHiCutOff) & ~isnan(tmp.MTFLoCutOff)    
    tmp.MTFCutoffGoodness  = sum(valnorm(indL:indH)>thr10dB)/(indH-indL+1);
    tmp.MTFCutOffBW = tmp.MTFHiCutOff - tmp.MTFLoCutOff;
elseif ~isnan(tmp.MTFHiCutOff) & isnan(tmp.MTFLoCutOff)
    tmp.MTFCutoffGoodness  = sum(valnorm(1:indH)>thr10dB)/(indH);
    tmp.MTFCutOffBW = tmp.MTFHiCutOff - fs(1);
elseif isnan(tmp.MTFHiCutOff) & ~isnan(tmp.MTFLoCutOff)
    tmp.MTFCutoffGoodness  = sum(valnorm(indL:end)>thr10dB)/(length(valnorm)-indL);
    tmp.MTFCutOffBW = fs(end) - tmp.MTFLoCutOff;
end;


% Absolute threshold - applied without normalisation.
if nargin>3
    % Low cut-off
    abovethr = vals>absthr;
    indL = min(find(abovethr));
    if ~isempty(indL) & indL>1
        % Linearly interpolate to find a cut-off frequency. 
        tmp.MTFLoFAbsThr = fs(indL-1) + ...
            (fs(indL)-fs(indL-1))*(absthr-vals(indL-1))/(vals(indL)-vals(indL-1));
    else
        tmp.MTFLoFAbsThr = nan;
    end;    

    % High side corner.
    indH = max(find(abovethr));
    if ~isempty(indH) && indH<length(fs)
        % Linearly interpolate to find a cut-off frequency. 
        tmp.MTFHiFAbsThr = fs(indH) + ...
            (fs(indH+1)-fs(indH))*(absthr-vals(indH))/(vals(indH+1)-vals(indH));
    else
        tmp.MTFHiFAbsThr = nan;
    end;    

    %  calculate a statistic of "goodness" of the passband for the cut-off. 
    if isnan(tmp.MTFHiFAbsThr) & isnan(tmp.MTFLoFAbsThr)
        tmp.MTFAbsThrGoodness =  sum(valnorm>absthr)./length(valnorm);
        tmp.MTFAbsThrBW = nan;
    elseif ~isnan(tmp.MTFHiFAbsThr) & ~isnan(tmp.MTFLoFAbsThr)    
        tmp.MTFAbsThrGoodness  = sum(valnorm(indL:indH)>absthr)/(indH-indL+1);
        tmp.MTFAbsThrBW = tmp.MTFHiFAbsThr - tmp.MTFLoFAbsThr;
    elseif ~isnan(tmp.MTFHiFAbsThr) & isnan(tmp.MTFLoFAbsThr)
        tmp.MTFAbsThrGoodness  = sum(valnorm(1:indH)>absthr)/(indH);
        tmp.MTFAbsThrBW = tmp.MTFHiFAbsThr - fs(1);
    elseif isnan(tmp.MTFHiFAbsThr) & ~isnan(tmp.MTFLoFAbsThr)
        tmp.MTFAbsThrGoodness  = sum(valnorm(indL:end)>absthr)/(length(valnorm)-indL);
        tmp.MTFAbsThrBW = fs(end) - tmp.MTFLoFAbsThr;
    end;
end;


