function statsmat = collate_byDataset(unitCIandVS,unitWR,varargin)
% collate_byDataset brings together two separate sets of analysis.

% This option allows us to specific the classifier time-constant manually.
if nargin>2&&  strcmp(varargin{1},'classifier_tau') 
    % Select specified tau for each set. 
    unitWR = getWROutputsWithTau(unitWR,varargin{2});    
else
    % Select best tau for each set.
    unitWR = getBestWROutputs(unitWR);     
end;

% This option allows us to subsitute another performance measure for the
% default. Allows easy swapping of measures for comparison.
if nargin>2&&  strcmp(varargin{1},'performance') 
    % Select the measure of performance.
    perf_measure = varargin{2};
else
    % Go with the default value.
    perf_measure = 'dprime';
end;
% I did look at introducing different names, but there is a lot to unpick
% here as we started with the d' measure. 

numDatasets = min([length(unitCIandVS), length(unitWR)] );

fprintf('Collating data...');
fprintf('\nPercent complete:');

% ----------------------  Load up  data  ------------------------

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
    

    % ----- retrieve the AM type and double check the data  ----------------
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

    % ------------ Summary stats for Rayleigh --------------
    rmt = nanmax( statsmat.Rayleigh{i} );
    if ~isempty(rmt)
        statsmat.Rayleigh_max{i} = rmt;
        statsmat.Rayleigh_mean{i} = nanmean( statsmat.Rayleigh{i} );
    else
        statsmat.Rayleigh_max{i} = nan;
        statsmat.Rayleigh_mean{i} = nan;
    end;       
    
    % Mean SPP below 1k. 
    below1kInds = find( statsmat.usedamfreqs{i}<1000 );   
    if ~isempty(below1kInds)
        statsmat.Rayleigh_meanBelow1k{i} = nanmean(statsmat.Rayleigh{i}(below1kInds));
    else
        statsmat.Rayleigh_meanBelow1k{i} = nan;
    end;
    % Low frequency SPP only
    closeto150Inds = find( (statsmat.usedamfreqs{i}==150 | statsmat.usedamfreqs{i}==200));  
    if ~isempty(closeto150Inds)
        statsmat.Rayleigh_150ish{i} = statsmat.Rayleigh{i}(closeto150Inds(1));
    else
        statsmat.Rayleigh_150ish{i} = nan;
    end;



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
                        
  % -------------------- Corelation Index ------------------------
    
    statsmat.CI{i} =  arrayfun(@(x) x.CIvalue_0lag,unitCIandVS(i).CI)';
    statsmat.CI_p{i} =  arrayfun(@(x) x.p,unitCIandVS(i).CI)';
    statsmat.CI_max{i} =  arrayfun(@(x) x.CIvalue_max,unitCIandVS(i).CI)';   

    if isempty(statsmat.CI_max{i})
        statsmat.CI_max{i} = nan;
    end;
    % Mean CI below 1k. 
    below1kInds = find(statsmat.usedamfreqs{i}<1000 & statsmat.CI_p{i}>.99);   
    if ~isempty(below1kInds)
        statsmat.CI_meanBelow1k{i} = nanmean(statsmat.CI{i}(below1kInds));
    else
        statsmat.CI_meanBelow1k{i} = nan;
    end;
    % Low frequency CI only
    closeto150Inds = find( (statsmat.usedamfreqs{i}==150 | statsmat.usedamfreqs{i}==200) & statsmat.CI_p{i}>.99);  
    if ~isempty(closeto150Inds)
        statsmat.CI_150ish{i} = statsmat.CI{i}(closeto150Inds(1));
    else
        statsmat.CI_150ish{i} = nan;
    end;

    % -------------------- SAC peaks -----------------------------

    % Main statistic required here is the number of peaks in the SAC.
    % We limit to those lags <1.1 x the period. The expected value is 3.
    statsmat.numSACpeaks_p0001{i} = ...
            arrayfun(@(x) sum(abs(x.setofsigpeaklags_p0001)<x.am_period*1.1), ...
                        unitCIandVS(i).SACpeaks)';

    % Set empty data (no peak analysis done) to nan
    nani= arrayfun(@(x) isempty(x.setofsigpeaklags_p0001), ...
                        unitCIandVS(i).SACpeaks)';
    statsmat.numSACpeaks_p0001{i}(nani) = nan;

    % Find the man, mean and values close to 150Hz.
    statsmat.numSACpeaks_p0001_max{i} = FnOrNan(statsmat.numSACpeaks_p0001{i},@max);
    
    % Mean CI below 1k. 
    statsmat.numSACpeaks_p0001_meanBelow1k{i} = ...
        FnOrNan(statsmat.numSACpeaks_p0001{i},@mean,statsmat.usedamfreqs{i}<1000);

     % Low frequency only
     closeto150Inds = find( (statsmat.usedamfreqs{i}==150 | statsmat.usedamfreqs{i}==200) & statsmat.CI_p{i}>.99);
     statsmat.numSACpeaks_p0001_150ish{i} = FnOrNan(statsmat.numSACpeaks_p0001{i},@max,closeto150Inds);

     % Extracting the saliences of the peaks.
     % We only look at the postive peaks   
     tmpsals = arrayfun(@(x) x.setofsigsaliences_p0001( ...
                      x.setofsigpeaklags_p0001>x.am_period*0.1 & x.setofsigpeaklags_p0001<x.am_period*0.9), ...
                   unitCIandVS(i).SACpeaks, 'uni',false)';    

     tmplags = arrayfun(@(x) x.setofsigpeaklags_p0001( ...
                      x.setofsigpeaklags_p0001>x.am_period*0.1 & x.setofsigpeaklags_p0001<x.am_period*0.9), ...
                   unitCIandVS(i).SACpeaks, 'uni',false)';    
           
     % Look separately for the peak closest to zero lag.           
     zerosal = arrayfun(@(x) x.setofsigsaliences_p0001( ...
                      x.setofsigpeaklags_p0001>x.am_period*-0.1 & x.setofsigpeaklags_p0001<=x.am_period*0.1), ...
                   unitCIandVS(i).SACpeaks, 'uni',false)';    

     zerolag = arrayfun(@(x) x.setofsigpeaklags_p0001( ...
                      x.setofsigpeaklags_p0001>x.am_period*-0.1 & x.setofsigpeaklags_p0001<=x.am_period*0.1), ...
                   unitCIandVS(i).SACpeaks, 'uni',false)';    
             
     % Look separately for the peak closest to the period lag.           
     periodsal = arrayfun(@(x) x.setofsigsaliences_p0001( ...
                      x.setofsigpeaklags_p0001>x.am_period*0.9 & x.setofsigpeaklags_p0001<=x.am_period*1.1), ...
                   unitCIandVS(i).SACpeaks, 'uni',false)';    

     periodlag = arrayfun(@(x) x.setofsigpeaklags_p0001( ...
                      x.setofsigpeaklags_p0001>x.am_period*0.9 & x.setofsigpeaklags_p0001<=x.am_period*1.1), ...
                   unitCIandVS(i).SACpeaks, 'uni',false)';    

               
     % --- 0th peak the first peak (which less than the period delay) ---
     % Intialise to zero for all modulation frequencies
     statsmat.SACpeak0sal_p0001{i} = zeros(1,length(zerosal));    
     statsmat.SACpeak0lag_p0001{i} = nan(1,length(zerosal));    
     % Find which modulation frequencies have a peak.   
     tmpzerosalind = find(cellfun(@(x) length(x)>=1,  zerosal));
     statsmat.SACpeak0sal_p0001{i}(tmpzerosalind) = cellfun(@(x) x(1),zerosal(tmpzerosalind));
     statsmat.SACpeak0lag_p0001{i}(tmpzerosalind) = cellfun(@(x) x(1),zerolag(tmpzerosalind));
               
     % --- For the first peak (which should be close to the period delay) ---
     % Intialise to zero for all modulation frequencies
     statsmat.SACpeak1sal_p0001{i} = zeros(1,length(tmpsals));    
     statsmat.SACpeak1lag_p0001{i} = nan(1,length(tmpsals));    
     % Find which modulation frequencies have 1 peak or more.   
     tmpsalinds = find(cellfun(@(x) length(x)>=1,  periodsal));
     % Pop those in the array   
     statsmat.SACpeak1sal_p0001{i}(tmpsalinds) = cellfun(@(x) x(1),periodsal(tmpsalinds));
     statsmat.SACpeak1lag_p0001{i}(tmpsalinds) = cellfun(@(x) x(1),periodlag(tmpsalinds));
     % Any non-significant values are then zero.

     % --- Second peak the first peak (which less than the period delay) ---
     % Intialise to zero for all modulation frequencies
     statsmat.SACpeak2sal_p0001{i} = zeros(1,length(tmpsals));    
     statsmat.SACpeak2lag_p0001{i} = nan(1,length(tmpsals));    
    % Find which modulation frequencies have 1 peak or more besides the period peak.   
     tmpsalinds = find(cellfun(@(x) length(x)>=1,  tmpsals));
     statsmat.SACpeak2sal_p0001{i}(tmpsalinds) = cellfun(@(x) x(1),tmpsals(tmpsalinds));
     statsmat.SACpeak2lag_p0001{i}(tmpsalinds) = cellfun(@(x) x(1),tmplags(tmpsalinds));
     
     % --- Third peak the first peak (which less than the period delay) ---
     % Intialise to zero for all modulation frequencies
     statsmat.SACpeak3sal_p0001{i} = zeros(1,length(tmpsals));    
     statsmat.SACpeak3lag_p0001{i} = nan(1,length(tmpsals));    
    % Find which modulation frequencies have 2 peak or more besides the period peak. 
     tmpsalinds = find(cellfun(@(x) length(x)>=2,  tmpsals));
     statsmat.SACpeak3sal_p0001{i}(tmpsalinds) = cellfun(@(x) x(2),tmpsals(tmpsalinds));
     statsmat.SACpeak3lag_p0001{i}(tmpsalinds) = cellfun(@(x) x(2),tmplags(tmpsalinds));

     % --- 4th peak the first peak (which less than the period delay) ---
     % Intialise to zero for all modulation frequencies
     statsmat.SACpeak4sal_p0001{i} = zeros(1,length(tmpsals));    
     statsmat.SACpeak4lag_p0001{i} = nan(1,length(tmpsals));    
    % Find which modulation frequencies have 3 peaks or more besides the period peaks.   
     tmpsalinds = find(cellfun(@(x) length(x)>=3,  tmpsals));
     statsmat.SACpeak4sal_p0001{i}(tmpsalinds) = cellfun(@(x) x(3),tmpsals(tmpsalinds));
     statsmat.SACpeak4lag_p0001{i}(tmpsalinds) = cellfun(@(x) x(3),tmplags(tmpsalinds));
     
     % --- 5th peak the first peak (which less than the period delay) ---
     % Intialise to zero for all modulation frequencies
     statsmat.SACpeak5sal_p0001{i} = zeros(1,length(tmpsals));    
     statsmat.SACpeak5lag_p0001{i} = nan(1,length(tmpsals));    
     % Find which modulation frequencies have 4 peaks or more besides the period peak  
     tmpsalinds = find(cellfun(@(x) length(x)>=4,  tmpsals));
     statsmat.SACpeak5sal_p0001{i}(tmpsalinds) = cellfun(@(x) x(4),tmpsals(tmpsalinds));
     statsmat.SACpeak5lag_p0001{i}(tmpsalinds) = cellfun(@(x) x(4),tmplags(tmpsalinds));

     % Simple sum of all the peaks. 
     statsmat.SACpeakTotalSal_p0001{i} = statsmat.SACpeak0sal_p0001{i} + statsmat.SACpeak1sal_p0001{i} + ...
                                         statsmat.SACpeak2sal_p0001{i} + statsmat.SACpeak3sal_p0001{i} + ...
                                         statsmat.SACpeak4sal_p0001{i} + statsmat.SACpeak5sal_p0001{i};     
     
     
     % Low frequency saliences only
     closeto150Inds = find( (statsmat.usedamfreqs{i}==150 | statsmat.usedamfreqs{i}==200));  
     if ~isempty(closeto150Inds)
        statsmat.SACpeak0lag_150ish_p0001{i} = statsmat.SACpeak0lag_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak1lag_150ish_p0001{i} = statsmat.SACpeak1lag_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak2lag_150ish_p0001{i} = statsmat.SACpeak2lag_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak3lag_150ish_p0001{i} = statsmat.SACpeak3lag_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak4lag_150ish_p0001{i} = statsmat.SACpeak4lag_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak5lag_150ish_p0001{i} = statsmat.SACpeak5lag_p0001{i}(closeto150Inds(1));        
 
        statsmat.SACpeak0sal_150ish_p0001{i} = statsmat.SACpeak0sal_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak1sal_150ish_p0001{i} = statsmat.SACpeak1sal_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak2sal_150ish_p0001{i} = statsmat.SACpeak2sal_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak3sal_150ish_p0001{i} = statsmat.SACpeak3sal_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak4sal_150ish_p0001{i} = statsmat.SACpeak4sal_p0001{i}(closeto150Inds(1));
        statsmat.SACpeak5sal_150ish_p0001{i} = statsmat.SACpeak5sal_p0001{i}(closeto150Inds(1));        
        statsmat.SACpeakTotalSal_150ish_p0001{i} = statsmat.SACpeakTotalSal_p0001{i}(closeto150Inds(1));        
     else
        statsmat.SACpeak0lag_150ish_p0001{i} = nan;
        statsmat.SACpeak1lag_150ish_p0001{i} = nan;
        statsmat.SACpeak2lag_150ish_p0001{i} = nan;
        statsmat.SACpeak3lag_150ish_p0001{i} = nan;
        statsmat.SACpeak4lag_150ish_p0001{i} = nan;
        statsmat.SACpeak5lag_150ish_p0001{i} = nan;        

        % N.B. Nan's not zeros, because the data does not exist.  
        % Rather than there being no peaks.
        statsmat.SACpeak0sal_150ish_p0001{i} = nan;
        statsmat.SACpeak1sal_150ish_p0001{i} = nan;
        statsmat.SACpeak2sal_150ish_p0001{i} = nan;
        statsmat.SACpeak3sal_150ish_p0001{i} = nan;
        statsmat.SACpeak4sal_150ish_p0001{i} = nan;
        statsmat.SACpeak5sal_150ish_p0001{i} = nan;
        statsmat.SACpeakTotalSal_150ish_p0001{i} = nan;
     end;     

     
    % ---------------------- Z-scores -----------------------------
    
    % Raw Z-score
    statsmat.Z{i} = arrayfun(@(x) x.stats.vals(5),unitCIandVS(i).ML);
    statsmat.Z_NHPP{i} = arrayfun(@(x) x.Z_NHPP,unitCIandVS(i).ML);   % NHPP based Z-score. 

    % The two "tests" of mode-locking from Laudanski et al. 2010.
    statsmat.ModeLockingTest_Z{i} = arrayfun(@(x) x.ModeLockingTest_Z,unitCIandVS(i).ML);
    statsmat.ModeLockingTest_ISIshuf{i} = arrayfun(@(x) x.ModeLockingTest_ISIshuf,unitCIandVS(i).ML);
    % Surrogate phase locking statistics - after interval shuffling. 
    statsmat.Rayleigh_Isur{i} = arrayfun(@(x) x.stats.vals(4),unitCIandVS(i).ML);
    statsmat.VS_Isur{i} = arrayfun(@(x) x.stats.vals(2),unitCIandVS(i).ML);

    % Single statistics from each record for population statistical models.
    zmt = nanmax( statsmat.Z{i} );
    % Mean and max values.
    if ~isempty(zmt)
        statsmat.Z_max{i} = zmt;
        statsmat.Z_mean{i} = nanmean( statsmat.Z{i} );
    else
        statsmat.Z_max{i} = nan;
        statsmat.Z_mean{i} = nan;
    end;       
    
    % Single statistics from each record for population statistical models.
    zmt = nanmax( statsmat.Z_NHPP{i} );
    % Mean and max values.
    if ~isempty(zmt)
        statsmat.Z_NHPP_max{i} = zmt;
        statsmat.Z_NHPP_mean{i} = nanmean( statsmat.Z_NHPP{i} );
    else
        statsmat.Z_NHPP_max{i} = nan;
        statsmat.Z_NHPP_mean{i} = nan;
    end;       

    % ISI shuffled Rayleigh statistics from each record for population statistical models.
    rmt = nanmax( statsmat.Rayleigh_Isur{i} );
    % Mean and max values.
    if ~isempty(rmt)
        statsmat.Rayleigh_Isur_max{i} = rmt;
        statsmat.Rayleigh_Isur_mean{i} = nanmean( statsmat.Rayleigh_Isur{i} );
    else
        statsmat.Rayleigh_Isur_max{i} = nan;
        statsmat.Rayleigh_Isur_mean{i} = nan;
    end;       
    
    % ISI shuffled VS from each record.
    vmt = nanmax( statsmat.VS_Isur{i} );
    % Mean and max values.
    if ~isempty(rmt)
        statsmat.VS_Isur_max{i} = vmt;
        statsmat.VS_Isur_mean{i} = nanmean( statsmat.VS_Isur{i} );
    else
        statsmat.VS_Isur_max{i} = nan;
        statsmat.VS_Isur_mean{i} = nan;
    end;       

    % Means below 1k. 
    below1kInds = find( statsmat.usedamfreqs{i}<1000 );   
    if ~isempty(below1kInds)
        statsmat.Z_meanBelow1k{i} = nanmean(statsmat.Z{i}(below1kInds));
        statsmat.Z_NHPP_meanBelow1k{i} = nanmean(statsmat.Z_NHPP{i}(below1kInds));
        statsmat.Rayleigh_Isur_meanBelow1k{i} = nanmean(statsmat.Rayleigh_Isur{i}(below1kInds));
        statsmat.VS_Isur_meanBelow1k{i} = nanmean(statsmat.VS_Isur{i}(below1kInds));
    else
        statsmat.Z_NHPP_meanBelow1k{i} = nan;
        statsmat.Z_meanBelow1k{i} = nan;
        statsmat.Rayleigh_Isur_meanBelow1k{i} = nan;
        statsmat.VS_Isur_meanBelow1k{i} = nan;
    end;
 
    % Low frequency Z only
    closeto150Inds = find( (statsmat.usedamfreqs{i}==150 | statsmat.usedamfreqs{i}==200));  
    if ~isempty(closeto150Inds)
        statsmat.Z_150ish{i} = statsmat.Z{i}(closeto150Inds(1));
        statsmat.Z_NHPP_150ish{i} = statsmat.Z_NHPP{i}(closeto150Inds(1));
        statsmat.Rayleigh_Isur_150ish{i} = statsmat.Rayleigh_Isur{i}(closeto150Inds(1));
        statsmat.VS_Isur_150ish{i} = statsmat.VS_Isur{i}(closeto150Inds(1));
    else
        statsmat.Z_150ish{i} = nan;
        statsmat.Z_NHPP_150ish{i} = nan;
        statsmat.Rayleigh_Isur_150ish{i} = nan;
        statsmat.VS_Isur_150ish{i} = nan;
    end;



    % --------------------- spikes per period --------------------
    statsmat.spikesperperiod{i} = arrayfun(@(x) x.spikecounts.spikesperperiod, ...
                                                  unitCIandVS(i).ML);
    smt = nanmax( statsmat.spikesperperiod{i} );
    if ~isempty(smt)
        statsmat.spikesperperiod_max{i} = smt;
        statsmat.spikesperperiod_mean{i} = nanmean( statsmat.spikesperperiod{i} );
    else
        statsmat.spikesperperiod_max{i} = nan;
        statsmat.spikesperperiod_mean{i} = nan;
    end;       
    
    % Mean SPP below 1k. 
    below1kInds = find( statsmat.usedamfreqs{i}<1000 );   
    if ~isempty(below1kInds)
        statsmat.spikesperperiod_meanBelow1k{i} = nanmean(statsmat.spikesperperiod{i}(below1kInds));
    else
        statsmat.spikesperperiod_meanBelow1k{i} = nan;
    end;
    % Low frequency SPP only
    closeto150Inds = find( (statsmat.usedamfreqs{i}==150 | statsmat.usedamfreqs{i}==200));  
    if ~isempty(closeto150Inds)
        statsmat.spikesperperiod_150ish{i} = statsmat.spikesperperiod{i}(closeto150Inds(1));
    else
        statsmat.spikesperperiod_150ish{i} = nan;
    end;
                                              
    % --------------------- Gai and Carney statistics --------------------

    if ~isempty(unitCIandVS(i).GC)

        twindows = unique([unitCIandVS(i).GC.temporalWindows]);
    
        % Low frequency index.
        closeto150Inds = find( (statsmat.usedamfreqs{i}==150 | statsmat.usedamfreqs{i}==200));  
    
        % Loop through each analysis window.
        for wi = 1:length(twindows)
    
            % Make a string to call the new field for reliability. 
            tstr =   num2str(twindows(wi));
            tstr(strfind(tstr,'.')) = 'p';
            datstr = ['reliability_R_' tstr 'ms' ];
    
            % Store reliability. 
            statsmat.(datstr){i} = arrayfun(@(x) x.reliability_R(wi), ...
                                                     unitCIandVS(i).GC);
    
            % Also store the reliability at 150Hz.
            dat150str = ['reliability_R150ish_' tstr 'ms' ];
            if ~isempty(closeto150Inds)
                statsmat.(dat150str){i} =  statsmat.(datstr){i}(closeto150Inds(1));
            else
                statsmat.(dat150str){i} = nan;
            end;
    
            % Store envelope fluctuations.
            datstr = ['envFluct_' tstr 'ms' ];
            statsmat.(datstr){i} = arrayfun(@(x) x.envFluct(wi), ...
                                                     unitCIandVS(i).GC);
            % Also store the reliability at 150Hz.
            dat150str = ['envFluct_R150ish_' tstr 'ms' ];
            if ~isempty(closeto150Inds)
                statsmat.(dat150str){i} =  statsmat.(datstr){i}(closeto150Inds(1));
            else
                statsmat.(dat150str){i} = nan;
            end;
        end;

    else
         % Loop through each analysis window to store dummy values.
        for wi = 1:length(twindows)
    
            % Store dummies for each one.
            % N.B. Only works because temporal windows are the same for all
            % units and the first one is not empty!

            % Make a string to call the new field for reliability. 
            tstr =   num2str(twindows(wi));
            tstr(strfind(tstr,'.')) = 'p';
            datstr = ['reliability_R_' tstr 'ms' ];
    
            % Store reliability. 
            statsmat.(datstr){i} = [];
    
            % Also store the reliability at 150Hz.
            dat150str = ['reliability_R150ish_' tstr 'ms' ];
            statsmat.(dat150str){i} = nan;
    
            % Store envelope fluctuations.
            datstr = ['envFluct_' tstr 'ms' ];
            statsmat.(datstr){i} = [];
            % Also store the reliability at 150Hz.
            dat150str = ['envFluct_R150ish_' tstr 'ms' ];
            statsmat.(dat150str){i} = nan;
        end;

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
    statsmat.classifier_tau{i} = nan;
   
   if ~isempty(unitWR(i).bestClassifier)
        statsmat.dprimes{i} = unitWR(i).wroutput.(perf_measure);
        
        % Store the tau for the best classifier
        statsmat.classifier_tau{i} = unitWR(i).wroutput.classifier.tau;
        
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

function value = FnOrNan(inputvalues,fun,criteria)
% Efficient function to find a particular value but return nan if empty.

if nargin<3
    criteria = ones(size(inputvalues));
end;
inputvalues = inputvalues(criteria);

if ~isempty(inputvalues)
    value = fun(inputvalues);
else
    value = nan;
end;
