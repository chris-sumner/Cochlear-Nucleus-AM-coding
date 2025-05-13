function dataout = stack_Datasets(datain,varargin)
% stack_Datasets - makes cell vector of all conditions.
%  dataout = stack_Datasets(datain,option)
% N.B. Have to decide which of the AM lists you are stacking with.
%   option: 'allamfreqs', 'usedamfreqs','usedamfreqs_dprime'

% Fields which do not vary with AM frequency.
fieldList_common = {'unitid','modLevel','amToneDur','carrFreq','depthMod', ...
                    'cf','cf_thr','type', 'typeNum', 'CV','NumFreqMoreThan3Peaks', ...
                    'NumFreqthatPhaseLock','fullDataUnitIndex','rationalisedType', ...
                    'BMF','BMFdp' , ...
                    'VS_LoCornerHz', 'VS_HiCornerHz', 'VS_LoCutOffHz', ...
                    'VS_HiCutOffHz', 'VS_MTFType', 'VS_MTFCornerGoodness', ...
                    'VS_MTFCutOffGoodness', 'VS_MTFCornerBW', 'VS_MTFCutOffBW', ...
                    'dprime_LoCornerHz', 'dprime_HiCornerHz', 'dprime_LoCutOffHz', ...
                    'dprime_HiCutOffHz', 'dprime_MTFType', 'dprime_MTFCornerGoodness', ...
                    'dprime_MTFCutOffGoodness', 'dprime_MTFCornerBW', 'dprime_MTFCutOffBW', ...
                     'VS_HiFAbsThr','VS_LoFAbsThr','VS_MTFAbsThrGoodness','VS_MTFAbsThrBW', ...
                    'dprime_HiFAbsThr','dprime_LoFAbsThr','dprime_MTAbsThrGoodness','dprime_MTFAbsThrBW', ...
                    'VS_max','VS_meanBelow1k','VS_150ish', ...
                    'dprimes_max', ...
                    'CI_max','CI_meanBelow1k','CI_150ish', ...
                    'Z_max','Z_mean','Z_meanBelow1k','Z_150ish', ...
                    'numSACpeaks_p01_max','numSACpeaks_p001_max','numSACpeaks_p0001_max', ...
                    'numSACpeaks_p01_meanBelow1k','numSACpeaks_p001_meanBelow1k','numSACpeaks_p0001_meanBelow1k', ...
                    'numSACpeaks_p01_150ish','numSACpeaks_p001_150ish','numSACpeaks_p0001_150ish', ...
                    'Z_NHPP_max','Z_NHPP_mean','Z_NHPP_meanBelow1k','Z_NHPP_150ish', ...
                    'spikesperperiod_max','spikesperperiod_mean','spikesperperiod_meanBelow1k','spikesperperiod_150ish' ... 
                    'Z_fMod1and2','envFluct_0p48ms_fMod1and2','reliability_R_0p24ms_fMod1and2','numSACpeaks_p0001_fMod1and2', ...
                    'SACpeak0sal_150ish_p0001','SACpeak1sal_150ish_p0001','SACpeak2sal_150ish_p0001', ...
                    'SACpeak3sal_150ish_p0001','SACpeak4sal_150ish_p0001','SACpeak5sal_150ish_p0001'
                   };
% Fields which do not vary contain all AM frequency.
fieldList_all = {'allamfreqs','allamfreqs_reBMF','allamfreqs_reBMFdp'};
 % Fields which have only 'used' AM frequency.   
fieldList_used = {'usedamfreqs','VS','Rayleigh','usedamfreqs_reBMF','usedamfreqs_reBMFdp','CI','CI_p','Z', ...
                  'spikesperperiod','ref0s','MLPoissTest','numberofpeaks','sigpeaklags','VS_reBestVS','VS_pcBestVS'};
 % Fields which have the dprime set of AM frequencies.   
fieldList_dprime = {'usedamfreqs_dprime','usedamfreqs_dprime_reBMF','usedamfreqs_dprime_reBMFdp','dprimes', ...
    'dprimes_reBestDP','dprimes_pcsBestDP','dprimeReDp1'};
% Fields which are not stackable.
fieldList_dontStack = {'type','rationalisedType','amType','dprime_MTFType','VS_MTFType'};
% 'sigpeaklags',

COMMON = 0;
ALL = 1;
USED = 2;
DPRIME = 3;

% Find the AM frequency option
if strcmp('allamfreqs',varargin)
    amstackoption = ALL; 
    amstackfield = 'allamfreqs'; 
elseif strcmp('usedamfreqs',varargin)
    amstackoption = USED; 
    amstackfield = 'usedamfreqs'; 
elseif strcmp('dprime',varargin)
    amstackoption = DPRIME;
    amstackfield = 'usedamfreqs_dprime'; 
end;

% Get all the fields.
fieldlist = fieldnames(datain);

% Loop through the fields.
for fi = 1:length(fieldlist)
    % Empty temporary output. 
    tmpout = cell(  size(datain.(fieldlist{fi}) ) );
    
    % Which AM frequency list is the field from? (length).
    if sum(strcmp(fieldlist{fi},fieldList_all)>0) 
        amtype = ALL; 
        amfield = 'allamfreqs';
    elseif sum(strcmp(fieldlist{fi},fieldList_used)>0)
        amtype = USED; 
        amfield = 'usedamfreqs';
    elseif sum(strcmp(fieldlist{fi},fieldList_dprime)>0)
        amtype = DPRIME;
        amfield = 'usedamfreqs_dprime';
    elseif sum(strcmp(fieldlist{fi},fieldList_common)>0)
        amtype = COMMON;
        amfield = [];
    end;
  
    if sum(strcmp(fieldlist{fi},fieldList_dontStack)>0)
        stackflag = false;
    else 
        stackflag = true;
    end;
    
    
    % Loop through the datasets for a single field.
    for di = 1:length(datain.(fieldlist{fi}))
 
        % The list am frequencies for stacking against.
        amstacklist = datain.(amstackfield){di};
        
        % if this list is not empty...
        if  ~isempty(amstacklist)
            % Initialise to nans - might not have values. 
            if stackflag
                tmpout{di} = nan(size(amstacklist))';
            else
                tmpout{di} = cell(size(amstacklist))';
            end;

            if ~isempty(amfield)  
                % Where there is actually condition for each AM frequency.
                amlist = datain.(amfield){di}; 
                % The list am frequencies for stacking against.
                amstacklist = datain.(amstackfield){di};

                % Work out where the values go. 
                [lia locb] = ismember( amlist, amstacklist );
                nlia =  sum(lia)>0;
                ldata = length(datain.(fieldlist{fi}){di});

                % stick em in there. 
                if stackflag && nlia>0 && ldata>=nlia
                    % most fields
                    onetmp = datain.(fieldlist{fi}){di}(find(lia));
                    if ~iscell(onetmp)                    
                        tmpout{di} = nan(size(amstacklist))';  % Initialise to nans
                        tmpout{di}(locb(lia)) = onetmp;
                    else
                        tmpout{di} = cell(size(amstacklist))';  % Initialse to empty cells
                        
                        % Complicated if each entry is a cell array itself.
                        % This is the case for 'sigpeaklags'                        
                        emptyentries = find(cellfun(@(x) isempty(x),  onetmp));
                        onetmp( emptyentries ) = arrayfun(@(x) nan, emptyentries,'uni',false);
                        
                        tmpout{di}(locb(lia)) = onetmp;
                        
                    end;
                elseif ~stackflag  && nlia>0 && ldata>=nlia                    
                    tmpout{di} = cell(size(amstacklist))';  % Initialse to empty cells

                    % one field cannot be reduced this way - so just
                    % inherit it.
                    tmpout{di}(locb(lia)) = {datain.(fieldlist{fi}){di}(find(lia))};
                end;  

            else
                if stackflag
                    tmpout{di} = datain.(fieldlist{fi}){di} * ones( size( amstacklist ) )';
                else 
                    for ami = 1:length(amstacklist)
                        tmpout{di}{ami} = datain.(fieldlist{fi}){di}; % * ones( size( amstacklist ) )';
                    end;
                end;
            end;
        end;
    end;
      
    % Pump them all out into one long vector. 
    dataout.(fieldlist{fi}) = [ tmpout{:} ];
end;


