function [dataout] = select_Datasets(datain, varargin)
% select_Datasets   select out those conditions and units of desire. 
%
% options:
%   'fn' generalised test function for inclusion.
%   'amfn' function to select out specific modulation frequencies. 
%   'oneperunit' condition to decide which 1 condition to take for each
%                neuron.

% Adds root path to any files.   
global bngFn VbngFn;
dataout = struct;

% Check parameters for conditions to select out.
ind = find(strcmp('fn',varargin));
n_field_fns = 0; field_fns = {};
if ~isempty(ind)
    for i=1:length(ind)
        n_field_fns = n_field_fns +1;
        field_fns{n_field_fns} = varargin{ind(i)+1};
    end;
end;

% Check parameters for which modulation conditions to select out.
ind = find(strcmp('amfn',varargin));
n_amfield_fns = 0; amfield_fns = {};
if ~isempty(ind)
    for i=1:length(ind)
        n_amfield_fns = n_amfield_fns +1;
        amfield_fns{n_amfield_fns} = varargin{ind(i)+1};
    end;
end;

% This is the condition for selecting a single dataset for each unit.
ind = find(strcmp('oneperunit',varargin));
if ~isempty(ind)
    oneper_fn =  varargin{ind+1};
end;

% Load up the unit list (again) - bit messy. 
load('rawdata\unitList');
colfn =  @(s) find(strcmp(EXPLOGLIST,s));

if 1
    % ----------------- selection criteria -------------- 
    
    counterOfIncl = 0;   % Counter of how many of the datasets are included. 

    % Check the AM type (clumsy loading of Scholes' stuff).
    %datain.AMtype = getAMTypes(datain); pause;
    
    % Now go through each dataset and determine whether to include it. 
    for uniti=1:length(datain.unitid)
        fprintf('Data in #%d:',uniti);
        
        % Section of different conditions specified by functions.
        % Things except AM frequency determine which datasets to keep
        chk = check_FieldList_common(datain,uniti,field_fns);
                
        % If this dataset passes all the inclusion criteria... 
        if chk
            counterOfIncl = counterOfIncl + 1;
            fprintf('output #%d.',counterOfIncl);
             
            % Select out the AM frequencies we are going to output - for all relevant fields.     
            % N.B. Takes all frequencies unless there is an option to do otherwise
            dataout = selectAMFreqsIndDataSet(datain,uniti,counterOfIncl,amfield_fns,dataout);
            
           % Store the unit index in the complete list (just in case).
            dataout.fullDataUnitIndex{counterOfIncl} = datain.fullDataUnitIndex{uniti};
        else
           % Otherwise this unit is not included in the output at all.
           fprintf(' Not included.');
        end;
                
    fprintf('\n');
    end;
    
end;



if sum(strcmp('oneperunit',varargin))>0
    % ------------ select out one dataset for each unit ------------
    % This is intended to ape the principle of Schole's but redoing it so
    % it works properly and doesn't take 200% modulation data if you don't
    % want it to. 
    
    % This is applied AFTER any other criteria. 
    
    % Find a list of unique units.
    fullunitlist = [dataout.unitid{:}]; 
    uniqueunits = unique(fullunitlist); 
    
    % List of which datasets to keep.
    datasetlist = [];
    
    % Loop through each unique unit...
    for ui = 1:length(uniqueunits)
        
        candidate_ui = find(fullunitlist == uniqueunits(ui));
        
        % Determine which one we select.
    
        % Find the field name in the test function.  
        fname = fieldnames(dataout);
        fieldchk = find(cellfun(@(x) ~isempty(findstr(oneper_fn,x)), fname));
        
        if  length(fieldchk)>0
            
            % Always select out the longest field name which matches.
            [maxfield maxi] = max(cellfun(@(x) length(x),fname(fieldchk)));
            fieldchk = fieldchk(maxi);
            
            % Get just that field and only the entries for the candiate datasets.
            eval([fname{fieldchk} ' = dataout.' fname{fieldchk} '(candidate_ui);']);
            
            % Evaluate the function.
            % N.B. this needs to output a unique index - which to keep.
            di = eval(oneper_fn);
            
            % If more than one meets the criteria - just take the first
            % one.
            if length(di)>1
                di = di(1);
            end;
            
        end;

        % Add the dataset to keep to the list of datasets to keep
        datasetlist = [datasetlist candidate_ui(di)];
    end;
    
    % Selectout those datasets to be kept. 
    % For every field. 
    dataout = structfun(@(x) x(datasetlist),dataout,'uni',false);
    
end;

% Added during revision
dpnotempty =  cellfun(@(x) ~isempty(x), dataout.dprimes);
dataout.dprimeReDp1(dpnotempty) = cellfun(@(x) x/x(1), dataout.dprimes(dpnotempty),'uni',false);



% Function to select out the required frequencies. 
function dataout = selectAMFreqsIndDataSet(datain,indexin,indexout,field_fns,inherited_dataout)

% Fields which do not vary with AM frequency.
fieldList_common = defnFieldList_common;

% Fields which contains the AM frequencies used for the dprime analysis.
fieldList_dprime = defnAMFieldList_dprime;

% Fields which do vary with all AM frequencies .
fieldList_all = defnAMFieldList_all;
    
% Fields which do vary with AM frequencies used for the other analyses.
fieldList_used = defnAMFieldList_used;

% Ineherit everything so far.
dataout = inherited_dataout;

% amList can be specified as an input parameter OR as criteria in field_fns
if isnumeric(field_fns)
    amList = field_fns;
else
    amList = [];
end;

% Carry over all the frequencies which  do not vary with modulation
% frequency.
for i=1:length(fieldList_common)
    dataout.(fieldList_common{i})(indexout) =  datain.(fieldList_common{i})(indexin);    
end;


% Work out the selected list of AM frequencies to include here. 
if isempty(amList)
    amList_all = check_AMList_Freq(datain,indexin,field_fns, ...
                    datain.allamfreqs{indexin},@defnAMFieldList_all);
else
    amList_all = amList;    
end;
            
% Carry over all the frequencies which use the allamfreqs list. 
for i=1:length(fieldList_all)
    dataout.(fieldList_all{i}){indexout} = selectAMInField( ...
                                amList_all, datain.allamfreqs{indexin}, ...
                                datain.(fieldList_all{i}){indexin} );    
end;

% Work out the list of AM frequencies to include here. 
if isempty(amList)
    amList_used = check_AMList_Freq(datain,indexin,field_fns, ...
                datain.usedamfreqs{indexin},@defnAMFieldList_used);
else
    amList_used = amList;    
end;

% Carry over the selected frequencies which use the usedamfreqs list. 
for i=1:length(fieldList_used)
    dataout.(fieldList_used{i}){indexout} = selectAMInField( ...
                                amList_used, datain.usedamfreqs{indexin}, ...
                                datain.(fieldList_used{i}){indexin} );    
end;

% Work out the list of AM frequencies to include here. 
if isempty(amList)
    amList_dprime = check_AMList_Freq(datain,indexin,field_fns, ...
                    datain.usedamfreqs_dprime{indexin},@defnAMFieldList_dprime);
else
    amList_dprime =  amList;    
end;

% Carry over the selected frequencies which use the usedamfreqs list. 
for i=1:length(fieldList_dprime)
    dataout.(fieldList_dprime{i}){indexout} = selectAMInField( ...
                                amList_dprime, datain.usedamfreqs_dprime{indexin}, ...
                                datain.(fieldList_dprime{i}){indexin} );    
end;

                
                
function fieldout = selectAMInField(AMlist, AMref, fieldin)
% Selects out the correct entries given a list of AM frequencies. 

commonVals = intersect(AMlist,AMref);
refInds = find(ismember(AMref,commonVals));
if ~isempty(fieldin)
    fieldout = fieldin(refInds);
else
    fprintf('Unexpected empty field. Odd.\n');
    fieldout = [];
end;


function thisdata = selectOneDataset(datain,indexin)
% Simply extracts out the individual data for indexed dataset.
% Used for testing inclusion criteria.

fieldlist = fieldnames(datain);

for i=1:length(fieldlist)    
    thisdata.(fieldlist{i}) = datain.(fieldlist{i}){indexin};    
end;

function fieldList_common = defnFieldList_common()

% Fields which do not vary with AM frequency.
fieldList_common = {'unitid','modLevel','amToneDur','carrFreq','depthMod', ...
                    'cf','cf_thr','type', 'typeNum', 'CV', ...
                    'rationalisedType','amType', ...
                    'modLvl_rationalised','BMF','BMFdp', ...
                    'VS_LoCornerHz', 'VS_HiCornerHz', 'VS_LoCutOffHz', ...
                    'VS_HiCutOffHz', 'VS_MTFType', 'VS_MTFCornerGoodness', ...
                    'VS_MTFCutOffGoodness', 'VS_MTFCornerBW', 'VS_MTFCutOffBW', ...
                    'dprime_LoCornerHz', 'dprime_HiCornerHz', 'dprime_LoCutOffHz', ...
                    'dprime_HiCutOffHz', 'dprime_MTFType', 'dprime_MTFCornerGoodness', ...
                    'dprime_MTFCutOffGoodness', 'dprime_MTFCornerBW', 'dprime_MTFCutOffBW', ...
                    'VS_HiFAbsThr','VS_LoFAbsThr','VS_MTFAbsThrGoodness','VS_MTFAbsThrBW', ...
                    'dprime_HiFAbsThr','dprime_LoFAbsThr','dprime_MTFAbsThrGoodness','dprime_MTFAbsThrBW', ...
                    'VS_max','VS_meanBelow1k','VS_150ish' , ...
                    'Rayleigh_max','Rayleigh_meanBelow1k','Rayleigh_150ish' , ...
                    'dprimes_max', ...
                    'Z_max','Z_mean','Z_meanBelow1k','Z_150ish', ...
                    'Z_NHPP_max','Z_NHPP_mean','Z_NHPP_meanBelow1k','Z_NHPP_150ish', ...
                    'Rayleigh_Isur_max','Rayleigh_Isur_mean','Rayleigh_Isur_meanBelow1k','Rayleigh_Isur_150ish' , ...
                    'VS_Isur_max','VS_Isur_mean','VS_Isur_meanBelow1k','VS_Isur_150ish' , ...
                    'CI_150ish','spikesperperiod_150ish', ...
                    'numSACpeaks_p0001_max', ...
                    'numSACpeaks_p0001_meanBelow1k', ...
                    'numSACpeaks_p0001_150ish', ...
                    'reliability_R150ish_0p08ms','reliability_R150ish_0p16ms', 'reliability_R150ish_0p24ms', ...
                    'reliability_R150ish_0p48ms', 'reliability_R150ish_0p64ms', 'reliability_R150ish_1p1ms', ...
                    'reliability_R150ish_1p6ms', ...
                    'envFluct_R150ish_0p08ms','envFluct_R150ish_0p16ms','envFluct_R150ish_0p24ms', ...
                    'envFluct_R150ish_0p48ms','envFluct_R150ish_0p64ms','reliability_R150ish_1p1ms', ...
                    'envFluct_R150ish_1p6ms', ...
                    'SACpeak0sal_150ish_p0001','SACpeak1sal_150ish_p0001','SACpeak2sal_150ish_p0001', ...
                    'SACpeak3sal_150ish_p0001','SACpeak4sal_150ish_p0001','SACpeak5sal_150ish_p0001', ...
                    'SACpeakTotalSal_150ish_p0001', ...
                    'SACpeak0lag_150ish_p0001','SACpeak1lag_150ish_p0001','SACpeak2lag_150ish_p0001', ...
                    'SACpeak3lag_150ish_p0001','SACpeak4lag_150ish_p0001','SACpeak5lag_150ish_p0001', ...    
                    'classifier_tau'
                    };

 function fieldList_dprime = defnAMFieldList_dprime()
% Fields which do  vary with AM frequencies used for the dprime analysis.
               
fieldList_dprime = {'usedamfreqs_dprime','usedamfreqs_dprime_reBMF','usedamfreqs_dprime_reBMFdp','dprimes', ...
                    'usedamfreqs_dprime_rationalised','usedamfreqs_dprime_reBMF_rationalised', ...
                    'usedamfreqs_dprime_reBMFdp_rationalised','dprimes_reBestDP','dprimes_pcBestDP'};
                
 function fieldList_all = defnAMFieldList_all()
% Fields which do vary with all AM frequencies .
               
fieldList_all = {'allamfreqs','allamfreqs_reBMF','allamfreqs_reBMFdp', ...
    'allamfreqs_rationalised','allamfreqs_reBMF_rationalised','allamfreqs_reBMFdp_rationalised'};
                
function fieldList_used = defnAMFieldList_used()
% Fields which do vary with AM frequencies used for the other analyses.
               
fieldList_used = {'usedamfreqs','VS','Rayleigh','usedamfreqs_reBMF','usedamfreqs_reBMFdp', ...
                  'usedamfreqs_rationalised','usedamfreqs_reBMF_rationalised','usedamfreqs_reBMFdp_rationalised', ...
                  'VS_reBestVS','VS_pcBestVS', ...
                  'Z', 'ModeLockingTest_Z', 'ModeLockingTest_ISIshuf','VS_Isur','Rayleigh_Isur',  ...
                  'Z_NHPP', ...
                  'CI','CI_p', ...
                  'numSACpeaks_p0001', ...
                  'reliability_R_0p08ms','reliability_R_0p16ms','reliability_R_0p24ms', ...
                  'reliability_R_0p48ms','reliability_R_0p64ms', 'reliability_R_1p1ms','reliability_R_1p6ms', ...
                  'envFluct_0p08ms','envFluct_0p16ms','envFluct_0p24ms','envFluct_0p48ms', ...
                  'envFluct_0p64ms', 'envFluct_1p1ms','envFluct_1p6ms', ...
                  'spikesperperiod', ...
                  'SACpeak1sal_p0001','SACpeak2sal_p0001','SACpeak3sal_p0001', ...
                  'SACpeak4sal_p0001','SACpeak5sal_p0001','SACpeak0sal_p0001', ...
                  'SACpeakTotalSal_p0001', ...
                  'SACpeak1lag_p0001','SACpeak2lag_p0001','SACpeak3lag_p0001', ...
                  'SACpeak4lag_p0001','SACpeak5lag_p0001','SACpeak0lag_p0001' ...
                  };
% Removed: 'numSACpeaks_p01','numSACpeaks_p001',                
                
function chk = check_FieldList_common(datain,indexin,field_fns)

% Get this one data set out. 
thisdata = selectOneDataset(datain,indexin);

% Check the inclusion criteria for this dataset.
f = 1; fcount = 0;
if length(field_fns)>0
    % Loop through all field_fns
    fi = 1;
    while fi<=length(field_fns) && prod(f)>0
        % Test if this one is a "common" field.
        fname = fieldnames(datain);
        %fname = fieldList_common;
        fieldchk = find(cellfun(@(x) ~isempty(findstr(field_fns{fi},x)), fname));
        
        %cellfun(@(x) ~isempty(findstr(field_fns{fi},x)), fieldList_common);
        if  length(fieldchk)>0
            [maxfield maxi] = max(cellfun(@(x) length(x),fname(fieldchk)));
            fieldchk = fieldchk(maxi);
            %fname =  fieldList_common{find(fieldchk)};
            ind = findstr(field_fns{fi},fname{fieldchk});                
            if ind==1
                cmdstr = ['thisdata.' field_fns{fi}];
            else
                cmdstr = [field_fns{fi}(1:ind-1) 'thisdata.' field_fns{fi}(ind:end)];
            end;
            
            fcount = fcount+1;            
            f(fcount) = eval(cmdstr);
        end;
        fi = fi + 1;
    end;
end;

% Passes if all the inclusion criteria pass.
chk = prod(f);



                
function amList = check_AMList_Freq(datain,indexin,field_fns,fullAMList,defnFieldList_Freq)

% Define those fields which vary with AM frequency.
fieldList_Freq = defnFieldList_Freq();

% Get this one data set out. 
thisdata = selectOneDataset(datain,indexin);

% Check the inclusion criteria for this dataset.
f = ones(1,length(fullAMList)); 
fcount = 0;
if length(field_fns)>0
    % Loop through all field_fns
    for fi=1:length(field_fns)
        % Test if this one is a "AM frequency" field.
        %fname = fieldnames(datain);
        fname = fieldList_Freq;
        % Get the list of valid expressions in the field_fn.
        field_fn_exprs = regexp( field_fns{fi},'\w+','match');
        % See if any of those match with the fieldnames (perfectly)
        fieldchk = [];
        for rexpi = 1:length(field_fn_exprs);
            tmpi = find(strcmp(field_fn_exprs{rexpi},fname));
            if ~isempty(tmpi)
                fieldchk = tmpi;
            end;
        end;
        if  length(fieldchk)>0 
            ind = findstr(field_fns{fi},fname{fieldchk});                
            if ind==1
                cmdstr = ['thisdata.' field_fns{fi}];
            else
                cmdstr = [field_fns{fi}(1:ind-1) 'thisdata.' field_fns{fi}(ind:end)];
            end;
            
            fcount = fcount+1;            
            fprintf('Testing: %s\n',cmdstr);
            if eval(['~isempty(thisdata.' fname{fieldchk} ')'])
                f(fcount,1:length(fullAMList)) = eval(cmdstr)';
            else
                f(fcount,1:length(fullAMList)) = 0;
            end;
        end;
    end;
end;

% Passes if all the inclusion criteria pass.
chki = find(prod(f,1));
amList = fullAMList(chki);





