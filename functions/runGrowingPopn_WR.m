function output = runGrowingPopn_WR(statsmat,unitCIandVSoutputs,unitoutputs,allSpikeTimes88and91,options)
% Takes a population by type and runs it in order of decreasing individual performance.

% This probably needs to be specified for each level and type.
%options.sampleN_bins = [1 2:2:10 15:5:40 ];  

% --------------------- Run each unit separately ----------------------

% Doing a good job on this is complex. It should be based on the performance 
% of individual neurons when tested with a comparable set of frequencies.
% This is different to the original analysis, where the classifier was run
% on all available frequencies.

% Therefore we run the individual neurons again, restricting the modulation
% frequencies to 100Hz steps at 600Hz or below.
% We run this for ALL the neurons individually. 

if ~isfield(options,'direction')
    options.direction = 'descend';
end;

options1u = options;  % Take any inherited options.
options1u.sampleN = 1;
options1u.sampleIter = 1;
options1u.sampleMode = 'one';
%options1u.combineoptions = {'amfn','allamfreqs<=900'}; 


for ui=1:length(statsmat.fullDataUnitIndex)
   fprintf('\n--------------------- UNIT %d ----------------------\n',ui); 
   
   fdi_tmp= statsmat.fullDataUnitIndex{ui};
   statsmat_tmp = select_Datasets(statsmat, ...
       'fn',['fullDataUnitIndex == ' num2str(fdi_tmp) ]  );
   if ~isempty(statsmat_tmp)
       tmpop = runPopn_WR(statsmat_tmp,unitCIandVSoutputs, ...
                                 unitoutputs,allSpikeTimes88and91,options1u);

                             
       % Store the d-primes for sorting. 
       if isfield(tmpop,'WRi')
            % This is how we specify the measure used to order the 
            % population.
            dprimes_standardset{ui} =  tmpop.WRi{1}.(options.measure); %dprime;       
            
           % Store that output for returning.
           output.eachUnit{ui}.WRi{1} = tmpop.WRi{1};
           % Also store the unitinfo (for checking etc.).
           output.eachUnit{ui}.unitinfo = unitoutputs(fdi_tmp).unitinfo;
           output.eachUnit{ui}.unitinfo.fullDataUnitIndex = fdi_tmp;

       else 
           dprimes_standardset{ui} =  nan;
           % Store that output for returning.
           output.eachUnit{ui}.WRi{1} = [];
           % Also store the unitinfo (for checking etc.).
           output.eachUnit{ui}.unitinfo = unitoutputs(fdi_tmp).unitinfo;
           output.eachUnit{ui}.unitinfo.fullDataUnitIndex = fdi_tmp;
       end;
   end;
end;
clear statsmat_tmp tmpop;

% --------------- Sort out the ordering ------------------

% Check which units have enough modulation frequencies. 
% This is the indexes of the units that do. 
if isfield(options,'minNumFreqs')
    unitsWithNufFreqs = find(cellfun(@(x) length(x)>options.minNumFreqs, dprimes_standardset));
else
    unitsWithNufFreqs = find(ones(size(dprimes_standardset)));
end;

% Calculate the  order in which to run.
meandprimes(unitsWithNufFreqs) = cellfun(@(x) mean(x(1:min(length(x),options.numFreqs4Mean))),dprimes_standardset(unitsWithNufFreqs));
[meanDPrimes, sampleOrder] = sort(meandprimes,options.direction);

% Remove any that did not satisfy the frequency condition.
idx2keep = find( ismember(sampleOrder,unitsWithNufFreqs) );
meanDPrimes = meanDPrimes(idx2keep);
sampleOrder = sampleOrder(idx2keep);

% In some cases we also take the best N neurons. This is used for the
% ascending case. 
if isfield(options,'bestN')
    if length(sampleOrder)>options.bestN
        sampleOrder = sampleOrder(1:options.bestN);
        meanDPrimes = meanDPrimes(1:options.bestN);
    end;
end;

% ------------------- Run the populations -------------------

% Loop for different size populations. 
for si  = 1:length(options.sampleN_bins)
        
    % ---------------- Run this population -------------------
    
    % Options - ordered sampling of first N units in order.
    % Order is set by the performance of the individual units as above.
    optionsMU = options;                                    % Take any inherited options.
    optionsMU.sampleN = options.sampleN_bins(si);           % How many neurons to select.
    optionsMU.sampleIter = 1;                               % Number of sample iterations.
    optionsMU.sampleMode = 'ordered';                       % Number of sample iterations
    optionsMU.sampleOrder = sampleOrder;                     % sampling: best first.    
    
    % Run across all unit types.
    [output.popns{si}, tmpstats, tmpips] = runPopn_WR(statsmat,unitCIandVSoutputs, ...
        unitoutputs,allSpikeTimes88and91,optionsMU);
    
    % --------------- Re-run the last unit added ----------------
    % for the range of frequencies tested in the population.
    % so it can be compared with the population. 
    % This is largely redundant. 
    if isfield(output.popns{si},'WRi')

        ui = sampleOrder(options.sampleN_bins(si));
        fdi_tmp= statsmat.fullDataUnitIndex{ui};
        statsmat_tmp = select_Datasets(statsmat, ...
           'fn',['fullDataUnitIndex == ' num2str(fdi_tmp) ]  );

       % Construct the combine option to select the correct frequencies. 
       ftxt = sprintf('%d,', output.popns{si}.WRi{1}.stimulusconditions);
       cmdstr = ['find(ismember(allusedamfreqs,['  ftxt(1:end-1) ']))'];   
       options1u.combineoptions = {'amfn',cmdstr};
                
       tmpop = runPopn_WR(statsmat_tmp,unitCIandVSoutputs, ...
                                 unitoutputs,allSpikeTimes88and91,options1u);
                             
       % Store that output for returning.
       output.addedUnit{si}.WRi{1} = tmpop.WRi{1};
       % Also store the unitinfo (for checking etc.).
       output.addedUnit{si}.unitinfo = unitoutputs(fdi_tmp).unitinfo;
       output.addedUnit{si}.unitinfo.fullDataUnitIndex = fdi_tmp;
        
    end;
    
end;

output.sampleN_bins = options.sampleN_bins;
output.options = options;
output.options.sampleOrder = sampleOrder; 









