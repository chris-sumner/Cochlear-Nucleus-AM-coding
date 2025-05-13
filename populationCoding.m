%% ---------------------------------------------------------------------
%                population coding

% Paths for spike train analysis code. Mex - requires compiling for correct processor.
addpath('functions\Wohlgemuth Ronacher\Matlab');
addpath('functions\Wohlgemuth Ronacher\C');

 % Load up the spike trains.
allSpikeTimes88and91 = loadSpikeTimes;

% N.B. It is important to use 'allamfreqs' for selection as 
%      'usedamfreqs' has been selected with additional criteria which 
%      should not apply here. 
default_criteria = { ...
    'fn','amToneDur>=100', ...
    'fn','carrFreq>3000', ...
    'fn','cf>3000', ...
    'fn','strcmp(amType,''AM'')', ...
    'fn','mean(diff(allamfreqs))<=100', ...
    'fn','depthMod==1', ...
    'fn','fullDataUnitIndex~=1526', ...     % Exlcude this PLN which has very noisy performance.
    };


%% ------------------------------------------------------------------------

% Looking at the best neurons by type of increasing population size.
% 
% Split by neuron type and level. 
%    - run that subset for each unit type separately. 
%    - run for an increasing population. 
% 
% Neurons can only be grouped if they share the set of frequencies. 
% Practical decision to be made about constraints on this to enable us to
% make clusters. 
% A minimum number of 8 (so minimum 800Hz) frequencies works well. 
% An upper limit of 9 (so 900Hz) means that it is a similar set for all 
% unit types. 

levels = [30,50,70];
typeList = {'ChS','ChT','PL','PLN','PBU','On'}

% Run for all main types and levels. Takes about X mins. 
clear growingPopnByTypeLvl;
for ti = 1:6
    for li = 1:3
        fprintf('------------- TYPE: %s, LEVEL:  %d  -------------\n', typeList{ti},levels(li));
        tmpstatsmat  = select_Datasets(statsmat, ...
            default_criteria{:}, ...
            'fn',['modLvl_rationalised==' num2str(levels(li))], ...
            'fn',['strcmp(rationalisedType,''' typeList{ti}   ''')']);
        clear options;
        options.sampleN_bins = [1:40];
        options.minNumFreqs = 8;
        options.combineoptions = {'amfn','allusedamfreqs<=900'}; 
        options.numFreqs4Mean = 6;
        options.measure = populationCodingMeasure;
        
        % Do this by adding units of decreasing individual performance.
        options.direction = 'descend';        
        tmpop = ...
            runGrowingPopn_WR(tmpstatsmat,unitVSoutputs, ...
                unitoutputs,allSpikeTimes88and91,options);
        growingPopnByTypeLvl(li,ti).best1st = tmpop;
        
        % Do this by adding units of increasing individual performance.
        options.direction = 'ascend';
        % This option makes sure that we take the SAME set as used for the 
        % descending set. 
        options.bestN = length(growingPopnByTypeLvl(li,ti).best1st.addedUnit);
        tmpop = ...
            runGrowingPopn_WR(tmpstatsmat,unitVSoutputs, ...
                unitoutputs,allSpikeTimes88and91,options);
        growingPopnByTypeLvl(li,ti).worst1st = tmpop;
        
    end;
end;
 
clear tmpop tmpstatsmat allSpikeTimes88and91 default_criteria;
%save('datafiles\popnAnalByTypeLvl_F1','growingPopnByTypeLvl','options','typeList','levels');
