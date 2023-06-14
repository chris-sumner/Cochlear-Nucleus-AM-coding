% 9/9/17. This is trying to figure out the way that Scholes selected out
% which data to analyse.
% This was pressing given that he had clearly included data which was 200%
% modulated in his analysis. 

%% ------------------------- loading data --------------------------

% Path required for some functions.
addpath functions;

% Adds root path to any files.
bngFn = @(s) ['D:\bng\' s];   
%drivefn = @(s) ['V:' s];

% Load up:
load(bngFn('rhode_modelocking\fromrhode\cn88data\justSpikesTimes88and91'),'allSpikeTimes88and91');    % Load up the spike trains.
load(bngFn('rhode_modelocking\fromrhode\cn88data\justUnitLists88and91PBUsorted'),'wholeList171212');  % A table of the data
load(bngFn('rhode_modelocking\fromrhode\cn88data\justUnitLists88and91PBUsorted'),'EXPLOGLIST');       % A key to that table,
load(bngFn('rhode_modelocking\fromrhode\cn88data\justUnitLists88and91PBUsorted'),'typeNames');        % type names.
load(bngFn('rhode_modelocking\fromrhode\cn88data\justUnitLists88and91PBUsorted'),'allAMtypes88and91');% AM type



%% ------------------------- Selecting data ---------------------------

% Logical values - which principally select out units of 30dB above
% threshold.
load(bngFn('rhode_modelocking\fromrhode\cn88data\justUnitLists88and91PBUsorted'),'bestCondInds');
load(bngFn('rhode_modelocking\fromrhode\cn88data\justUnitLists88and91PBUsorted'),'bestCondIndsAllFreqs');
%bestCondInds = bestCondIndsAllFreqs;

% However... these two sets are different. Scholes switched to using
% bestCondIndsAllFreqs at some point but it's not clear why.

% This was code to figure out the difference between bestCondInds & bestCondIndsAllFreqs
indsunits = unique(wholeList171212(logical(bestCondInds),strcmp(EXPLOGLIST,'uniqueID')))
indsunitsallf = unique(wholeList171212(logical(bestCondIndsAllFreqs),strcmp(EXPLOGLIST,'uniqueID')))
% They are different lists of units. 195 vs. 359.

% The key difference is in the carrier - the extra units did not have a
% carrier frequency which I think indicates they had a noise carrier. 
notinboth = find(~ismember( indsunitsallf,indsunits ))
bob = wholeList171212(notinboth,:)
hist(bob(:,find(strcmp(EXPLOGLIST,'carrFreq'))))
hist(wholeList171212(:,find(strcmp(EXPLOGLIST,'carrFreq'))),100) 

% The code which *might* have prodcued this, taken from 'SfNAbsDataAnal' is
% here encapuslated in a function, which returns a similar array. 
condsNr30dBSL = selectDataByThreshold(wholeList171212,30);
% However, this yields even more units!
newindsunits = unique(wholeList171212(logical(condsNr30dBSL),strcmp(EXPLOGLIST,'uniqueID')))

% So exactly how these lists are produced remains a mystery. It is possible
% that additional selection criteria have been applied


%% --------------- Criteria for selecting data --------------------

% These two functions summarise the various ranges of variables and
% quantities we care about.
dataSummary(wholeList171212,EXPLOGLIST);
dataSummary(allAMtypes88and91);

% These lead one to a number of simple ways the data should be pared before
% even starting any analysis:
% amToneDur   - nearly all 100ms. Little point analysing the others.
% depthMod    - mostly 1.0 (33k) 2.0 (20k). Two populations probably worth
%               analysing separately
% AM type     - many more plain AM (47k) than the other 3 (and I don't
%               know what they are). 
% modFreq     - actually a very large range of these - up to 5k but many up
%               to 2k. I think there is a problem with the different range
%               for PLs. 
% carrFreq    - > 3000 nans. Not clear what this means (noise carrier?). 
%               We don't want to analyse anything with carrier frequency or
%               CF <3k. 
% cf & cf_thr - >1500 nans. Not clear what this means. Surely cannot
%              analysse them?

% These things may underlie the differences in Chris' dataset selections.

% We need a standard repeatable way of controlling which data we look at.
% - A function with several ways of extracting subsets of data.
%      : 30dB above threshold.
%      : rejecting a standard set of unanalysable things as per above.
%      : selecting out specific things such as a range of sound levels.

% But then: How do these variables depend upon each other?
% We have always ignored high modulation frequencies as they are applied
% very unevenly across the datasets. 
% - Look at the range of these values across different unit types. 
% This is the only criteria that really matters for how we run the
% classifier. 

% We need a robust way of dividing up data into units for the Wohlemuth
% classifier (check existing method).

% It appears that sometimes the same unit is unit with the same conditions
% multiple times. 

% This shows that the range of f mods is actually moderately similar across
% the major unit types.
typeinds = @(type) find( wholeList171212(:,find(strcmp(EXPLOGLIST,'TypeNum')))  == find(strcmp(typeNames,type))-1 )  ;
dataSummary(wholeList171212(typeinds('ChS'),:),EXPLOGLIST);
dataSummary(wholeList171212(typeinds('ChT'),:),EXPLOGLIST);
dataSummary(wholeList171212(typeinds('PL'),:),EXPLOGLIST);
dataSummary(wholeList171212(typeinds('PLN'),:),EXPLOGLIST);
dataSummary(wholeList171212(typeinds('PBU'),:),EXPLOGLIST);


