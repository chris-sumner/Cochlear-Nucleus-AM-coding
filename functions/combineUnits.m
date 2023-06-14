function [popSpkTSets, unitinfoout, units2keep] = combineUnits(spkTSets,unitinfo,varargin)
% combineUnits combines units into a homogeneous set. 

nunits = length(unitinfo); 
amfield = 'usedamfreqs';            % Determines which am field we use. 

% ------------ Checking that the units match each other ------------

% Collect the levels,depth and duration - to checthte stimulus 
% conditions are the same. 

% If using rationalised values.
if sum(strcmp(varargin,'userationalised'))>0
    rLvls = [30 50 70];    
    fprintf('Rationalising sound levels to 30,50,70dB SPL\n');
    modLevels =  unique(arrayfun(@(x) max(rLvls(findnearest(x.modLevel,rLvls))),unitinfo));
else
    modLevels = unique( arrayfun(@(x) x.modLevel, unitinfo) );    
end;

depthMod = unique( arrayfun(@(x) x.depthMod, unitinfo) );
amToneDur = unique( arrayfun(@(x) x.amToneDur, unitinfo) );

% If any of these have more than one then they do not match. 
errorCode = false;
if  sum(strcmp(varargin,'ignorelevel'))==0 &&  length(modLevels) >1  
    fprintf('combineUnits: sound levels do not match\n');
    errorCode = true; 
end;
 
if  length(depthMod)>1 
    fprintf('combineUnits: modulation depths do not match\n');
    errorCode = true; 
end;

% if length(amToneDur)>1 
%     fprintf('combineUnits: tone durations do not match\n');
%     errorCode = true; 
% end;
    
% Next check that no unit is used more than once    
unitid =  arrayfun(@(x) x.unitid, unitinfo) ;
uniqueunitid =  unique(unitid);

if sum(strcmp(varargin,'ignoreunits'))==0  && length(unitid) > length(uniqueunitid)
    fprintf('combineUnits: units are not unique\n');
    errorCode = true; 
end;   

if errorCode
    error('Stimulus conditions do not match, or neurons are not unique.\nAnalysis stopped.');
end;


% ------------- Processing the frequencies -----------

%     % Rationalised AM frequencies - rounded to 100Hz.
%     statsmat.allamfreqs_rationalised{i} =  100*ceil(statsmat.allamfreqs{i}/100);
%     statsmat.usedamfreqs_rationalised{i} =  100*ceil(statsmat.usedamfreqs{i}/100);
%     statsmat.usedamfreqs_dprime_rationalised{i} =  100*ceil(statsmat.usedamfreqs_dprime{i}/100);

usedamfreqs = arrayfun(@(x) x.(amfield), unitinfo,'Uni',false);
if sum(strcmp(varargin,'userationalised')>0)
    fprintf('Rationalising modulation frequency to multiples of 100Hz\n');
    usedamfreqs = cellfun(@(x) 100*ceil(x/100), usedamfreqs,'Uni',false);
end;    
allusedamfreqs = unique([ usedamfreqs{:} ]);

% Check for additional constraints on which frequencies to use. 
% Pretty crude - you need to use a command which contains "allusedamfreqs"
amselectflag = find(strcmp(varargin,'amfn'));
if ~isempty(amselectflag)
    amselectfn = varargin{amselectflag+1};
    amselectinds = eval(amselectfn);    
    allusedamfreqs = allusedamfreqs(amselectinds);
end;

freqmatches = cellfun(@(x) ismember(allusedamfreqs,x)',usedamfreqs,'Uni',false);
freqmatches = [freqmatches{:}]';

commonfreqmatches = prod(freqmatches,1);
commonfreqs = allusedamfreqs(find(commonfreqmatches));

fprintf('Common AM frequencies:\n\t');
fprintf('%d,',commonfreqs);
fprintf('\n');

% Select out the common frequencies.
for i=1:length(unitinfo)    
    finds = find(ismember( ceil(unitinfo(i).(amfield)), commonfreqs));
    if length(finds)~=length(commonfreqs)
        finds = find(ismember( unitinfo(i).(amfield)+50, commonfreqs));
    end;
    if length(finds)~=length(commonfreqs)
        finds = find(ismember( 100*ceil(unitinfo(i).(amfield)/100), commonfreqs));
    end;
    popSpkTSets{i} = spkTSets{i}(:,finds);
 
    nSweeps{i} = size(popSpkTSets{i},1);
end;
    
% Compare the number of sweeps. 
uniqueNSweeps = unique([nSweeps{:}]);
maxNSweeps = max(uniqueNSweeps);

% If they do not match -  remove those that do not. 
% Remove units with less than the maximum. 
if length(uniqueNSweeps)>1
    fprintf('Not all units have the same number of sweeps');    
    units2remove = find([nSweeps{:}] ~= maxNSweeps);
    fprintf('Removing the %dth unit\n',units2remove);
    units2keep = find([nSweeps{:}] == maxNSweeps);
        
    popSpkTSets = popSpkTSets(units2keep);
    unitinfo = unitinfo(units2keep);
else
    units2keep = [1:length(unitinfo)];
end;


% Make a new unitinfo that encompasses all the data. 
unitinfoout.unitid = arrayfun(@(x) x.unitid, unitinfo);
unitinfoout.modLevels = modLevels;
unitinfoout.depthMod = depthMod;
unitinfoout.amToneDur = amToneDur;
unitinfoout.carrFreq = arrayfun(@(x) x.carrFreq, unitinfo);
unitinfoout.cf = arrayfun(@(x) x.carrFreq, unitinfo);
unitinfoout.cf_thr = arrayfun(@(x) x.cf_thr, unitinfo);
unitinfoout.allamfreqs = commonfreqs;
unitinfoout.usedamfreqs = commonfreqs;
unitinfoout.TypeNum = arrayfun(@(x) x.TypeNum, unitinfo);
unitinfoout.TypeName = arrayfun(@(x) x.TypeName, unitinfo,'Uni',false);
unitinfoout.wholeListInds = arrayfun(@(x) x.wholeListInds, unitinfo,'Uni',false);
unitinfoout.scholesWholeListFirstEntry = arrayfun(@(x) x.scholesWholeListFirstEntry, ...
    unitinfo,'Uni',false);


















