function stats = summaryDatasets(datain,varargin)
% Extract statistics about the dataset.
% Look for duplicate sets.
% Display nicely


stats.unitlist = unique([datain.unitid{:}]);
stats.unitTable = cell2table({int32(0) nan nan '' nan nan nan nan});
stats.unitTable.Properties.VariableNames = {'UnitID','CF','CFthreshold','Type','CV', ...
                                            'NoofDatasets','NoofLevels','NoofModDepths'};

% Option to select which am frequency field to use.                                         
if nargin == 1
    amfreqstr = 'allamfreqs';
else
    amfreqstr = varargin{1};
end;
                                        
% Loop through every unit.
for ui=1:length(stats.unitlist)

    uinds = find([datain.unitid{:}]==stats.unitlist(ui));
        
    % Get the different conditions.
    modLevel =  [datain.modLevel{uinds}]';
    carrFreq =  [datain.carrFreq{uinds}]';
    amType =    cellfun(@(x) x{1},datain.amType(uinds),'uni',false)';
    amToneDur = [datain.amToneDur{uinds}]';
    depthMod =  [datain.depthMod{uinds}]';

    % Modulation frequency statistics.
    allamfreqs = datain.(amfreqstr)(uinds)';
    allamfreqs(find(cellfun(@(x) isempty(x),allamfreqs))) = {nan};
    
    minModFreq = cellfun(@(x) min(x), allamfreqs );
    maxModFreq = cellfun(@(x) max(x), allamfreqs );
    numModFreq = cellfun(@(x) length(x), allamfreqs );
    stepModFreq = cellfun(@(x) mean(diff(x)), allamfreqs );
    
    % Make a table of the conditions for that unit.
    stats.conditionTable{ui} = table(amType,carrFreq,amToneDur,depthMod,modLevel, ...
        minModFreq,maxModFreq,numModFreq,stepModFreq);
    
    
    % Check for duplicate rows....
    unqConds = unique(stats.conditionTable{ui});
    if size(unqConds,1)<size(stats.conditionTable{ui})
        fprintf('Unit %d has duplicate datasets\n',stats.unitlist(ui));
        display(stats.conditionTable{ui});
    end;
    
    % Find all characteristics for that unit.
    cf = unique( [datain.cf{uinds}] );
    notnancfs = find(~isnan(cf));
    cf = cf(notnancfs);
    if isempty(cf)
        cf = nan;
    end;
    if length(cf)>1
        ambig_flag = 1;
        fprintf('WARNING unit %d has more than one CF\n',stats.unitlist(ui));
    end
       
    
    cf_thr = unique( [datain.cf_thr{uinds}] );
    notnan = find(~isnan(cf_thr));
    cf_thr = cf_thr(notnan);
    if isempty(cf_thr)
        cf_thr = nan;
    end;
    if length(cf_thr)>1
        ambig_flag = 1;
        fprintf('WARNING unit %d has more than one CF thr\n',stats.unitlist(ui));
    end;

    type = unique( datain.type(uinds) );
    if length(type)>1
        ambig_flag = 1;
        fprintf('WARNING unit %d has more than one type\n',stats.unitlist(ui));
    end;

    CV = unique( [datain.CV{uinds}] );
    if length(CV)>1
        ambig_flag = 1;
        fprintf('WARNING unit %d has more than one CV value:\n\t',int32(stats.unitlist(ui)));
        fprintf('%g, ',CV);
        fprintf(' Taking (nan)mean.\n');
        CV = nanmean(CV);
    end;
    
    noofDatasets = length(uinds);
    noofLevels = length(unique(modLevel));
    noofModDepths = length(unique(depthMod));
         
    % Make the unit table.. 
    stats.unitTable(ui,:) = {stats.unitlist(ui) cf cf_thr type CV noofDatasets noofLevels noofModDepths};
        
end;
    
fprintf('\nData summary:\n');
fprintf('\t%d units\n',length(stats.unitlist));

stats.typelist = unique( stats.unitTable.Type );
for ti = 1:length(stats.typelist)
    % indexes to units of that type:
    inds  = find(arrayfun(@(x) x==1, strcmp(stats.unitTable.Type,stats.typelist{ti} )));
    fprintf('\n\nType: %s, n=%d units\n',stats.typelist{ti},length(inds));
    
    % CF range:
    fred = stats.unitTable.CF(inds);
    fprintf('\tCFs: %g - %g Hz, mean= %g Hz\n',min(fred),max(fred),mean(fred));
   
    % Threshold range:
    % CF range:
    fred = stats.unitTable.CFthreshold(inds);
    fprintf('\tCF thresholds: %g - %g (dB), mean= %g (dB)\n',min(fred),max(fred),mean(fred));
    
    % Threshold range:
    % CF range:
    fred = stats.unitTable.CV(inds);
    fprintf('\tCV: %g - %g , mean= %g \n',min(fred),max(fred),mean(fred));

    % Carrier range:
    fred = cellfun(@(x) x.carrFreq',stats.conditionTable(inds),'uni',false);
    fredlist = [fred{:}];
    fprintf('\tCarrier Frequencies: %g - %g Hz, mean= %g Hz\n',min(fredlist),max(fredlist),mean(fredlist));

    % Duation range:
    fred = cellfun(@(x) x.amToneDur',stats.conditionTable(inds),'uni',false);
    fredlist = [fred{:}];
    fprintf('\tAM Tone Durations: %g - %g ms, mean= %g ms\n',min(fredlist),max(fredlist),mean(fredlist));

    % Modulation frequencies.
    fprintf('\nModulation freqencies (%s):\n\t\t',amfreqstr);
    
    % Minimum
    fred = cellfun(@(x) x.minModFreq',stats.conditionTable(inds),'uni',false);
    fredlist = [fred{:}];
    fredlist = fredlist(~isnan(fredlist));
    unifred = unique(fredlist);
    fprintf('\n\tMinimums:\n\t\t');
    fprintf('%3g, ',unique(fredlist)); fprintf('\n\t\t');
    if length(unifred)==1
       fprintf('Always same value')
    else
       fredh = hist(fredlist,unifred);
       fprintf('%3g, ',fredh); 
    end;

     % Maximum
    fred = cellfun(@(x) x.maxModFreq',stats.conditionTable(inds),'uni',false);
    fredlist = [fred{:}];
    fredlist = fredlist(~isnan(fredlist));
    unifred = unique(fredlist);
    fprintf('\n\tMaximum:\n\t\t');
    fprintf('%3g, ',unique(fredlist)); fprintf('\n\t\t');
    if length(unifred)==1
       fprintf('Always same value')
    else
       fredh = hist(fredlist,unifred);
       fprintf('%3g, ',fredh); 
    end;

   
     % Number of modulation frequencies.
    fred = cellfun(@(x) x.numModFreq',stats.conditionTable(inds),'uni',false);
    fredlist = [fred{:}];
    unifred = unique(fredlist);
    fredlist = fredlist(~isnan(fredlist));
    fprintf('\n\tNumber of modulation frequencies:\n\t\t');
    fprintf('%3g, ',unique(fredlist)); fprintf('\n\t\t');
    if length(unifred)==1
       fprintf('Always same value')
    else
       fredh = hist(fredlist,unifred);
       fprintf('%3g, ',fredh); 
    end;

    
     % Modulation frequency step size.
    fred = cellfun(@(x) x.stepModFreq',stats.conditionTable(inds),'uni',false);
    fredlist = [fred{:}];
    fredlist = fredlist(~isnan(fredlist));
    unifred = unique(fredlist);
    fprintf('\n\tFrrequency step sizes:\n\t\t');
    fprintf('%3g, ',unique(fredlist)); fprintf('\n\t\t');
    if length(unifred)==1
       fprintf('Always same value')
    else
       fredh = hist(fredlist,unifred);
       fprintf('%3g, ',fredh); 
    end;
    
    
    % Modulation levels.
    fred = cellfun(@(x) x.modLevel',stats.conditionTable(inds),'uni',false);
    fredlist = [fred{:}];
    fredh = hist(fredlist,unique(fredlist));
    nfreds = cellfun(@(x) length(unique(x)), fred);
    nh = hist(nfreds,[0:10]);
    fprintf('\n\tModulation levels:\n\t\t');
    fprintf('%3g, ',unique(fredlist)); fprintf('\n\t\t');
    fprintf('%3g, ',fredh); 
    fprintf('\n\tNumber of levels for each unit:\n\t');
    fprintf('\t%g',[0:10]);fprintf('\n\t');
    fprintf('\t%g',nh);fprintf('\n');

    % Modulation depths.
    fred = cellfun(@(x) x.depthMod',stats.conditionTable(inds),'uni',false);
    fredlist = [fred{:}];
    nfreds = cellfun(@(x) length(unique(x)), fred);
    nh = hist(nfreds,[0:10]);
    fprintf('\n\tModulation depths:\n\t\t');
    fprintf('%g, ',unique(fredlist));
    fprintf('\n\tNumber of depths for each unit:\n\t')
    fprintf('\t%g',[0:10]);fprintf('\n\t');
    fprintf('\t%g',nh);fprintf('\n');
    
end;

fprintf('\n');
% Number of duplicate datasets.
% (display duplicate datasets).

% Summary of number of units of each type.
unitCounts = [];
for ti = 1:length(stats.typelist)
    % indexes to units of that type:
    inds  = find(arrayfun(@(x) x==1, strcmp(stats.unitTable.Type,stats.typelist{ti} )));
    fprintf('\nType: %s, n=%d units',stats.typelist{ti},length(inds));
    unitCounts(ti) = length(inds);
end;
fprintf('\nTotal number of units: %d.\n',sum(unitCounts));
