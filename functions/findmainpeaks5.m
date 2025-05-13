function [peaksets, trofsets salsets] = findmainpeaks5(x)
% findmainpeaks5
% A reasonably robust aplgorithm for picking sets of peaks and troughs
% from a vector. Designed for shuffled autocorrelation. 
%
% It returns sets, of increasing numbers of peaks and troughs. At every
% addition it seeks to add a single peak and trough which have the maximum 
% 'salience' of all the choices. 

global PLOTOPT;

peaklist = peakpick(x);
troflist = trofpick(x);

if isempty(peaklist) || isempty(troflist)
    peaksets = {};
    trofsets = {};
    salsets = {};    
    return;
end;
    
% Find the biggest peak:
[maxpeak maxpeaki] = max(x(peaklist));
% Find the lowest minima:
[mintrof mintrofi] = min(x(troflist));

% The first solultion is the global maximum and minimum troughs.
peaksets{1} =  peaklist(maxpeaki); 
trofsets{1} =  troflist(mintrofi);
salsets{1} = x(peaklist(maxpeaki))-x(troflist(mintrofi));
minsals(1) = salsets{1};

% Loop around until no more peaks are found.
iter = 1; itermore = true;
setn = 2;
while itermore
    
    [peakspickedout troughspickedout saliencesout] = peakandtroughdivider(x,peaklist,troflist, ...
        sort(peaksets{setn-1}),sort(trofsets{setn-1}));
    %N.B. peaks and troughs are sorted into ascending order for passing
    %back into peakandtroughdivider. 
    
    oldsetn = setn-1;    
    if ~isempty(peakspickedout) && ~isempty(troughspickedout)
        for newseti = 1:length(peakspickedout)
            
            % Sort the new trough and peak sets in order. 
            newpeakset = sort([peaksets{oldsetn} peakspickedout(1:newseti)]);
            newtrofset = sort([trofsets{oldsetn} troughspickedout(1:newseti)]);
            
            % peakandtrough divider only computes the saliences of the new
            % peaks it has just created, not the 'new' saliences of the
            % other peaks with the new troughs.
            % So for this set recompute the saliences. 
            [sals_sorted peaks_sorted trofs_sorted] =computesaliences(x,newpeakset,newtrofset);
                        
            % Store the set, in order of decreasing salience.
            salsets{setn} = sals_sorted;
            peaksets{setn} = peaks_sorted;
            trofsets{setn} = trofs_sorted;    
            % N.B. order of troughs is rather arbitrary since a salience is defined for a peak not a trough.
            minsals(setn) = nanmin( salsets{setn} );

           setn = setn+1;
         end;             
         iter = iter+1;
    else
        itermore = false;
    end;     
end;

% Remove any nans (not sure why they would be there).
peaksets = cellfun(@(x) x(~isnan(x)), peaksets,'uni',false);
trofsets = cellfun(@(x) x(~isnan(x)), trofsets,'uni',false);
salsets = cellfun(@(x) x(~isnan(x)), salsets,'uni',false);


function [peakspickedout troughspickedout saliencesout] = peakandtroughdivider(x,peaklist,troflist,peakspickedin,trofspickedin)
% peakandtroughdivider
% Takes a list of peaks and troughs already picked, and the complete list
% of minima and maxima, and adds a peak and a trough (or trough and peak)
% between each pair. Thus is approximately doubles the number of peaks and
% troughs. 

% Remove any nans
peakspickedin = peakspickedin(~isnan(peakspickedin));
trofspickedin = trofspickedin(~isnan(trofspickedin));

%fprintf('Sub-dividing peaks and troughs...\n');
% Make the total list of troughs and peaks. 
peaksandtroughsin  = sort([peakspickedin trofspickedin]);

% work out if this is a peak or a trough.
peaki = find(peakspickedin == peaksandtroughsin(1));
trofi =  find(trofspickedin == peaksandtroughsin(1));

% --- The first in the list is a special case ---
% Goes from the beginning of the series to the first peak or trough.

% Construct the sublist of peaks and troughs
subtroflist = troflist(troflist<=peaksandtroughsin(1));
subpeaklist = peaklist(peaklist<=peaksandtroughsin(1));   

% If this is the next thing is a peak...
if peakspickedin(1)<trofspickedin(1)
    [newpeaksout newtrofsout dummy1 dummy2 newsalsout] = splitwithtroughandpeak(x,subpeaklist,subtroflist,'beginning');
% If this is the next thing is a trough...
elseif peakspickedin(1)>trofspickedin(1)
    [newpeaksout newtrofsout dummy1 dummy2 newsalsout] = splitwithpeakandtrough(x,subpeaklist,subtroflist,'beginning');
end;
% N.B. whether the first thing is a trough or peak, which function
% you call is reversed compared with the other situations.


for i=1:length(peaksandtroughsin)-1
    newpeak = [];
    newtrof = [];
    newsal = [];
    
    % work out if this is a peak or a trough.
    peaki = find(peakspickedin == peaksandtroughsin(i));
    trofi =  find(trofspickedin == peaksandtroughsin(i));
    
    if ~isempty(peaki)                
        % Find the next picked trough.
        nextrough = min(trofspickedin(trofspickedin>peaksandtroughsin(i)));
        
        % Construct the sublist of peaks and troughs
        subtroflist = troflist(troflist>peaksandtroughsin(i) & troflist<=nextrough);
        subpeaklist = peaklist(peaklist>=peaksandtroughsin(i) & peaklist<=nextrough);   
        
        if ~isempty(subtroflist) & ~isempty(subpeaklist)
            % Subdivide the peak/trough section.
            [newpeak newtrof dummy1 dummy2 newsal] = splitwithpeakandtrough(x,subpeaklist,subtroflist);
        end;
        
        if newpeak<newtrof
            error('');
        end;
        
    elseif ~isempty(trofi)
        % Find the next picked peak.
        nextpeak = min(peakspickedin(peakspickedin>peaksandtroughsin(i)));
        
        % Construct the sublist of peaks and troughs
        subtroflist = troflist(troflist>=peaksandtroughsin(i) & troflist<=nextpeak);
        subpeaklist = peaklist(peaklist>=peaksandtroughsin(i) & peaklist<=nextpeak);   

        if ~isempty(subtroflist) & ~isempty(subpeaklist)
            % Subdivide the peak/trough section.
            [newpeak newtrof dummy1 dummy2 newsal] = splitwithtroughandpeak(x,subpeaklist,subtroflist);
        end;
        
        if newpeak>newtrof
            error('');
        end;
        
    end;
    
    if ~isempty(newpeak) && ~isempty(newtrof)
        if (isempty(newpeaksout) || isempty(find(newpeaksout==newpeak)) ) ... 
            & (isempty(newtrofsout) || isempty(find(newtrofsout==newtrof)))                     
            % Add them to the list of new peaks and troughs.
            newpeaksout = [newpeaksout newpeak];
            newtrofsout = [newtrofsout newtrof];
            newsalsout = [newsalsout newsal];       
        else
            error('Duplication of peaks or troughs should not happen');
        end;
    end;
end;

% ---- Special case: the end ----
% Considers the very last picked peak or trough to the end of the function

newpeak = [];
newtrof = [];
newsal = [];

% work out if this is a peak or a trough.
peaki = find(peakspickedin == peaksandtroughsin(end));
trofi =  find(trofspickedin == peaksandtroughsin(end));

% Construct the sublist of peaks and troughs
subtroflist = troflist(troflist>=peaksandtroughsin(end));
subpeaklist = peaklist(peaklist>=peaksandtroughsin(end));    

% If this is the next thing is a peak...
if ~isempty(peaki)
  [newpeak newtrof dummy1 dummy2 newsal] = splitwithpeakandtrough(x,subpeaklist,subtroflist,'ending');
% If this is the next thing is a trough...
elseif ~isempty(trofi)
  [newpeak newtrof dummy1 dummy2 newsal] = splitwithtroughandpeak(x,subpeaklist,subtroflist,'ending');
end;

% Add them to the list of new peaks and troughs.
% Note that this is at the end so it is possible to add only either a peak
% or trough.
if ~isempty(newpeak) && ~isempty(newtrof)
    if (isempty(newpeaksout) || isempty(find(newpeaksout==newpeak)) ) ... 
        & (isempty(newtrofsout) | isempty(find(newtrofsout==newtrof)))                     
        newpeaksout = [newpeaksout newpeak];
        newtrofsout = [newtrofsout newtrof];
        newsalsout = [newsalsout newsal];
    else
        error('Duplication of peaks or troughs should not happen');
    end;
end;

% --- Combine the new list of peaks with the old.
% N.B. Used to hand out the complete list, now only hand out the new ones.
% And in order of decreasing salience
if ~isempty(newsalsout)
    [saliencesout salorder] = sort(newsalsout,2,'descend');
    peakspickedout = newpeaksout(salorder(1));
    troughspickedout =newtrofsout(salorder(1));
    saliencesout = saliencesout(1);
else
    saliencesout =[];
    peakspickedout = [];
    troughspickedout = [];
end;



function [newpeak newtrof newpeakpos newtrofpos newminsalience] = splitwithpeakandtrough(x,peaklist,troflist,endcond)
% Function which aims to split a restricted bit of vector, which starts with a 
% peak and ends with a trough, adding a peak and a trough of the maximum salience. 

global CRITERION;

% Initialise outputs to nothing.
newpeak=[]; newtrof=[]; newpeakpos=[]; newtrofpos=[]; newminsalience= [];

if isempty(peaklist) | isempty(troflist)
    return;
end;

beginflag = 0; endflag = 0;
if nargin>3
    if strcmp(endcond,'beginning')
        beginflag = 1;
    elseif strcmp(endcond,'ending')
        endflag = 1;
    end;
end;
        
% Check that there are further peaks and troughs to find.
if (length(peaklist)<2 & length(troflist)<2) & ~beginflag & ~endflag
 %   fprintf('splitwithpeakandtrough: No peaks and troughs to find.\n');
    return;
end;

% If this is the beginning of the input, the first peak in the peaklist
% has not already been picked, and so there needs to be a dummy peak 
% at the start.
if beginflag
    peaklist = [nan peaklist];
    % Also, the first reversal could be a peak, in which case there also
    % needs to be an extra dummy trough so that the sequence is always
    % peak-trough-peak-trough... 
   if troflist(1)>peaklist(2)
        troflist = [nan troflist];
    end;
end;

% If this is the end of the input, then the last trough has not already been
% picked in which case there needs to be a final trough corresponding to the 
% one that would already have been picked.
if endflag 
    % If the last reversal is also a trough then you need to add a peak
    % as well so that the it ends trough-peak-trough.
    if troflist(end)>peaklist(end)
        peaklist = [peaklist nan];
    end;
    troflist = [troflist nan];
end;

% There are three saliences for every combination:
% peak-newtrough,newtrough-newpeak,newpeak-trough.
saliences = nan*zeros(length(troflist),length(peaklist),3);

% Compute all the peak-trough differences (saliences) between 
% Adjacent peaks and troughs for every combination of the choice
% of new peak and trough.
maxtrofi = length(troflist);

% This structure inherently assumes that there IS a peak and trough.
% However at a beginning and an end it is possible to find only one of
% those. 

% Loop through all the troughs:
for ti = 1:maxtrofi

    % Decide which peaks to loop through. 
    peakinds = ti:length(peaklist);    % This works and is needed but I'm not quite sure why.
    if ~isnan(troflist(ti))
        % Either peaks must follow this peak, or be an added 'nan'
        peakinds = find(peaklist>troflist(ti) | isnan(peaklist));
    end;
    
    % Loop through all the peaks: 
    for pi =peakinds
        
        % peak-newtrough
        if  ~isnan(peaklist(1)) && ~isnan(troflist(ti))  
        %if ~beginflag & peaklist(pi)>troflist(ti)            
            saliences(ti,pi,1) = x(peaklist(1)) - x(troflist(ti));
        end;

        % newpeak-newtrough
        % if peaklist(pi)>troflist(ti)        
        if  ~isnan(peaklist(pi)) && ~isnan(troflist(ti))
            saliences(ti,pi,2) = x(peaklist(pi)) - x(troflist(ti));
        end;
        
        % newpeak-trough
        % Only compute if the 
        %if ~endflag %& troflist(end)>peaklist(pi)
        if  ~isnan(peaklist(pi)) && ~isnan(troflist(end))            
            saliences(ti,pi,3) = x(peaklist(pi)) - x(troflist(end));        
        end;
    end;
end;   
  
% Find the minumum salience for each of the peak trough choices. 
minsaliences = nanmin(saliences,[],3);

% If there exists the option to pick other than trough alone,
% cut out the trough alone options. 
minsalience_notroughonly = minsaliences(:, ~isnan(peaklist));
if sum(~isnan(minsalience_notroughonly(:)))>0
    minsaliences(:,isnan(peaklist)) = nan;
end;

% Find the maximum of all those minimums: this is the new peak-trough
% combination that has the most prominent new peak.
[bestminsalience bestpeaktroughind] = nanmax(minsaliences(:));

% Remove any negatives - since they are not peak-trough pairs.
posinds = find(bestminsalience>CRITERION);
bestminsalience = bestminsalience(posinds);
bestpeaktroughind = bestpeaktroughind(posinds);

% Best min salience must be +ve or it this pair does not constutite a peak
% and trough.
if ~isempty(bestminsalience) % & bestminsalience>CRITERION
    % Find the position in the salence matrix of the most salient set of
    % peak and trough.
    [newtrofpos newpeakpos] = find(minsaliences == bestminsalience);

    % If the minimum salience is the same across a number of combinations.
    % This can happen if one minimum pair is shared across combinations. 
    if length(newtrofpos)>1
        % Compute the sum of the saliences. 
        sumsaliences = nansum(saliences,3);
        % Find the sums of saliences for the picked combinations.
        for i=1:length(newtrofpos)
            bestsumsaliences(i) = sumsaliences(newtrofpos(i),newpeakpos(i));
        end;
        [bestoverallsumsaliences bestoverallsumsaliencesi] = max(bestsumsaliences);
    
        newtrofpos = newtrofpos(bestoverallsumsaliencesi);
        newpeakpos = newpeakpos(bestoverallsumsaliencesi);
    end;
    
    newpeak = peaklist(newpeakpos);
    newminsalience = bestminsalience;
    newtrof = troflist(newtrofpos);
end;




function [newpeak newtrof newpeakpos newtrofpos newminsalience] = splitwithtroughandpeak(x,peaklist,troflist,endcond)
% Function which aims to split a restricted bit of vector which starts in a trough
% and ends in a peak, adding a peak and a trough of the maximum salience. 
% This is the opposite of the splittroughwithpeak code.

global CRITERION;

% Initialise outputs to nothing.
newpeak=[]; newtrof=[]; newpeakpos=[]; newtrofpos=[];newminsalience= [];

% If there are not peaks and troughs, there is nothing to do.
if isempty(peaklist) | isempty(troflist)
    return;
end;

% If this is not really trough-peak, but start-peak...
beginflag = 0; endflag = 0;
if nargin>3
    if strcmp(endcond,'beginning')
        beginflag = 1;
    elseif strcmp(endcond,'ending')
        endflag = 1;
    end;
end;

% Check that there are further peaks and troughs to find.
if length(peaklist)<2 & length(troflist)<2 & ~beginflag & ~endflag
  %  fprintf('splitwithpeakandtrough: No peaks and troughs to find.\n');
    return;
end;

% If this is the beginning of the input, the first trough in the peaklist
% has not already been picked, and so there needs to be a dummy trough 
% at the start.
if beginflag
    troflist = [nan troflist];
    % Also, the first reversal could be a trough, in which case there also
    % needs to be an extra dummy peak so that the sequence is always
    % trough-peak-trough-peak... 
   if troflist(2)<peaklist(1)
        peaklist = [nan peaklist];
    end;
end;

% If this is the end of the input, then the last peak has not already been
% picked in which case there needs to be a final peak corresponding to the 
% one already picked.
if endflag 
    % If the last reversal is also a peak then you need to add a trough
    % as well so that the it ends peak-trough-peak.
    if peaklist(end)>troflist(end)
        troflist = [troflist nan];
    end;
    peaklist = [peaklist nan];  
end;

% There are three saliences for every combination:
% newpeak-trough,newpeak-newtrough,peak-newtrough.
saliences = nan*zeros(length(peaklist),length(troflist),3);

maxpeaki = length(peaklist); % min(find(peaklist>troflist(ti)));

% Loop through all the troughs:
for pi = 1:maxpeaki
    % Loop through all the peaks: 
    trofinds = pi:length(troflist);
    if ~isnan(peaklist(pi))
        trofinds = find(troflist>peaklist(pi) | isnan(troflist));
    end;
        
     for ti =trofinds
        % newpeak-trough
        %if ~beginflag %& (troflist(ti)>peaklist(pi)
        if ~isnan(peaklist(pi)) && ~isnan(troflist(1))
            saliences(pi,ti,1) = x(peaklist(pi)) - x(troflist(1));
        end;

        % newpeak-newtrough
        %if troflist(ti)>peaklist(pi)
        if ~isnan(peaklist(pi)) && ~isnan(troflist(ti))
            saliences(pi,ti,2) = x(peaklist(pi)) - x(troflist(ti));
        end;
        
        % peak-newtrough
        %if ~endflag % & troflist(ti)>peaklist(pi)
        if ~isnan(peaklist(end)) & ~isnan(troflist(ti))
            saliences(pi,ti,3) = x(peaklist(end)) - x(troflist(ti));
        end;
        
    end;
end;   
  
% Find the minumum salience for each of the peak trough choices. 
minsaliences = nanmin(saliences,[],3);

% If there exists the option to pick other than trough alone,
% cut out the trough alone options. 
minsalience_notroughonly = minsaliences(~isnan(peaklist),:);
if sum(~isnan(minsalience_notroughonly(:)))>0
    minsaliences(isnan(peaklist),:) = nan;
end;

% Find the maximum of all those minimums: this is the new peak-trough
% combination that has the most prominent new peak.
[bestminsalience bestpeaktroughind] = nanmax(minsaliences(:));

% Remove any negatives - since they are not peak-trough pairs.
posinds = find(bestminsalience>CRITERION);
bestminsalience = bestminsalience(posinds);
bestpeaktroughind = bestpeaktroughind(posinds);


% Best min salience must be +ve or it this pair does not constutite a peak
% and trough.
if ~isempty(bestminsalience)
    [newpeakpos newtrofpos] = find(minsaliences == bestminsalience);
        
    % If the minimum salience is the same across a number of combinations.
    % This can happen if one minimum pair is shared across combinations. 
    if length(newtrofpos)>1
        % Compute the sum of the saliences. 
        sumsaliences = nansum(saliences,3);
        % Find the sums of saliences for the picked combinations.
        for i=1:length(newtrofpos)
            bestsumsaliences(i) = sumsaliences(newpeakpos(i),newtrofpos(i));
        end;
        [bestoverallsumsaliences bestoverallsumsaliencesi] = max(bestsumsaliences);
    
        newtrofpos = newtrofpos(bestoverallsumsaliencesi);
        newpeakpos = newpeakpos(bestoverallsumsaliencesi);
    end;    
    newpeak = peaklist(newpeakpos);

    newtrof = troflist(newtrofpos);
    newminsalience = bestminsalience;
end;










