function [sals_sorted peaks_sorted troughs_sorted] = computesaliences(x,peaklist,troflist)
% Function to compute the complete list of saliences.

% Remove any nans and sort the peaks and troughs. 
speaklist = sort(peaklist(~isnan(peaklist)));
stroflist = sort(troflist(~isnan(troflist)));
% Common minimum length 
minlen = min(length(speaklist),length(stroflist));

% The amplitudes of the peaks and troughs.
xpeaks = x(speaklist);
xtrofs = x(stroflist);

% If the peaks start before troughs, shove them 
% along one to change that. 
if min(speaklist)<min(stroflist)
    xtrofs = [nan xtrofs];
end;
   
% Make sure stroflist is one longer than the speaklist
xtrofs = [xtrofs nan(1,(length(xpeaks)-length(xtrofs) +1))];

% Compute the saliences on either side. 
sal_b4 = xpeaks - xtrofs(1:end-1);   
sal_after = xpeaks - xtrofs(2:end);

% Overall minimum salience for each pair. 
saliences = min([sal_b4;sal_after],[],1);

% Make sure that the lists out are all the same length
% and that the saliences correspond correctly with the peaks.
% This is done differently for peaks and troughs because 
% salience is defined for peaks, not troughs, and the
% presence of nans in the salience makes a mess of it. 

% If there are more peaks than troughs    
% This bit is done BEFORE the resorting. 
if length(stroflist)<length(speaklist)
    if min(speaklist)>min(stroflist)
        % If the sequence starts with a peak
        % add nan to the beginning of the troughs
        stroflist = [nan stroflist];
    else
        % Otherwise add nan to the end.
        stroflist = [stroflist nan];
    end;
end;

% Work out the order of the saliences, 
% BEFORE adding nans to pad. 
[sals_sorted sal_sortedi] = sort(saliences,2,'descend');


% Re-order peaks and troughs
peaks_sorted = speaklist(sal_sortedi);
troughs_sorted = stroflist(sal_sortedi);

% If there are more troughs than peaks
% then AFTER sorting add a nan to pad 
% salience and peaks. 
if length(stroflist)>length(speaklist)
    peaks_sorted = [peaks_sorted nan];
    sals_sorted = [sals_sorted nan];    
    troughs_sorted = [troughs_sorted stroflist(~ismember([1:length(stroflist)],sal_sortedi))];    
end;



