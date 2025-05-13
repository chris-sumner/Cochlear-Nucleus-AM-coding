function wrop2 =  WR_ConfusionMatrixStats(wroutput)
% WR_ConfusionMatrixStats - adds more statistics of the confusion matrix.
%
% wrop2 =  WR_ConfusionMatrixStats(wroutput)
%      
% Adds:    
%   hits
%   falsealarms
%   misses
%   correct rejections
% 

% We are just adding some more statistics to the structure. 
wrop2 = wroutput;

% check if there is any data...
if ~isfield(wroutput,'stimulusset')
    fprintf('No confusion matrix to analyse.\n');
    return;
end;

% From WR_classifier:
%  confusionmatrix - matrix of how many times each spike train of a given class j
%   (column index) was decided to belong to class i (row index). 

% ------- construct hit rate, fa and dprime -------

% This is based on code from WR_classifier.
uniquestimset = unique(wroutput.stimulusset(:));
noof_uniquestimuli = length(uniquestimset);
[sortedstimset sortedstimorder] = sort(wroutput.stimulusset(:));
sweepsperstimulus = hist(wroutput.stimulusset(:),uniquestimset);

confmat = wroutput.confusionmatrix;

hits =  sum(confmat.*eye(noof_uniquestimuli),1);
falsealarms =  sum(confmat,2) - hits';

% N.B. Hit rate needs to be corrected for comparisons that were not
% possible because they were with the same spike train.
hr_correctionfactor = sweepsperstimulus./(sweepsperstimulus-1);

% Calculate hit rate.
hitrate =  ( hits./ (sum(confmat)) );

% Calculate the false alarm rate. 
noof_nontargets = sum(confmat(:)) - sum(confmat);
farate = falsealarms./ noof_nontargets';

% New code starts here. 
misses = sum(confmat,1) - hits;       % Number of times that stimulus was tested when a different choice was made.
correctrejections = sum(confmat(:)) - sum(confmat,1) - sum(confmat,2)';
                                      % Number of times a stimulus was not
                                      % presented, and not chosen. 

precision = hits./(hits + falsealarms');
recall = hits./(hits + misses);
f1score = 2 * precision .* recall./(precision + recall);                                      

% Copy statistics to structure. 
wrop2.hits = hits;
wrop2.falsealarms = falsealarms';
wrop2.falsealarmrate = farate';
wrop2.misses = misses;
wrop2.correctrejections = correctrejections;
wrop2.precision = precision;
wrop2.recall = recall;
wrop2.hitrate = hitrate;
wrop2.f1score = f1score;
wrop2.classaccuracy = (hits + correctrejections)./(hits + correctrejections + falsealarms' + misses);
wrop2.overallaccuracy = sum(hits)/sum(confmat(:));
