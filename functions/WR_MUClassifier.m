function unitoutput = WR_MUClassifier(spktrainsets,stimulusset,classifieroptions)
% WR_MUClassifier   MU spike train classification based on Wohlegmuth and Ronacher 
%
% SYNOPSIS:
% Wohlegmuth and Ronacher (2007; J. Neurphysiol. 97:3082-3092) used this method
% to classify neuronal responses to AM of different frequency. It uses a spike
% distance metric similar van Rossum, except that it is based on alpha function
% windows instead of exponential. The classification algorithm is a single trial
% template match. For each iteration of the classifier one spike train is picked 
% for each class (e.g. AM frequency) at random. All remaining spike trains are 
% then classed as belonging to the template which is most similar. 
%
% THe MU version takes in a set of spiketrains from different units. It
% computes the Euclidean distance across units (sqrt of sum of squares). 
%
% USE:
% classdata =  WR_MUClassifier(spiketrains,stimulusset[,options])
%    classdata: output data structure containing the results of the classification.
%         wr_spktrdist - NxN symmetrical matrix of spike-train pair distances.
%         stimulusconditions - list of all the classes, in order of 
%               increasing value.
%         sweepspercondition - list of the number of spike trains for each class.  
%         stimulusset - the full list of which class each spike train belongs to. 
%         confusionmatrix - matrix of how many times each spike train of a given class j
%                (column index) was decided to belong to class i (row index). 
%         hitrate - proportion of true positive decisions for each class.
%         falsealarmrate - proportion of false positive decisions for each class.
%         dprime - d' (signal detection measure of sensitivity) for each class
%             Z(hitrate) - Z(falsealarmrate)
%         classifier - contains all the parameters used by the classifier:
%             iterations, store_distances, tau and alpha
%    spiketrains: a cell array or matrix of N spike trains. 
%    stimulusset: numeric matrix of class values (or designations) for each spike train.
%    options: a structure of options for the classifier which can contain any of:  
%         iterations - number of classifier iterations (default 100)
%         store_distances - flag of whether to output wr_spktrdist.
%         tau - time constant for the spike distance window (2.45/alpha).
%         alpha - decay constant for the spike distance window (default 1).
%
% Spike distance code: wr_metric, written by Rob Mill (Dec, 2012). Uses an
% analytic solution to efficiently calculate the W&R spike distance metric.
%
% C.J. Sumner and R.Mill Dec 2012. 



% Paths for Rob's code. Mex - requires compiling for correct processor. 
% addpath 'U:\robertm\Wohlgemuth Ronacher\current-temp\Matlab';
% addpath 'U:\robertm\Wohlgemuth Ronacher\current-temp\C';

% Check that the dimensionality of the stimulus set and spike train set are
% identical
if size(spktrainsets{1})~=size(stimulusset)
    error('WR_Classifier: Spike train and stimulus set must be the same size');
end;

% ---- Parse the classifier parameters ----
% Default iterations
classifier.iterations = 100;
classifier.sweeppermutations = 50;

% Add in any options, overriding the existing ones.
optionfields = fieldnames(classifieroptions);
for i=1:length(optionfields)
    classifier.(optionfields{i}) = classifieroptions.(optionfields{i});
end;

% Sets a flag to disable the sweep permuting. 
if classifier.sweeppermutations == 0;
    classifier.sweeppermute = false;
end;
if isfield(classifier,'sweeppermute') && ~classifier.sweeppermute
    classifier.sweeppermutations = 1;
end;
    
   
% Set up so that the coincidence window can be specified as either alpha or
% tau.
if ~isfield(classifier,'alpha') & ~isfield(classifier,'tau')
    classifier.alpha = 1;
    classifier.tau = 2.45./classifier.alpha;
elseif ~isfield(classifier,'alpha')
    classifier.alpha = 2.45./classifier.tau;
elseif ~isfield(classifier,'tau')
    classifier.tau = 2.45./classifier.alpha;
end;    

% Display some info about the classifier parameters:
fprintf('Iterations: %d\nTau:',classifier.iterations);
fprintf('%d,',classifier.tau);
fprintf('\nalpha:');
fprintf('%g,',classifier.alpha);
fprintf('\n');

% Output is initialised to just the classifier specs. 
unitoutput.classifier = classifier;
	
% ----------- work out structure for the random templates --------

% Group together conditions which are similar. 
uniquestimset = unique(stimulusset(:));
noof_uniquestimuli = length(uniquestimset);
[sortedstimset sortedstimorder] = sort(stimulusset(:));
sweepsperstimulus = hist(stimulusset(:),uniquestimset);


% If there are no sweeps then exit.
if isempty(sweepsperstimulus) | sum(sweepsperstimulus==0)~=0
    fprintf('Some conditions with no sweeps! Not running.\n');
    return;
end;

fprintf('# Spiketrains for each neuron: %d\n',length(stimulusset(:)));
fprintf('# neurons: %d\n',length(spktrainsets(:)));
fprintf('Stimulus conditions:'); fprintf('%.4g ',uniquestimset); fprintf('\n');
fprintf('# sweeps of each:   '); fprintf('%.4g ',sweepsperstimulus); fprintf('\n');
	
% ------- calculate the raw spike distances -------

% This calculates the raw spike distances between every spike train. 
% This is done in the 'sorted stimulus' order. 

% Initialise (huge) empty matrix for all the units separately. 
wr_MUspktrdist = zeros(length(sortedstimorder),length(sortedstimorder),length(spktrainsets));

% Preallocate the confusion matrix. 
confmat = zeros(noof_uniquestimuli); % preallocate

% Do this for every unit. 
% In order that the result does not reflect some oddity of the way
% sweeps are ordered (it shouldn't anyway). Calculate spike train
% distances multiple times whilst permutating sweeps differently
% between neurons. 
for si = 1:classifier.sweeppermutations    
    
    if classifier.sweeppermute    
        fprintf('Sweep permutation %d\n',si);
    else
        fprintf('Sweep permute is off\n');
    end;
    
    % Set the spike train distances to zero. 
    wr_MUspktrdist(:) = 0;
    
    for ui=1:length(spktrainsets)          
        % Randomly permute the sweeps for each unit. 
        % We assume sweep order does not need to be preserved. 
        if classifier.sweeppermute
            spktrainsets{ui} = permuteEachColumn(spktrainsets{ui});
            
            % Oversampling?            
        end;            
        % Calculate the spike train distance for a unit. 
        tmp_wr_spktrdist = wr_metric(spktrainsets{ui}(sortedstimorder),classifier.alpha(ui),'norm');
        % Put nans down the diagonal. 
        tmp_wr_spktrdist(find(eye(size(tmp_wr_spktrdist)))) = nan;

        % Distances for all units are stored here. 
        wr_MUspktrdist(:,:,ui) =  tmp_wr_spktrdist;        
    end;        
        
    % ANY WEIGHT OPTIMISATION WOULD GO HERE.
    % PERHAPS LIMIT THE NUMBER OF ITERATIONS FOR THE FIT.
    % EACH ITERATION, BASED ON A DIFFERENT PERMUTATION OF SPIKES.

    % random subsampling should be done here. 
    
    if length(size(wr_MUspktrdist))>2
        wr_spktrdist = sqrt(sum( wr_MUspktrdist.^2, 3));
    else
        wr_spktrdist = wr_MUspktrdist;
    end;
    
    % ---- perform the iterative template matching -----

    % Each 'class' will have 'iterations' of template choices, taken from the
    % full set of options. This generates the correct number of random numbers and then 
    % multiply them by the correct factors and round up. 
    % This the indexes telling for each iteration what the set of templates is. 
    templateset = ones(classifier.iterations,1)* (cumsum(sweepsperstimulus)- sweepsperstimulus(1)) + ...
        ceil( (ones(classifier.iterations,1)*sweepsperstimulus) ...
        .*rand(classifier.iterations,noof_uniquestimuli));
    
    % Iterate around and find the smallest distance in each case.
    stimuluschoice = zeros(classifier.iterations,size(wr_spktrdist,1));

    for i=1:classifier.iterations
        % Select out the spikedistances for the selected templates vs. every other spiketrain. 
        iterationdistances = wr_spktrdist(templateset(i,:),:);

        % For each MU spike train find the nearest template.     
        [mindist minind] = nanmin(iterationdistances);        
        stimuluschoice(i,:) = uniquestimset(minind)';

    end;

    % -------  Confusion matrix generation -----

    %confmat = zeros(noof_uniquestimuli); % preallocate
    for i=1:size(stimuluschoice,2)    
        stimuluscondition = find(uniquestimset==sortedstimset(i));  % Find the sweeps for that condition
        confmat(:,stimuluscondition) = confmat(:,stimuluscondition) + ...
            hist(stimuluschoice(:,i),uniquestimset)';
    end;
    
end;

% As the confustion matrix is the sum of multiple permuations of the
% sweeps from each neruon...
confmat = confmat/classifier.sweeppermutations;


% ------- construct hit rate, fa and dprime -------

% For every template - one spike train cannot be compared.
% Templates are chosen at random from the set of each stimuli. 
% Thus (1/noof sweeps) for that stimulus condition cannot be matched as the distance is nan.
% Therefore the hitrate must be adjusted by 1/noof_sweeps.

hits =  sum(confmat.*eye(noof_uniquestimuli),1);
falsealarms =  sum(confmat,2) - hits';

% N.B. Hit rate needs to be corrected for comparisons that were not
% possible because they were with the same spike train.
hr_correctionfactor = sweepsperstimulus./(sweepsperstimulus-1);

% Calculate hit rate.
hitrate =  hr_correctionfactor.*( hits./ (sum(confmat)) );

% Calculate the false alarm rate. 
noof_nontargets = sum(confmat(:)) - sum(confmat);
farate = falsealarms./ noof_nontargets';

dprime = zlookup(hitrate) - zlookup(farate)';
%pcmax = plookup(dprime);

% have not computed dprime.
mean_dprime = mean(dprime);

% Also compute MI?



% ------------- output to data structure ---------
if classifier.store_distances
    unitoutput.wr_spktrdist = wr_spktrdist;
end;
unitoutput.stimulusconditions = uniquestimset;
unitoutput.sweepspercondition = sweepsperstimulus;
unitoutput.stimulusset = stimulusset;
unitoutput.confusionmatrix = confmat;
unitoutput.hitrate = hitrate;
unitoutput.falsealarmrate = farate;
unitoutput.dprime = dprime;
unitoutput.mean_dprime = mean_dprime;
%unitoutput.pcmax = pcmax;

fprintf('Finished.\n');



% ---- z and p lookup - just copied from log2_detection ----

function zscore = zlookup(pvalue)

p = [0.999:-0.001:0.991  0.99:-0.01:0.01 0.009:-0.001:0.001 ];

zp = [3.090 2.878 2.748 2.652 2.576     2.512 2.457 2.409 2.366         2.326 2.054 1.881 1.751 1.645 ...
      1.555 1.476 1.405 1.341 1.282     1.227 1.175 1.126 1.080 1.036    0.994 0.954 0.915 0.878 0.842 ...
      0.806 0.772 0.739 0.706 0.674     0.643 0.613 0.583 0.553 0.524   0.496 0.468 0.440 0.412 0.385 ...
      0.358 0.332 0.305 0.279 0.253     0.228 0.202 0.176 0.151 0.126   0.100 0.075 0.050 0.025 0.000];
zp = [zp -fliplr(zp(1:end-1))];  

zscore = nan*zeros(size(pvalue));

for i=1:length(pvalue(:))
    pdiff = abs(p - pvalue(i));
	zscoretmp = zp(find(pdiff == min(pdiff)));
    if ~isempty(zscoretmp)
    	zscore(i) = zscoretmp(1);
    else
        zscore(i) = nan;
    end;
end;

function pvalue = plookup(zscore)

p = [0.999:-0.001:0.991  0.99:-0.01:0.01 0.009:-0.001:0.001 ];

zp = [3.090 2.878 2.748 2.652 2.576     2.512 2.457 2.409 2.366         2.326 2.054 1.881 1.751 1.645 ...
      1.555 1.476 1.405 1.341 1.282     1.227 1.175 1.126 1.080 1.036    0.994 0.954 0.915 0.878 0.842 ...
      0.806 0.772 0.739 0.706 0.674     0.643 0.613 0.583 0.553 0.524   0.496 0.468 0.440 0.412 0.385 ...
      0.358 0.332 0.305 0.279 0.253     0.228 0.202 0.176 0.151 0.126   0.100 0.075 0.050 0.025 0.000];
zp = [zp -fliplr(zp(1:end-1))];

for i=1:length(zscore)
    zdiff = abs(zp - zscore(i));
    pvaluetmp = p(find(zdiff == min(zdiff)));

    if ~isempty(pvaluetmp)
        pvalue(i) = pvaluetmp(1);
    else 
        pvalue(i) = nan;
    end;
end;


function spkset = permuteEachColumn(spktrainset)

setsize = size(spktrainset);
spkset = cell(setsize);

for j=1:setsize(2)    
   spkset(:,j) = spktrainset(randperm(setsize(1)),j);
end;
    
    







