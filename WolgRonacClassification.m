
%% ----------- Run the classifier  ------------------

% Paths for Rob's code. Mex - requires compiling for correct processor.
addpath('functions\Wohlgemuth Ronacher\Matlab');
addpath('functions\Wohlgemuth Ronacher\C');
addpath functions;

% Load up:

% Load up spike times. 
% Split across 4 files to fit on github.
allSpikeTimes88and91 = loadSpikeTimes;

% Loads unit info and other supporting data.
load('rawdata\unitList');

% Which measure to use to select the best classifier tau
%classifierSelectMeasure = 'dprime';      % Before 2nd revision. 
classifierSelectMeasure = 'softmax_Z';

%%

% -------------- SELECTING and indexing MTFs  ----------------

% There is one decision to be made here which is the range of
% frequencies submitted to the classifier - if not all of them. This is
% because this affects the classifier results at other frequencies.

% % indices of lists to keep: modFreq<=2000
keepConds = wholeList171212(:,find(strcmp(EXPLOGLIST,'modFreq')))<=2000; 
% Shorten the data to exclude MTFs we will not analyse. 
wholeList171212 = wholeList171212(keepConds,:);
allSpikeTimes88and91 = allSpikeTimes88and91(keepConds);
allAMtypes88and91 = allAMtypes88and91(keepConds);

% Dividing the data up - into over 3000 mtfs.
[newdatainds dividingdata] = divideData(wholeList171212,EXPLOGLIST,allAMtypes88and91);

tau_list = [1 2 5 10 20 50];   % List of time constants to use.

for thisdatastartind = 1:length(newdatainds)

    fprintf('\n\n--------------------------------------\nDataset:%d\n',thisdatastartind);
    
    % Find the rows this piece of data occupy.
    if thisdatastartind<length(newdatainds)
        endind = newdatainds(thisdatastartind+1)-1;
    else
        endind = size((newdatainds),1);
    end
    thisdatainds = [newdatainds(thisdatastartind):endind];
    pars = dividingdata(newdatainds(thisdatastartind),:);
    
    % Get all the spike trains.
    spktrainset = [ allSpikeTimes88and91{thisdatainds} ];
    % Cut down to within 20-100ms
    spktrainset = cellfun(@(x)  x(x>20 & x<=100),spktrainset,'UniformOutput',false);

    % Get the am frequencies, VS, Z scores
    amfreqs = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'modFreq')));
    VSlist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'VS1')));  %30);
    Rayleighlist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'Ray1')));  %32);
    Zlist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'newZs')));  %66);
    MLtest = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'newZaboveRefrac1')));  %70);
    SperP = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'spikesPerPeriod')));  %17);
    ref0s = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'newZaboveRefrac0')));  %69);
    cilist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'CI')));  %60);
    ciPlist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'CIpVal')));  %71);

    % Remove any conditions with no spikes in any sweep.
    spktrainset_chk = cellfun(@(x) ~isempty(x), spktrainset);
    conditionswithspikes = find(sum(spktrainset_chk,1));
    spktrainset = spktrainset(:,conditionswithspikes);
    
    if ~isempty(spktrainset)
        amfreqset = ones(size(spktrainset,1),1)*(amfreqs(conditionswithspikes)');
    else
        amfreqset = [];
    end;
    
    % Classifier parameters.
    clear classifier;
    classifier.iterations = 1000;
    classifier.store_distances = 0;

    clear wroutput;
    for taui = 1:length(tau_list)
        % Input should just be a list of spiketrains and list of conditions of the same dimensionality.
        classifier.tau = tau_list(taui);
        wroutput(taui) = WR_Classifier(spktrainset,amfreqset,classifier);
    end;

    % A load of useful information about the data.
    unitinfo.unitid = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'uniqueID')));  
    unitinfo.modLevel = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'modLevel')));   
    unitinfo.amToneDur = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'amToneDur')));  
    unitinfo.carrFreq =  wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'carrFreq')));   
    unitinfo.depthMod =  wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'depthMod')));  
    unitinfo.cf =  wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'cf')));   
    unitinfo.cf_thr =  wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'cf_thr')));   
    unitinfo.allamfreqs =  amfreqs;
    unitinfo.usedamfreqs =  amfreqs(conditionswithspikes);

    unitinfo.TypeNum = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'TypeNum')));  %3);
    unitinfo.TypeName = typeNames{wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'TypeNum')))+1};
    unitinfo.VS = VSlist;
    unitinfo.Rayleigh = Rayleighlist;
    unitinfo.Z = Zlist;
    unitinfo.ref0s = ref0s;
    unitinfo.ci = cilist;
    unitinfo.cipVal = ciPlist;
    unitinfo.MLPoissTest = MLtest;
    unitinfo.CV = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'newCV')));  %63);
    unitinfo.spikesperperiod = SperP;
    unitinfo.conditionswithspikes = conditionswithspikes;
    unitinfo.classifierTauList = tau_list;
    
    % This keeps the entire set of statistics for the FIRST condition.
    % Useful for double checking data unit identity later. 
    unitinfo.scholesWholeListFirstEntry = wholeList171212(newdatainds(thisdatastartind),:);
    
    % Store all the results.
    unitoutputs(thisdatastartind).wroutput = wroutput;
    unitoutputs(thisdatastartind).unitinfo = unitinfo;

    % Code to plot output for every dataset
    % Set to "save" - also requires a folder to put the figures in. 
    if isfield(wroutput(1),'stimulusconditions')
          [filename] = plotWRClassifierOverallPerf(unitoutputs(thisdatastartind),thisdatastartind,'dontsave');
          unitoutputs(thisdatastartind).unitinfo.figName = filename;
          close(gcf);
    else
          % If there is an empty set of data (there is one)
          % make an empty set of unitinfo. 
          unit1unitinfonames = fieldnames(unitoutputs(1).unitinfo);
          fn2 = cell(2,length(unit1unitinfonames));
          fn2(1,:) = unit1unitinfonames';
          unitoutputs(thisdatastartind).unitinfo = struct(fn2{:});
    end;
    
    % Save every 50 units.
    if mod(thisdatastartind,50)==0
        save  datafiles\WR_results_tmp unitoutputs;
    end;

    close all;
end;
save datafiles\WR_results_tmp unitoutputs;

% Apply a set of additional measures. 
% This was added in response to one of the reviewers' comments. 
for di = 1:length(unitoutputs)
    fprintf('%d....\n',di);
    nwroutputs = length(unitoutputs(di).wroutput);
    for ci = 1:nwroutputs
        wroutputtmp(ci) = WR_ConfusionMatrixStats( unitoutputs(di).wroutput(ci) );
    end;         
    unitoutputs(di).wroutput =  wroutputtmp;        
    clear wroutputtmp;
end;

% ----------------  Apply softmax anaysis ------------------------------

% This is a method of summarising the confusion matrix which was 
% developed for the 2nd revision of the paper. 

figure; % We plot the fit for one classifier tau for each dataset.  
% Fit softmax model.
wsmooth = 0.0001;  % This is a very weak constraint. 
wnoise = 0.0001;   % This amount of noise to the starting parameters works well to stabilise the fit.
niter = 10;        % Number of times to fit.  
Zmax = 8;          % Final limit on Z becuase it makes no sense to have very high values.  
plims = [0 .99];   % Similarly, proportion correct >0.99 is meaningless in these data. 

for di = 1:length(unitoutputs)
    fprintf('%d....\n',di);
    nwroutputs = length(unitoutputs(di).wroutput);
    if isfield(unitoutputs(di).wroutput,'stimulusset')
        for ci = 1:nwroutputs
            wroutput_in = unitoutputs(di).wroutput(ci);
            % This does the fitting...      
            [softmax_Z softmax_B softmax_err softmax_fit] = fitSoftMax2ConfMat(wroutput_in.confusionmatrix,wsmooth,wnoise,niter,plims,Zmax);

            wroutput_in.softmax_Z = softmax_Z;  % N.B. In the code the c' statistic is called Z. Easy to confuse with Z_isi!  
            wroutput_in.softmax_B = softmax_B;
            % Store the fit, but remove the predicted confusion matrices,
            % which are easy to regenerate. 
            wroutput_in.softmax_fit = rmfield(softmax_fit  ,{'pmodel_init','pmodel_final'} );
            % Give us something to look at. 
            if ci == 1
               close(gcf);
               plotSoftMaxFit2ConfMat(softmax_fit,wroutput_in(1).confusionmatrix);
               set(gcf,'name',['Dataset:' num2str(di)]);
               drawnow;
            end;
            wroutput_out(ci) =  wroutput_in;
        end;
       unitoutputs(di).wroutput =  wroutput_out;             
    end;       
    clear wroutput_out;
end;

% Select out the tau that gives the best overall performance.
for i=1:length(unitoutputs);
    if ~isempty(unitoutputs(i).wroutput) && isfield(unitoutputs(i).wroutput(1),'dprime')
        unitdprimes = reshape([unitoutputs(i).wroutput(:).(classifierSelectMeasure)],length(unitoutputs(i).unitinfo.usedamfreqs),6);
        unitotaldprimes = sum(unitdprimes);
        [maxdprime bestind] = max(unitotaldprimes);
        unitoutputs(i).bestClassifier.index = bestind;
        unitoutputs(i).bestClassifier.tau = unitoutputs(i).wroutput(bestind).classifier.tau;
    end;
end;

save datafiles\WR_results_tmp unitoutputs;
%save datafiles\WR_results unitoutputs;   % This would overwrite the supplied file. 


for i=1:length(unitoutputs)
    unitop;
    
    if ~isempty(unitop.wroutput) && isfield(unitop.wroutput(1),'dprime')       
        for wi = 1:length(unitop.wroutput)
            tmpwrs(wi) = rmfield(unitoutputs(i).wroutput(wi),fields2rm);
        end;
        unitoutputs(i).wroutput(wi) = tmpwrs;
        clear tmpwrs;
    end;
end;


% Make the fle smaller for github. No reason to do this if you are
% re-running it.
% fields2rm = {'dprime_mafc_5','dprime_mafc_6','dprime_mafc_7','dprime_mafc_f1_5','softmax_fit'}
% morefields = [fields2rm,'stimulusset']
% for di = 1:length(unitoutputs)
%     fprintf('%d....\n',di);
%     nwroutputs = length(unitoutputs(di).wroutput);
%     if ~isempty(unitoutputs(di).wroutput) && isfield(unitoutputs(di).wroutput(1),'dprime') 
%         for ci = 1:nwroutputs
%             if isfield(unitoutputs(di).wroutput(ci),'stimulusset')                
%                 wroutputtmp(ci) = rmfield(unitoutputs(di).wroutput(ci),morefields);           
%             else
%                 wroutputtmp(ci) = rmfield(unitoutputs(di).wroutput(ci),fields2rm);           
%             end;               
%         end;         
%      unitoutputs(di).wroutput =  wroutputtmp;        
%     clear wroutputtmp;  
%     end;
% 
% end;
% %save datafiles\WR_results_tmp unitoutputs;
% save datafiles\WR_results unitoutputs;
% N.B. There is one supplemental figure that will not work now. 



