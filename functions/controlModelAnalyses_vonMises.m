% ---------- Neural fluctuation --------------

% It was not clear what neural fluctuation measures. Specifically:
% What does it do as vector strength varies?

% This a simulation which mimicked neural responses by simulating phase locking
% using t von Mises distribution. It ran all the analysis
% methods. But we only looked at envelope fluctuation because that was the
% mystery. 

% It demonstrates essentially that the neural fluctuation measure rises 
% monotonically with VS. This is robust across different d.c. rates 
% (i.e. numbers of spikes). 

clear;

% Path required.
addpath functions;
addpath ..\SACcode;  % You need a path to the correlation index code. 

% load up the real data for reference
load('datafiles\SpikeStatsSets','unitVSoutputs');           % Load VS analysis.
load('datafiles\WR_results','unitoutputs');     % load WR classifier performance.

% Paths for Rob's code. Mex - requires compiling for correct processor.
addpath('functions\Wohlgemuth Ronacher\Matlab');
addpath('functions\Wohlgemuth Ronacher\C');

% The list of vector strengths we want to look at.
VS_list = [0:.1:1];
    
% Construct the input parameters for the control model. 
pars.amfreqs = [50:100:1050]';
pars.dcrates = [50:100:450];
pars.nsweeps =  100;                        % Sweeps higher than data for accurate measures. 
pars.dur_s = 0.08;
pars.classifier_amfreq_inds = [1:length( pars.amfreqs)];

% Which statistics to run
pars.calcGC = true;
pars.calcML = false;
pars.calcSACpeaks = false;

% Repeat for sytematically varying VS.
mi = 1;
for vsi = 1:length(VS_list)

    % Specify the VS and angle. Here we set them the same across all
    % frequencies.
    ang = pi + pars.amfreqs'*pi/30;  % We allow for some phase shift.
    vs =  VS_list(vsi)*ones(1,length( pars.amfreqs)); %  1  - pars.classifier_amfreq_inds/40;     
    pars.phaselocking = [vs' ang'];
    
    % And for varying dc rates. 
    for dci = 1:length(pars.dcrates)
        % Set the dc rate.
        pars.dcrate = pars.dcrates(dci)*ones(length(pars.amfreqs),1);
    
        % Find some classifier specs. Aribtarily take the first one in the dataset 
        bestInd = unitoutputs(1).bestClassifier.index;

        % Perform the simulation and all the calculations
        ctlModels{mi} = VSControlModel(pars,unitoutputs(1).wroutput(bestInd).classifier);
        ctlModels{mi}.specifiedPhaseLocking = pars.phaselocking;
        ctlModels{mi}.specifiedDCrate = pars.dcrate; 
        mi = mi+1; 
    end;
end;

% Reformat to extract out the envelope fluctuation and reliablity measures.
Twin = ones(length(VS_list)*length(pars.dcrates),1)*ctlModels{1}.GC{1}.temporalWindows;
envFluct_table = table;
reliability_table = table;

for ami = 1:length(pars.amfreqs )    

    tmp = cell2mat(cellfun(@(x) (x.GC{ami}.envFluct), ctlModels,'uni',false)');
    tmpvs = cell2mat(cellfun(@(x) (x.specifiedPhaseLocking(ami,1)), ctlModels,'uni',false)');
    tmpdc = cell2mat(cellfun(@(x) (x.specifiedDCrate(ami,1)), ctlModels,'uni',false)');

    VSmat = tmpvs*ones(1,length(ctlModels{1}.GC{1}.temporalWindows));
    DCmat = tmpdc*ones(1,length(ctlModels{1}.GC{1}.temporalWindows));

    tmptable = table(pars.amfreqs(ami)*ones(length(VS_list)*length(pars.dcrates)*length(ctlModels{1}.GC{1}.temporalWindows),1), ...
           DCmat(:),VSmat(:),Twin(:),tmp(:), ...
          'VariableNames',{'AM_freq','Spike_rate','VS','t_win','envFl'} );   
    %      'VariableNames',{'AM_freq','VS','t_win','envFl'} );   
    envFluct_table = [envFluct_table;tmptable];

    tmp = cell2mat(cellfun(@(x) (x.GC{ami}.reliability_R), ctlModels,'uni',false)');
    tmptable = table(pars.amfreqs(ami)*ones(length(VS_list)*length(pars.dcrates)*length(ctlModels{1}.GC{1}.temporalWindows),1), ...
           DCmat(:),VSmat(:),Twin(:),tmp(:), ...
          'VariableNames',{'AM_freq','Spike_rate','VS','t_win','reliability_R'} );   
    reliability_table = [reliability_table;tmptable];

end;

% Look at how envelope fluctuation depends on VS.
% It scales montonically with DC rate though.
figure; 
DC_choice = 250;
VS_selection = [0 0.4 0.7 1];
for pi = 1:4
    subplot(2,2,pi);
    tmptable = unstack(envFluct_table(envFluct_table.VS==VS_selection(pi) & envFluct_table.Spike_rate==DC_choice,:),'envFl','t_win');
    plot(tmptable.AM_freq, table2array( tmptable(:,4:end) ));
    title(['VS:' num2str(VS_selection(pi))]);
    xlabel('Modulation frequency');
    ylabel('Neural/Env fluctuation');
end;
legend(num2str(ctlModels{1}.GC{1}.temporalWindows'));
% Envelope fluctuation goes up monotonically with vector strength.
% It MUST also go up monotonically with extra peaks, so it is good
% describing the contributions of both, but cannot separate modes from
% phase locking. 

figure; 
plot(envFluct_table.Spike_rate, envFluct_table.envFl./envFluct_table.Spike_rate)

% Look at how reliablity depends on VS.
% This was not the point of the exercise so did not pursue it far. 
figure; 
DC_choice = 450;
VS_selection = [0 0.4 0.7 .9];
for pi = 1:4
    subplot(2,2,pi);
    tmptable = unstack(reliability_table(reliability_table.VS==VS_selection(pi) & ...
        reliability_table.Spike_rate==DC_choice,:),'reliability_R','t_win');
    plot(tmptable.AM_freq, table2array( tmptable(:,4:end) ));
    title(['VS:' num2str(VS_selection(pi))]);
    xlabel('Modulation frequency');
    ylabel('Reliability');
end;
legend(num2str(ctlModels{1}.GC{1}.temporalWindows'));
% This is a bit messier actually. 

% Remove unwanted variables.
clear tmpdc vs ami DCmat dci dcratemat dcrates tmp tmpdc tmpvs tmptable vsi VSmat ans mi p Twin;
clear unitoutputs unitVSoutputs;

% Save this.
%save datafiles\VMcontrolmodel;
