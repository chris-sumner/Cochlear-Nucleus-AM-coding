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
% unit types. Only PLN drops to 800Hz max at the highest level. 
% Making the maximum frequency 900Hz for most of the neurons, 
% and 800Hz for the PLNs should make the SU performance monotonic.

levels = [30,50,70];
typeList = {'ChS','ChT','PL','PLN','PBU','On'}

% Run for all main types and levels. Takes about X mins. 
clear growingPopnByTypeLvl;
for ti = 1:6
    for li = 1:3
        tmpstatsmat  = select_Datasets(statsmat, ...
            default_criteria{:}, ...
            'fn',['modLvl_rationalised==' num2str(levels(li))], ...
            'fn',['strcmp(rationalisedType,''' typeList{ti}   ''')']);
        clear options;
        options.sampleN_bins = [1:40];
        options.minNumFreqs = 8;
        options.combineoptions = {'amfn','allusedamfreqs<=900'}; 
        options.numFreqs4Mean = 6;
        
        % Do this by adding units of decreasing individual performance.
        options.direction = 'descend';        
        tmpop = ...
            runGrowingPopn_WR(tmpstatsmat,unitVSoutputs, ...
                unitoutputs,allSpikeTimes88and91,options);
        growingPopnByTypeLvl(li,ti).best1st = tmpop;
        
        % Do this by adding units of increasing individual performance.
        options.direction = 'ascend';
        % This option makes sure that we take the SAME set as used for the 
        % descending set. Can remove the real duffers. 
        options.bestN = length(growingPopnByTypeLvl(li,ti).best1st.addedUnit);
        tmpop = ...
            runGrowingPopn_WR(tmpstatsmat,unitVSoutputs, ...
                unitoutputs,allSpikeTimes88and91,options);
        growingPopnByTypeLvl(li,ti).worst1st = tmpop;
        
    end;
end;
 
save datafiles\popnAnalByTypeLvl;


%% -------------- Figure  ---------------

popnSize2Disp = 10; % The number of neurons for the population mtfs.
extraPopnSize2Disp = 30; % An additional one for the highest sound level.
                         % SHows how this differs between types. 
growXLim = [35 20 35 42];                         

typesToDisplay = {'ChS','ChT','PL','PLN'}
figh = figure('position',[100 100 1200 800],'paperposition',[.5 .5 24 16]);
typeColorMap;

for ci = 1:length(typesToDisplay)

    ti = find(strcmp(typeList, typesToDisplay{ci}));  % 1; % Type index.
    colorlist = min(1,[1.2 0.8 0.5]'*findRatTypeCol(typeList{ti}));  % Get the list of colours.

    mtf_h(ci) = subplot('position',[0.05+0.22*(ci-1) 0.63 0.18 0.27 ]);
    ngrow_h(ci) = subplot('position',[0.05+0.22*(ci-1) 0.2 0.18 0.27 ]);

    for li =1:3
        totalN = length(growingPopnByTypeLvl(li,ti).best1st.addedUnit);    
        totalNworst = length(growingPopnByTypeLvl(li,ti).worst1st.addedUnit);    

        % Some analysis of these results.     
        fullDataIndex_addedUnit = ...
            cellfun(@(x)  x.WRi{1}.fullDataMUIndexes, growingPopnByTypeLvl(li,ti).best1st.addedUnit(1:totalN));
        meanDPrime_addedUnit = ...
            cellfun(@(x) mean(x.WRi{1}.dprime(1:options.numFreqs4Mean)), growingPopnByTypeLvl(li,ti).best1st.addedUnit(1:totalN));
        meanDPrime_eachUnit = ...
            cellfun(@(x) mean(x.WRi{1}.dprime(1:options.numFreqs4Mean)),  ...
            growingPopnByTypeLvl(li,ti).best1st.eachUnit(  growingPopnByTypeLvl(li,ti).best1st.options.sampleOrder ));
        meanDPrime_addedUnit_rev = ...
            cellfun(@(x) mean(x.WRi{1}.dprime(1:options.numFreqs4Mean)), growingPopnByTypeLvl(li,ti).worst1st.addedUnit(1:totalNworst));
        meanDPrime_eachUnit_rev = ...
            cellfun(@(x) mean(x.WRi{1}.dprime(1:options.numFreqs4Mean)),  ...
            growingPopnByTypeLvl(li,ti).worst1st.eachUnit(  growingPopnByTypeLvl(li,ti).worst1st.options.sampleOrder ));
        fullDataMUIndex_lastUnit = ...
            cellfun(@(x) x.WRi{1}.fullDataMUIndexes(end), growingPopnByTypeLvl(li,ti).best1st.popns(1:totalN))
        nofFreqs_popn = ...
            cellfun(@(x)  length(x.WRi{1}.stimulusconditions), growingPopnByTypeLvl(li,ti).best1st.popns(1:totalN));
        maxFreq_popn = ...
            cellfun(@(x)  max(x.WRi{1}.stimulusconditions), growingPopnByTypeLvl(li,ti).best1st.popns(1:totalN));
        meanDPrime_popn = ...
            cellfun(@(x)  mean(x.WRi{1}.dprime(1:options.numFreqs4Mean)), growingPopnByTypeLvl(li,ti).best1st.popns(1:totalN));
        meanDPrime_popn_rev = ...
            cellfun(@(x)  mean(x.WRi{1}.dprime(1:options.numFreqs4Mean)), growingPopnByTypeLvl(li,ti).worst1st.popns(1:totalNworst));

        subplot(mtf_h(ci)); hold on;
        tmp_i = find(growingPopnByTypeLvl(li,ti).best1st.sampleN_bins(1:totalN)==popnSize2Disp);

        cellfun(@(x) plot( x.WRi{1}.stimulusconditions, x.WRi{1}.dprime, 'color',colorlist(li,:),'linewidth',1), ...
                    growingPopnByTypeLvl(li,ti).best1st.popns(tmp_i))
        
        if li == 3
          tmp_ii = find(growingPopnByTypeLvl(li,ti).best1st.sampleN_bins(1:totalN)==extraPopnSize2Disp);
          if ~isempty(tmp_ii)
                cellfun(@(x) plot( x.WRi{1}.stimulusconditions, x.WRi{1}.dprime,'--','color',colorlist(li,:),'linewidth',1), ...
                    growingPopnByTypeLvl(li,ti).best1st.popns(tmp_ii));
                        % PLot the sound level close to the line.
            stimf = growingPopnByTypeLvl(li,ti).best1st.popns{tmp_ii}.WRi{1}.stimulusconditions;                
            xperf = growingPopnByTypeLvl(li,ti).best1st.popns{tmp_ii}.WRi{1}.dprime;
            text(stimf(li)-30,xperf(1+li),['n=' num2str(extraPopnSize2Disp)],'fontweight','bold','color',colorlist(li,:));

           end;
        end;
                
        ylim([0 7]); xlim([50 950]); box off;
        xlabel('Modulation frequency (Hz)','fontweight','bold');
        ylabel('d''','fontweight','bold');
        if li == 1
            text(700,6.5,typeList{ti},'fontweight','bold','fontsize',14);
        end;

        % PLot the sound level close to the line.
        stimf = growingPopnByTypeLvl(li,ti).best1st.popns{tmp_i}.WRi{1}.stimulusconditions;                
        xperf = growingPopnByTypeLvl(li,ti).best1st.popns{tmp_i}.WRi{1}.dprime;
        text(stimf(2+li)-30,xperf(2+li),[num2str(levels(li))],'fontweight','bold','color',colorlist(li,:));

        % Growth of performance with the the number of neurons. 
        subplot(ngrow_h(ci));      
        ngrow = growingPopnByTypeLvl(li,ti).best1st.sampleN_bins(1:totalN);
        plot(ngrow,meanDPrime_popn(1:totalN), 'color',colorlist(li,:) ,'linewidth',1); hold on; 
        plot(ngrow,meanDPrime_eachUnit(1:totalN),':', 'color',colorlist(li,:),'linewidth',1);

        % PLot the sound level close to the line.
        text(ngrow(5+2*li)-.5,meanDPrime_popn(5+2*li),[num2str(levels(li))],'fontweight','bold','color',colorlist(li,:));

        xlabel('Number of neurons','fontweight','bold');
        ylabel('Ensemble d''','fontweight','bold');
        ylim([0 7]); xlim([0 growXLim(ci)]); box off;
        if li == 1
            text(7/9.5*max(xlim),6.5,typeList{ti},'fontweight','bold','fontsize',14);
        end;
               
    end;        
end;

subplot(  mtf_h(1) );
text(0,8,'a. Modulation transfer functions for small neural populations (n=10 or 30)','fontweight','bold','fontsize',10);
subplot(  ngrow_h(1) );
text(0,8,'b. Average performance for f_m_o_ds up to 600Hz as a function of population size','fontweight','bold','fontsize',10);

%print -dtiff -r600 figures\Figure8.tif
%hgsave(figh,'figures\Figure8');







