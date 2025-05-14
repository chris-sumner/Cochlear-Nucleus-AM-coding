%% ----------------------------------------------------------------------
%                   This the main analysis script 
%
% All the code is supplied to go from raw data (spike times) to figures. 
%
% To produce the intermediate data files required to run this script:
% 1. WolgRonacClassification - produces datafiles\WR_results.mat (~2 days to run). 
%        - No need to run this. File supplied on Github (22MB).
%        - Confusion matrices are generated via bootstrapping, and softmax 
%          statistic is a fit to this, so both will vary slightly if re-run.
% 2. SpikeStatsCalc - produces datafiles\SpikeStats.mat (1.3GB).
%        - Intermediate processing stage. Not required to run the code
%          below. This file is not stored on GitHub due to size, and takes a 
%          couple of days to generate (1.3GB). Output from (3) is provided. 
%        - Analysing the SAC peaks is part of this code but is set not to
%          run by default. It produces datafiles\SACpeaksA/B/C.mat (>100MB & 
%          additional weeks to run; split to fit on Github). Also not
%          needed for the code below
% 3. groupSpikeStatsByMTF - prodcues datafiles SpikeStatsSets.mat (700MB or 15MB).
%        - The data is reorganised into MTF "datasets" (SpikeStats.mat is
%          organised condition by condition), and economised slightly.
%        - Only takes a few minutes to run but requires SpikeStats.mat. 
%        - On github there is a smaller version, which reduces this to
%          16MB,by ruthlessly removing anything not needed subsequently.          
%

%% -----------------------------------------------------------------------
clear

% Load up and collate all the data into a format where each measure is an
% array of cells for each dataset.

% Path required.
addpath functions;

% Load analysis for everything except the classifier.
% Small version of file supplied on github.
load('datafiles\SpikeStatsSets_small','unitVSoutputs');     

% load WR classifier performance. File supplied on Github.
global unitoutputs
load('datafiles\WR_results','unitoutputs');    

% This integrates all the information into a single set of cell arrays.
% 'performance' parameter is what method of classifier scoring is used. 
% This is (due to changes in review), called 'dprime' (we used to use d'!). 
% In the early stags of processing it is referred to as Z . 
% In the paper is referred to as c', by analogy to d' and to avoid confusion with Z_isi.
% At this level in the code, "dprime" in statsmat is whatever your classifier statistic is 
% (so probably c').
statsmat = collate_byDataset(unitVSoutputs,unitoutputs,'performance','softmax_Z');

% This is for alternative measures of the classifier performance.  
statsmat_HR = collate_byDataset(unitVSoutputs,unitoutputs,'performance','hitrate');
statsmat_F1 = collate_byDataset(unitVSoutputs,unitoutputs,'performance','f1score');

% If you wish to replace a measure for all stages of analysis...
%statsmat = statsmat_F1;
%statsmat = statsmat_HR;
 
% Uncomment this to do the entire analysis with a specified resolution for
% the classifier. By default tau is chosen to opimise d' for each dataset. 
%statsmat = collate_byDataset(unitVSoutputs,unitoutputs,'classifier_tau',10);

% Summary of all the data.
statsmat_summary = summaryDatasets(statsmat);


%% -------------- SELECTION OF DATA ---------------

% Selection of data. Includes all levels, and modulation depth. 
% ~1000 datasets (each one varying only in modulation frequency).
% We remove any low-CFs, noise carriers, very short stimuli,
% & any datasets with very few AM frequencies. 
selection_criteria = { ...
    'fn','amToneDur>=100', ...
    'fn','carrFreq>3000', ...
    'fn','strcmp(amType,''AM'')', ...
    'fn','cf>3000', ...
    'fn','length(usedamfreqs)>5'...
    };

% This is all the data included in the paper.
statsmat_selected = select_Datasets(statsmat, selection_criteria{:});
% Summary of the units and stimulus conditions. 
statsmat_summary = summaryDatasets(statsmat_selected);

% STACK - reorganises data to a long single column format for all fields. 
% Remember: "dprime" is whatever your classifier statistic is. 
stats_stacked = stack_Datasets(statsmat_selected,'dprime');

% Also do this separately for different modulation depths. 
% Note: 'dprime' option here means we use the modulation frequencies for which
% the classifier analysis was ran.  
stats_stacked_mod1 = stack_Datasets(select_Datasets(statsmat_selected,'fn','depthMod==1'),'dprime');
stats_stacked_mod2 = stack_Datasets(select_Datasets(statsmat_selected,'fn','depthMod==2'),'dprime');
% Similar selection but with a modulation depth <1. Need both formats here
statsmat_modless1 = select_Datasets(statsmat_selected, ...
    'fn','depthMod<1');
summaryDatasets(statsmat_modless1)
stats_stacked_modless1 = stack_Datasets(select_Datasets(statsmat_selected,'fn','depthMod<1'),'dprime')

% Select out one dataset per unit for some of the analyses. 
% Looks for the closest one to 30 dB SPL, which is the lowest for most. 
% Modulation depth must be 1.
statsmat_oneperunit = select_Datasets(statsmat_selected, ...
    'fn','depthMod==1', ...
    'oneperunit','findnearest(30,[modLevel{:}])');
stats_oneperunit_stacked = stack_Datasets(statsmat_oneperunit,'dprime')

% Select the hit rate and F1 versions of the data. 
statsmat_selected_HR = select_Datasets(statsmat_HR, selection_criteria{:});
statsmat_selected_F1 = select_Datasets(statsmat_F1, selection_criteria{:});

% Some core criteria for some statistics & displays. 
coreoptions = { 'minn',4, ...
                'mean','sd','median'};
            
% Make this for plotting.
dPrime_reDp1 = makePopnTable(stats_stacked,'rationalisedType', ...
    'usedamfreqs_dprime_rationalised','dprimeReDp1', ...
    coreoptions{:},'fn','Rayleigh>13.8','roworder',[1 2 4 5 6 3]);     

            
%% ------------------------------------------------------------------------
%       Data for MTF comparisons - VS vs. Classifier statistic 

% This extracts BMF statistics for MTFs - so only keeps values at the BMF. 
% One for each dataset, or neuron.
% Needs to be done seprately for c' and VS meaurements. 
% N.B. It is run separately for 100 & 200% modulation.

% BMF-VS
amselection_criteria_BMF = ...
    {'amfn','usedamfreqs_reBMF==0', ...
     'amfn','allamfreqs_reBMF==0', ...
     'amfn','usedamfreqs_dprime_reBMF==0', ...
     'fn','VS_MTFCornerGoodness>=0.75'};
stats_selected_BMF_stacked_mod1 = stack_Datasets(select_Datasets(statsmat_selected, amselection_criteria_BMF{:},'fn','depthMod==1'),'usedamfreqs');  
stats_selected_BMF_stacked_mod2 = stack_Datasets(select_Datasets(statsmat_selected, amselection_criteria_BMF{:},'fn','depthMod==2'),'usedamfreqs');  
stats_oneperunit_BMF_stacked = stack_Datasets(select_Datasets(statsmat_oneperunit, amselection_criteria_BMF{:}),'usedamfreqs');

% BMF-d' - i.e. the chosen classifier statistic (c'). 
amselection_criteria_BMFdp = ...
    {'amfn','usedamfreqs_reBMFdp==0', ...
     'amfn','allamfreqs_reBMFdp==0', ...
     'amfn','usedamfreqs_dprime_reBMFdp==0', ...
     'fn','dprime_MTFCornerGoodness>=0.75'};
stats_selected_BMFdp_stacked_mod1 = stack_Datasets(select_Datasets(statsmat_selected, amselection_criteria_BMFdp{:},'fn','depthMod==1' ),'dprime'); 
stats_selected_BMFdp_stacked_mod2 = stack_Datasets(select_Datasets(statsmat_selected, amselection_criteria_BMFdp{:},'fn','depthMod==2' ),'dprime'); 
stats_oneperunit_BMFdp_stacked = stack_Datasets(select_Datasets(statsmat_oneperunit, amselection_criteria_BMFdp{:}),'usedamfreqs');
            

%% ------------------------------------------------------------------------
%                       MTF statistical analysis

% N.B. In the paper some of these are applied to the oneperunit dataset and
% some to all units. 

% Statements about VS in paper:

% Over all unit types - oneperunit
fprintf('Mean BMF-VS (all units): %g  s.d: %g\n', mean(stats_oneperunit_BMF_stacked.BMF), std(stats_oneperunit_BMF_stacked.BMF));
fprintf('%% of BMF-VS 200-800Hz: %g\n',100*sum(stats_oneperunit_BMF_stacked.BMF<=800 & stats_oneperunit_BMF_stacked.BMF>=200)/length(stats_oneperunit_BMF_stacked.BMF));
fprintf('%% of BMF-VS <200Hz: %g\n',100*sum(stats_oneperunit_BMF_stacked.BMF<200 )/length(stats_oneperunit_BMF_stacked.BMF));
%Mean BMF-VS (all units): 288.649  s.d: 188.041
% % of BMF-VS 200-800Hz: 61.0811
% % of BMF-VS <200Hz: 36.7568

% The range of sound levels in the "oneperunit" selection
fprintf('Sound levels for "oneperunit" selection: %g +/-%g\n', mean([statsmat_oneperunit.modLevel{:}]),std([statsmat_oneperunit.modLevel{:}]));
fprintf('\tRange %g - %g dB SPL\n', min([statsmat_oneperunit.modLevel{:}]),max([statsmat_oneperunit.modLevel{:}]));
fprintf('\t <=30dB SPL: %g %%\n', 100*sum([statsmat_oneperunit.modLevel{:}]<=30)/length([statsmat_oneperunit.modLevel{:}]) );
fprintf('\t 31-40dB SPL: %g %%\n', 100*sum([statsmat_oneperunit.modLevel{:}]>30 & [statsmat_oneperunit.modLevel{:}]<41)/length([statsmat_oneperunit.modLevel{:}]));
fprintf('\t 41-50dB SPL: %g %%\n', 100*sum([statsmat_oneperunit.modLevel{:}]>40 & [statsmat_oneperunit.modLevel{:}]<51)/length([statsmat_oneperunit.modLevel{:}]));
fprintf('\t >50dB SPL: %g %%\n', 100*sum([statsmat_oneperunit.modLevel{:}]>50 )/length([statsmat_oneperunit.modLevel{:}]));
figure;
hist([statsmat_oneperunit.modLevel{:}],[10:5:80]);
xlabel('Sound level (dB SPL)');ylabel('Dataset count')

% Peak VS vs. type (one per unit):
kruskalwallis(stats_selected_BMF_stacked_mod1.VS,stats_selected_BMF_stacked_mod1.rationalisedType)
% p=0. On>ChS~ChT>PLN>PBU>PL. 

% Peak VS vs. level. - decreases
kruskalwallis(stats_selected_BMF_stacked_mod1.VS,stats_selected_BMF_stacked_mod1.modLevel)

% BMF-VS.
kruskalwallis(stats_selected_BMF_stacked_mod1.BMF,stats_selected_BMF_stacked_mod1.rationalisedType)
%p=1e-4, Not a vast difference in median values. The difference is in the variances.

% Statements about classifier metric in paper (may differ slightly if classifier is re-run):

% Over all unit types
fprintf('Mean BMF-c'' (all units): %g  s.d: %g\n', mean(stats_oneperunit_BMFdp_stacked.BMFdp), std(stats_oneperunit_BMFdp_stacked.BMFdp));
fprintf('%% of BMF-c'' <200Hz: %g\n',100*sum( stats_oneperunit_BMFdp_stacked.BMFdp<200)/length(stats_oneperunit_BMFdp_stacked.BMFdp));
fprintf('%% of BMF-c'' = 50Hz: %g\n',100*sum( stats_oneperunit_BMFdp_stacked.BMFdp==50)/length(stats_oneperunit_BMFdp_stacked.BMFdp));
% Mean BMF-c' (all units): 147.701  s.d: 264.729
% % of BMF-c' <200Hz: 84.4828
% % of BMF-c' = 50Hz: 63.7931

% Peak d' vs. type:
kruskalwallis(stats_selected_BMFdp_stacked_mod1.dprimes,stats_selected_BMFdp_stacked_mod1.rationalisedType)
% p=0. ChS~ChT>On>PLN>PBU>PL. 

% d' vs. sound level.
kruskalwallis(stats_selected_BMFdp_stacked_mod1.dprimes,stats_selected_BMFdp_stacked_mod1.modLvl_rationalised)
% p=0, decreasing with level.
% Source       SS        df       MS      Chi-sq   Prob>Chi-sq
% ------------------------------------------------------------
% Groups   1.53749e+06     2   768743.9   77.32    1.62646e-17
% Error    8.14699e+06   485    16797.9                       
% Total    9.68448e+06   487                                                                 


%% ------------------------------------------------------------------------
%      Statistical modelling - predicting d' at different frequencies

% Here we fit a single statistical model across all the datasets in order
% to explore what other statistics and variables predict the classification
% performance. 

% Models of softmax_Z - takes Z (i.e. c') as the variable to predict. 
% [linearsoftmaxmodels dPrime_reBMFdp dPrime_reBMFdp_reBestDP stats_stacked_4models]= ...
%      statsModels_general(statsmat_selected,coreoptions,'Z','linear');
%save('datafiles\linearsoftmaxmodels','linearsoftmaxmodels'); % Takes 40 mins or so. 
load datafiles\linearsoftmaxmodels;   % Or just load it. 
nlm = linearsoftmaxmodels.nlm;
dPrime_reBMFdp = linearsoftmaxmodels.dPrime_reBMFdp;
dPrime_reBMFdp_reBestDP = linearsoftmaxmodels.dPrime_reBMFdp_reBestDP;
stats_stacked_4models = linearsoftmaxmodels.stats_stacked_4models;

% Models of hitrate.
% [HRmodels] = statsModels_general(statsmat_selected_HR,coreoptions,'HR','linear');
% save('datafiles\HRmodels','HRmodels'); 
% load datafiles\HRmodels; % Not provided.
% nlm = HRmodels.nlm;
% dPrime_reBMFdp = HRmodels.dPrime_reBMFdp;
% dPrime_reBMFdp_reBestDP = HRmodels.dPrime_reBMFdp_reBestDP;
% stats_stacked_4models = HRmodels.stats_stacked_4models;

% Models of F1
% [F1models] = statsModels_general(statsmat_selected_F1,coreoptions,'F1','linear');
% save('datafiles\F1models','F1models'); 
% load datafiles\F1models; % Not provided.
% dPrime_reBMFdp = F1models.dPrime_reBMFdp;
% dPrime_reBMFdp_reBestDP = F1models.dPrime_reBMFdp_reBestDP;
% stats_stacked_4models = F1models.stats_stacked_4models;
% nlm = F1models.nlm;
% 

% The % of data at 30,50 or 70 dB SPL (for a comment in the paper).
% NOTE: These were based on the entire dataset, which is an error!  
tmplvls = [statsmat.modLevel{:}];
100*sum(tmplvls==30 | tmplvls==50 | tmplvls==70)/length(tmplvls)

tmpdpths = [statsmat.depthMod{:}];
100*sum(tmpdpths==.5 | tmpdpths==1 | tmpdpths==2)/length(tmpdpths)



%% ------------------ Paper ready figures 4-7 -----------------

% Figures 1-3
% See ExampleUnitFigures

xlsxfilename = 'datavalues.xlsx';   % This is the file that the required data values are written to.

% Figure 4 Collected examples of MTFs.
% To produce this figure for the alternative classifier statistics:
%scorename ='Classifier hit rate';  scorelim = 1.1; filesuffix =  '_hitrate';
%scorename ='Classifier F1';  scorelim = 1.1; filesuffix =  '_F1';
scorename = 'Classifier peformance (c'')'; scorelim = 8.5; filesuffix =  '_Zsoftmax';
DPvsVSegs_paperfig_v3(statsmat_oneperunit,scorelim,scorename);
%print('-r600','-dtiff',['figures\Figure4' filesuffix '.tif']);
%saveas(gcf,['figures\Figure4' filesuffix],'fig')

% Figure 5 showing how the transfer function is similar across all neurons
%measure = 'hit rate'; maxscore = 1.05;  filesuffix =  '_hitrate'; mmin = -.02; mmax = 1.02; pmin = -0.02; pmax = 1.02;
%measure = 'F1'; maxscore = 1.05;  filesuffix =  '_F1'; mmin = -.02; mmax = 1.02; pmin = -0.02; pmax = 1.02;
measure = 'metric c'''; maxscore = 8.3;  filesuffix =  '_Zsoftmax'; mmin = -2; mmax = 8.1; pmin = -0.15; pmax = 10.5; xlsxsheet = 'Fig5B';
statsModels_paperfig; 
%saveas(gcf,['figures\Figure5' filesuffix ],'fig')
%print('-dtiff','-r600',['figures\Figure5' filesuffix '.tif']);

% Figure 6 showing low-level data - classifier statistic and VS on one figure. 
% measure = 'F1'; scorelims = [0 1.05];  filesuffix =  '_F1';xlsxnames = {'datavalues_FigS10.xlsx', 'FigS10B'};
%measure = 'hitrate'; scorelims = [0 1.05];  filesuffix =  '_hitrate'; xlsxnames = {'datavalues_FigS9.xlsx', 'FigS9B'};
measure = 'metric c'''; scorelims = [-0.5 7.8];  filesuffix =  '_Zsoftmax'; xlsxnames = {xlsxfilename, 'Fig6'};
%measure = 'metric c'' (tau = 10ms)'; scorelims = [-0.5 7.8];  filesuffix =  '_Zsoftmax'; xlsxnames = {'datavalues_FigS13.xlsx', 'FigS13'};
% N.B. For Figure 13, S13A is the same as 6A - added manually to datavalues_FigS13
VS_DP_lowLvl_paperfig_v3(stats_oneperunit_stacked, stats_oneperunit_BMFdp_stacked, ...
    stats_oneperunit_BMF_stacked, coreoptions,scorelims,['classifier ' measure]); %,xlsxnames);
%print('-dtiff','-r600',['figures\Figure6' filesuffix '.tif']);
%saveas(gcf,['figures\Figure6' filesuffix ])
% For the figure with fixed tau:
%print('-r600','-dtiff',['figures\Figure6' filesuffix '_tau10ms.tif']);
%saveas(gcf,'Figure6_tau10ms','fig')

% Figure 7 showing the chnages in modulation coding with sound level.
% This is for depth = 1 only. 
%measure =  'hit rate'; scorelims = [-0.05 1.1];  filesuffix =  '_hitrate';
%measure =  'F1'; scorelims = [-0.05 1.1];  filesuffix =  '_F1';
measure =  'metric c''';; scorelims = [-0.5 8];  filesuffix =  '_Zsoftmax';
DP_VS_3Lvls_paperfig_v3(stats_stacked_mod1, stats_selected_BMFdp_stacked_mod1, ...
    stats_selected_BMF_stacked_mod1, coreoptions,scorelims,measure);
set(gcf,'Name','depthMod = 1 only');
%print('-dtiff','-r600',['figures\Figure7' filesuffix '.tif']);
%saveas(gcf,['figures\Figure7' filesuffix],'fig')
% For the figure with fixed tau:
%saveas(gcf,'Figure7_tau10ms','fig')
%print('-dtiff','-r600',['figures\Figure7' filesuffix '_tau10ms.tif']);



%% -------------------- Population coding --------------------

% This code runs the population coding analysis and displays figure 8.
% ~2 hours on a typical PC. 
% Instead, load up from data files and paste in code for the figures. 
populationCodingMeasure = 'softmax_Z';
%populationCoding;
%save('datafiles\popnAnalByTypeLvl_Zsoftmax','growingPopnByTypeLvl','options','typeList','levels');

% Or for the alternative measures. Requires setting statistic using lines 61-62.  
%populationCodingMeasure = 'f1score'; populationCoding;
%save('datafiles\popnAnalByTypeLvl_F1','growingPopnByTypeLvl','options','typeList','levels');
%populationCodingMeasure = 'hitrate';populationCoding;
%save('datafiles\popnAnalByTypeLvl_hitrate','growingPopnByTypeLvl','options','typeList','levels');

% Or load instead.
load datafiles\popnAnalByTypeLvl_Zsoftmax; 
%load datafiles\popnAnalByTypeLvl_F1;
%load datafiles\popnAnalByTypeLvl_hitrate;

% Make the figure
%maxscore = 1.05; filesuffix = '_hitrate';
%maxscore = 1.05; filesuffix = '_F1'; scorename = '(F1)'
maxscore = 8.2; filesuffix = '_Zsoftmax'; scorename = '(C)'
populationCodingFigure(growingPopnByTypeLvl,options,populationCodingMeasure,maxscore,typeList,levels,scorename)
%print('-dtiff','-r600',['figures\Figure8' filesuffix '.tif']);
%saveas(gcf,['figures\Figure8_' filesuffix] ,'fig')

% One reviewer asked whether it was valid to group different neurons with
% different CFs for the population analysis. d' does not depend on CF,
% and the CFs in populations did not vary systematically. CFs are high too,
% so there is no concern about phase-locking to the carrier.
plotRoleOfCF;
% print -dtiff -r300 figures\R2_effectOfCF.tiff


%% ----------- Mode-locking and reliability analysis -----------

% This produces Figure 9.
% Need to run statsModels before this.
%maxscore = 1.05; scorename = 'Envelope classification (hit rate)';xlsxsheet = 'Fig9BCEF_hitrate';
maxscore = 8.1; scorename = 'Envelope classification (c'')'; xlsxsheet = 'Fig9BCEF';
%maxscore = 1.05; scorename = 'Envelope classification (F1)';xlsxsheet = 'Fig9BCEF_F1';

plotSpikeCodingAnalysis;
%print -dtiff -r600 figures\Figure9_hitrate.tif;
%print -dtiff -r600 figures\Figure9_Zsoftmax.tif;
%print -dtiff -r600 figures\Figure9_F1.tif;
%saveas(gcf,'figures\Figure9_Zsoftmax','fig')

% Writes all the model predictions and variables to an Excel file.
%makeFigure9AS6S8DataSheet('datavalues_Fig9A.xlsx',stats_stacked_4models,nlm,'9A');
%makeFigure9AS6S8DataSheet('datavalues_FigS9.xlsx',stats_stacked_4models,nlm,'9A');  %To write the Hit rate numbers
%makeFigure9AS6S8DataSheet('datavalues_FigS10.xlsx',stats_stacked_4models,nlm,'9A');  %To write the f1 numbers

%% ------------------ Supplementary materials ------------------

% N.B. This will not work without re-running SpikeStatsCalc (days) and 
% groupSpikeStatsByMTF (minutes) because it requires details of softmax fits 
% removed to save space in the github files. 
figure('Position',[100 100 700 800]);
chs_ind = 270; pl_ind = 582;
ahs = plotEgConfMatRow(unitoutputs(chs_ind),[0.8 0.15],'a.');
title(ahs(1),'Data'); title(ahs(2),'Softmax fits'); title(ahs(3),'Fitted parameters');
ahs = plotEgConfMatRow(unitoutputs(pl_ind),[0.55 0.15],'b.');
ahs = plotEgConfMatRow(unitoutputs(334),[0.3 0.15],'c.');
ahs = plotEgConfMatRow(unitoutputs(524),[0.05 0.15],'d.');
%print -dtiff -r300 figures\FigureS1_softmaxexamples.tif
%saveas(gcf,'figures\FigureS1_softmaxexamples','fig')

% CheckingSoftMaxFits gives us the quality of the fits 
% Requires softmax fits as above.

% Figure S2 - reproduction of Figure 7 with 200% modulation data.
DP_VS_3Lvls_paperfig_v3(stats_stacked_mod2,  stats_selected_BMFdp_stacked_mod2, ...
    stats_selected_BMF_stacked_mod2, coreoptions,[-0.5 10],'metric c''');
set(gcf,'Name','depthMod = 2 only');
%print('-r600','-dtiff','figures\FigureS2_200pcmod_Zsoftmax.tif');
%saveas(gcf,'figures\FigureS2_200pcmod_softmax','fig')

% Figure S3 - just plotting the MTF-d's for the few neurons we have.
plotShallowModulationDepths;
%print('-r600','-dtiff','figures\FigureS3_shallowmods.tif');
%saveas(gcf,'figures\FigureS3_shallowmods','fig')

% Figure S4 - just plotting the MTF-d's for the neurons with steps<50Hz.
maxscore = 8.3; plotSmallModFreqSteps;
% print('-r300','-dtiff','figures\FigureS4_smallamsteps.tif');
% saveas(gcf,'figures\FigureS4_smallamsteps','fig')

% Figure S5 - more details about ZISI and why the formal "mode-locking"
%             test fails. 
plotMLTestStatistics;
%print -dtiff -r300 figures\FigureS5_MLstats.tiff
%saveas(gcf,'figures\FigureS5_MLstats','fig')

% Neural fluctuation is a monotonic function of vector strength
% This simulation code supports the statement in the paper. 
controlModelAnalyses_vonMises;

% Figure S6 and S8 - showing a greater range of models and different SAC
%                    peaks
plotMoreStatisticalModels;
%print -dtiff -r300 figures\FigureS6_S8_moremodels.tiff;
%saveas(gcf,'figures\FigureS6_S8_moremodels','fig')
% Writes all the model predictions and variables to an Excel file.
makeFigure9AS6S8DataSheet('datavalues_FigS6_chk.xlsx',stats_stacked_4models,nlm,'S6');
makeFigure9AS6S8DataSheet('datavalues_FigS8.xlsx',stats_stacked_4models,nlm,'S8');

% Figure S7 - the proportions of SAC peaks at each frequency in each type. 
SACPeaks_numByFmodAndType(stats_oneperunit_stacked,'datavalues_S7.xlsx');
%print -dtiff -r300 figures\FigS7_SACpeaksByType.tif    

% Figures S9 and 10 are collages made from the hit rate and F1 figures in
% powerpoint. 

% Plot showing that(more or a less) taking just a statistic at ~150Hz is OK.
% This plot does not appear in the Supplmental, but some of the statistics
% are quoted.
plotStatisticSelfCorrelations;

% Figure S11 - Pairwise correlations between all the major statistics, including d'
PlotCorrelationTwixtStats;
%print -dtiff -r300 figures\FigureS11_statistic_correlations.tiff;
%saveas(gcf,'figures\FigureS11_statistic_correlations','fig')

% Figure S12 - Temporal resolution of the "best" classifier for each
%             dataset.
plotClassifierTaus;
% Also writes the values to an Excel spreadsheet.
%print -dtiff -r300 figures\FigureS12_bestclassifiertaus.tif
%saveas(gcf,'figures\FigureS12_bestclassifiertaus','fig')

% Figure S10 & 11 - The effect of fixing the classifier at tau = 10.
% Run this code uncommenting the "statsmat = " where tau is fixed at 10.
% The code Figure 6 and Figure 7 will produce the bottom row in each of S10 & S11.
% Crop, cut and paste to make Figure S10 & Figure S11.

