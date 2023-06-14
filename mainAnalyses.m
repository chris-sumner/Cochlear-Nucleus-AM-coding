%% ----------------------------------------------------------------------

% Need to run WolgRonacClassification and VSCalc first. 


%% -----------------------------------------------------------------------
clear

% Load up and collate all the data into a format where each measure is an
% array of cells for each dataset.

% Path required.
addpath functions;

load('datafiles\VS','unitVSoutputs');           % Load VS analysis.
load('datafiles\WR_results','unitoutputs');     % load WR classifier performance.

% This integrates all the information into a single set of cell arrays.
statsmat = collate_byDataset(unitVSoutputs,unitoutputs,'dontcheck')
% Summary of all the data.
statsmat_summary = summaryDatasets(statsmat);


%% -------------- SELECTION OF DATA ---------------

% Selection of data. Includes all levels, and modulation depth. 
% ~1000 datasets.
selection_criteria = { ...
    'fn','amToneDur>=100', ...
    'fn','carrFreq>3000', ...
    'fn','strcmp(amType,''AM'')', ...
    'fn','cf>3000', ...
    'fn','length(usedamfreqs)>5'...
    };

% This is all the data included in the paper currently.
statsmat_selected = select_Datasets(statsmat, selection_criteria{:});

% Summary of the units and stimulus conditions. 
statsmat_summary = summaryDatasets(statsmat_selected);
% Allows you to check you have the stuff you want. 
% Can also specify which frequency units to use. 
summaryDatasets(statsmat_selected,'usedamfreqs_rationalised')

% Select out one dataset per unit for some of the analyses. 
% Looks for the closest one to 30 dB SPL. 
% Modulation depth must be 1.
statsmat_oneperunit = select_Datasets(statsmat_selected, ...
    'fn','depthMod==1', ...
    'oneperunit','findnearest(30,[modLevel{:}])');
 
% Similar seletion of one per unit but with a modulation depth of 200%.
statsmat_oneperunit_mod2 = select_Datasets(statsmat_selected, ...
    'fn','depthMod==2', ...
    'oneperunit','findnearest(30,[modLevel{:}])');

% STACK - different ways to stack - which frequencies to use. 
stats_stacked = stack_Datasets(statsmat_selected,'dprime');
stats_oneperunit_stacked = stack_Datasets(statsmat_oneperunit,'dprime');
stats_oneperunit_mod2_stacked = stack_Datasets(statsmat_oneperunit_mod2,'dprime');

% Some core criteria for some statistics. 
coreoptions = { 'minn',4, ...
                'mean','sd','median'};

            
%% ------------------------------------------------------------------------
%               Quantitative MTF comparisons - VS vs. d''

% This exttracts BMF statistics MTFs - VS and d. One for each neuron.

% -------------------- BMF-VS vs. BMF-d' --------------------

% TO PERFORM ANALYSIS ON A DIFFERENT MODULATION DEPTH - CHANGE depthMod IN ALL BOTH PLACES.

% Select out BMF datasets.  Makes constructing tables easy.     
% Also:
%   - only take units where the characterisation of the MTF was resasonable. 

% BMF-VS
amselection_criteria_BMF = ...
    {'amfn','usedamfreqs_reBMF==0', ...
     'amfn','allamfreqs_reBMF==0', ...
     'amfn','usedamfreqs_dprime_reBMF==0', ...
     'fn','VS_MTFCornerGoodness>=0.75', ...
     'fn','depthMod==1'};
statsmat_selected_BMF = select_Datasets(statsmat, selection_criteria{:}, amselection_criteria_BMF{:} );
stats_selected_BMF_stacked = stack_Datasets(statsmat_selected_BMF,'usedamfreqs');
statsmat_oneperunit_BMF = select_Datasets(statsmat_oneperunit, selection_criteria{:}, amselection_criteria_BMF{:} );
stats_oneperunit_BMF_stacked = stack_Datasets(statsmat_oneperunit_BMF,'usedamfreqs');

% BMF-d'
amselection_criteria_BMFdp = ...
    {'amfn','usedamfreqs_reBMFdp==0', ...
     'amfn','allamfreqs_reBMFdp==0', ...
     'amfn','usedamfreqs_dprime_reBMFdp==0', ...
     'fn','dprime_MTFCornerGoodness>=0.75', ...
     'fn','modDepth==1'}
statsmat_selected_BMFdp = select_Datasets(statsmat, selection_criteria{:}, amselection_criteria_BMFdp{:} );
stats_selected_BMFdp_stacked = stack_Datasets(statsmat_selected_BMFdp,'dprime');
statsmat_oneperunit_BMFdp = select_Datasets(statsmat_oneperunit, selection_criteria{:}, amselection_criteria_BMFdp{:} );
stats_oneperunit_BMFdp_stacked = stack_Datasets(statsmat_oneperunit_BMFdp,'usedamfreqs');

            

%% ------------------------------------------------------------------------
%                       MTF statistical analysis

% N.B. In the paper some of these are applied to the oneperunit dataset and
% some to all units. 

% Statements about VS in paper:

% Over all unit types - oneperunit.
fprintf('Mean BMF-VS (all units): %g  s.d: %g\n', mean(stats_oneperunit_BMF_stacked.BMF), std(stats_oneperunit_BMF_stacked.BMF));
fprintf('%% of BMF-VS 200-800Hz: %g\n',100*sum(stats_oneperunit_BMF_stacked.BMF<=800 & stats_oneperunit_BMF_stacked.BMF>=200)/length(stats_oneperunit_BMF_stacked.BMF));
fprintf('%% of BMF-VS <200Hz: %g\n',100*sum(stats_oneperunit_BMF_stacked.BMF<200 )/length(stats_oneperunit_BMF_stacked.BMF));
%Mean BMF-VS (all units): 288.649  s.d: 188.041
% % of BMF-VS 200-800Hz: 61.0811
% % of BMF-VS <200Hz: 36.7568

% Range of BMFs in a given type - all selected units.
thisType = 'ChS';
ChI = find(strcmp(stats_selected_BMF_stacked.rationalisedType,thisType) );
fprintf('Mean BMF-VS (%s): %g\n',thisType, mean(stats_selected_BMF_stacked.BMF(ChI)));
fprintf('%% of BMF-VS below 200Hz: %g\n',100*sum(stats_selected_BMF_stacked.BMF(ChI)<200)/length(ChI));
fprintf('%% of BMF-VS below 100Hz: %g\n',100*sum(stats_selected_BMF_stacked.BMF(ChI)<100)/length(ChI));
fprintf('%% of BMF-VS 200-500Hz: %g\n',100*sum(stats_selected_BMF_stacked.BMF(ChI)<500 & stats_selected_BMF_stacked.BMF(ChI)>200)/length(ChI));
% Mean BMF-VS (ChS): 302.381
% % of BMF-VS below 200Hz: 22.2222
% % of BMF-VS below 100Hz: 5.55556
% % of BMF-VS 200-500Hz: 71.4286

% Peak VS vs. type (one per unit):
kruskalwallis(stats_selected_BMF_stacked.VS,stats_selected_BMF_stacked.rationalisedType)
% p=0. On>ChS~ChT>PLN>PBU>PL. 

% Peak VS vs. level. - decreases
kruskalwallis(stats_selected_BMF_stacked.VS,stats_selected_BMF_stacked.modLevel)

% Corner frequency bandwidth (-3dB) vs. type.
kruskalwallis(stats_selected_BMF_stacked.VS_MTFCornerBW,stats_selected_BMF_stacked.rationalisedType)
% p=0. On>PLN~PL>ChT>PBU~ChS

% Upper corner frequency vs. type.
kruskalwallis(stats_selected_BMF_stacked.VS_HiCornerHz,stats_selected_BMF_stacked.rationalisedType)
%p=0, On>PL~PLN>ChT>PBU~Chs.

% BMF-VS.
kruskalwallis(stats_selected_BMF_stacked.BMF,stats_selected_BMF_stacked.rationalisedType)
%p=1e-4, Not a vast difference in median values. The difference is in the variances.

% Statements about d' in paper:

% Over all unit types
fprintf('Mean BMF-d'' (all units): %g  s.d: %g\n', mean(stats_oneperunit_BMFdp_stacked.BMFdp), std(stats_oneperunit_BMFdp_stacked.BMFdp));
fprintf('%% of BMF-dp <200Hz: %g\n',100*sum( stats_oneperunit_BMFdp_stacked.BMFdp<200)/length(stats_oneperunit_BMFdp_stacked.BMFdp));
fprintf('%% of BMF-dp = 50Hz: %g\n',100*sum( stats_oneperunit_BMFdp_stacked.BMFdp==50)/length(stats_oneperunit_BMFdp_stacked.BMFdp));
%Mean BMF-d' (all units): 138.92  s.d: 278.166
% % of BMF-dp <200Hz: 86.9318
% % of BMF-dp = 50Hz: 70.4545

% Peak d' vs. type:
kruskalwallis(stats_selected_BMFdp_stacked.dprimes,stats_selected_BMFdp_stacked.rationalisedType)
% p=0. ChS~ChT>On>PLN>PBU>PL. 

% BMF-d'.
kruskalwallis(stats_selected_BMFdp_stacked.BMFdp,stats_selected_BMFdp_stacked.rationalisedType)
%p=0.0, Howwever, few above 200Hz

% d' corner frequency vs. type.
kruskalwallis(stats_selected_BMFdp_stacked.dprime_MTFCornerBW,stats_selected_BMFdp_stacked.rationalisedType)
% p=0. 

% d' vs. sound level.
kruskalwallis(stats_selected_BMFdp_stacked.dprimes,stats_selected_BMFdp_stacked.modLevel)
% p=0, decreasing with level.

%% ------------------------------------------------------------------------
%      Statistical modelling - predicting d' at different frequencies

% And here we fit a single model across all the datasets - accounting for 80% 
% of the variance of >13000 d' measurements with 6 parameters and a linear
% model.
statsModels;

% The % of data not at 30,50 or 70 dB SPL (for a comment in the paper).
tmplvls = [statsmat.modLevel{:}];
100*sum(tmplvls==30 | tmplvls==50 | tmplvls==70)/length(tmplvls)

tmpdpths = [statsmat.depthMod{:}];
100*sum(tmpdpths==.5 | tmpdpths==1 | tmpdpths==2)/length(tmpdpths)



%% ------------------ Paper ready figures 4-7 -----------------

% Figures 1-3
% See ExampleUnitFigures

% Figure 4 Collected examples of MTFs.
DPvsVSegs_paperfig_v3(statsmat_oneperunit)
%print('-r600','-dtiff','figures\Figure4.tif');

% Figure 5 showing how the transfer function is similar across all neurons
statsModels_paperfig; 
%print -dtiff -r600 figures\Figure5.tif;

% Figure 6 showing low-level data - presenting d' and compating with VS on
% one figure. 
VS_DP_lowLvl_paperfig_v2(stats_oneperunit_stacked, ...
    stats_oneperunit_BMFdp_stacked,  ....
    stats_oneperunit_BMF_stacked, coreoptions);
%print('-r600','-dtiff','figures\Figure6.tif');

% Figure 7 showing the chnages in modulation coding with sound level. 
% This selects a modulation depth of 1 within the figure. 
DP_VS_3Lvls_paperfig_v2(stats_stacked,  stats_selected_BMFdp_stacked, stats_selected_BMF_stacked, coreoptions)
%print('-r600','-dtiff','figures\Figure7.tif');


%% -------------------- Population coding --------------------

% This code runs the population coding analysis and displays figure 8.
populationCoding;





