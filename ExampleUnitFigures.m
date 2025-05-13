%% ------------------------------------------------------------------------
%              Template for making the example figures
%  ------------------------------------------------------------------------

clear;

% Load up and collate all the data into a format where each measure is an
% array of cells for each dataset.

% Path's required.
addpath functions;
addpath('functions\Wohlgemuth Ronacher\Matlab');
addpath('functions\Wohlgemuth Ronacher\C');

% Processed data files required. 
%load('datafiles\VS','unitVSoutputs');
load('datafiles\SpikeStatsSets','unitVSoutputs');
load datafiles\WR_results_corrected_extra2;
load('rawdata\unitList','EXPLOGLIST');       

% Extract unit outputs and a unitid (indentifier with which to find the
% unit).
unitinfo = [unitoutputs(:).unitinfo];
unitids = [unitinfo(:).unitid];

 
%% PL ***************************************************************
unitType = 'PL';
unit2run = 88340053;
unitDataSetIndex = 1;
levCondInds = 27:52; 
psthChoice = 2;
modStrChoice=1;
isiMax = 12;
phMax=35;
theseNSweeps = 25;
freqConds = 50:100:350;
freqgap = 600;
scatAxEnd = [22 22 22] ;

plotConfMat = 1;
plotIHPH = 1;
plotIH = 0;
newfig = true;

runmodel =  true; % Only needs to be run once
runclassifier = true;

metric_scalefactor = 8;
%UnitFigTemplate_main_Feb22
UnitFigTemplate_main_Mar25

%print -dtiff -r600 figures\Figure3_PL_exFig.tif
%saveas(gcf,'figures\Figure3_PL_exFig','fig')

  
%% ChS ***************************************************************
unitType = 'ChS';
unit2run = 88299021;
levCondInds = 79:117;
unitDataSetIndex = 3;
psthChoice = 2; modStrChoice=3;
isiMax = 45;
phMax=60;
theseNSweeps = 10;
freqConds = [75 125 175 400 ];
freqgap = 200;
plotConfMat = 1;
plotIHPH = 1; 
plotIH = 0;
newfig = true;

runmodel = true; 
runclassifier = true;

metric_scalefactor = 8;
UnitFigTemplate_main_Mar25;

%print -dtiff -r600 figures\Figure2_ChS_exFig.tif
%saveas(gcf,'figures\Figure2_ChS_exFig','fig')
% print -dtiff -r600 ChS_exFig_tau2ms.tif

% Method figure for the paper.
% Must be run AFTER the previous function. 
WRMethodFig_Feb22;
%print -dtiff -r600 figures\Figure1_WRmethodFig.tif
%saveas(gcf,'figures\Figure1_WRmethodFig','fig')


 
 
 
 

