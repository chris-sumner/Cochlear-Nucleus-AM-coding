%% --------------------------- Figure 9 ----------------------------------

% Most panels are at one fmod. 
stackinds = (stats_oneperunit_stacked.usedamfreqs== 150 | stats_oneperunit_stacked.usedamfreqs== 100) & ...
              cellfun(@(x) ~isempty(x),stats_oneperunit_stacked.rationalisedType);

% We use a simple sum of the 1st for peaks for illustration here. 
stats_oneperunit_stacked.SACpeaks1234sum_150ish_p0001 = ...
    stats_oneperunit_stacked.SACpeak0sal_150ish_p0001 + ...
    stats_oneperunit_stacked.SACpeak2sal_150ish_p0001 + ...
    stats_oneperunit_stacked.SACpeak3sal_150ish_p0001 + ...
    stats_oneperunit_stacked.SACpeak4sal_150ish_p0001;
                    
% ----------- Make the figure ----------

xSize = 28; ySize = 18; xLeft = (21-xSize)/2;yTop = (30-xSize)/2;
figh=figure;set(figh,'PaperUnits','centimeters');
set(figh,'paperposition',[xLeft yTop xSize ySize],'units','centimeters','position',[2 4 xSize ySize],'DefaultAxesFontSize',8);%,'DefaultAxesColorOrder',colorOrder
typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;

% a. Summary of correlations across measures
jj=1; ah(jj) = axes('position',[.05 .64 .35 .31]);hold on;
% Need to run the models before this. 
RegressionAllModelsPanel;
text(min(xlim),max(ylim)*1.05,'a. Statistical models of envelope classification in all data','fontweight','bold');

% c. Example SACs.
jj=4; ah(jj) = subplot('position',[.05 .07 .35 .3]); hold on; 
plotPanel_EgSACs;

% This was the draft correlations figure for the paper at one point. 
plotDPrimeCorrelations8(stats_oneperunit_stacked,stackinds,maxscore,scorename,xlsxsheet)

