%% --- How might choice of summary statisic affect models? ---

% - All statistics drop at modulation frequency (Z less than others), so we take the lowest. 
% - But which modulation frequencies should we include in the average? 
% The conclusion is it probably doesn't matter very much. 

% ----- Looking at how adjacent values are correlated -----

% We look at the correlations between adjacent values. 
% This plot goes a bit overboard. 
% Column 1 shows the correlations for all pairs below ~1kHz fmod (quoted in Supplemental). 
% Column 2 shows correlations between frequencies 1-9 and 3-11. etc.
% Correlations are strong for all thse values - though gradually dropping
% as the frequency different increases.
% The 2nd peak amplitude of the SAC is the worst but that is not surprising
% because so many neurons don't have one - or only at low frequencies. 

xSize = 20; ySize = 28; xLeft = (21-xSize)/2;yTop = (30-xSize)/2;
figh=figure;set(figh,'PaperUnits','centimeters');
set(figh,'paperposition',[xLeft yTop xSize ySize],'units','centimeters','position',[2 4 xSize ySize],'DefaultAxesFontSize',8);%,'DefaultAxesColorOrder',colorOrder

subplot(5,4,1); 
autoCorrStatsAndPlot(statsmat_oneperunit.Z,1,10,'Z');
subplot(5,4,2); 
autoCorrStatsAndPlot(statsmat_oneperunit.Z,1,9,'Z',1);
subplot(5,4,3); 
autoCorrStatsAndPlot(statsmat_oneperunit.Z,1,8,'Z',2);
subplot(5,4,4); 
autoCorrStatsAndPlot(statsmat_oneperunit.Z,1,7,'Z',3);

subplot(5,4,5); 
autoCorrStatsAndPlot(statsmat_oneperunit.envFluct_0p48ms,1,10,'NuFl');
subplot(5,4,6); 
autoCorrStatsAndPlot(statsmat_oneperunit.envFluct_0p48ms,1,9,'NuFl',1);
subplot(5,4,7); 
autoCorrStatsAndPlot(statsmat_oneperunit.envFluct_0p48ms,1,8,'NuFl',2);
subplot(5,4,8); 
autoCorrStatsAndPlot(statsmat_oneperunit.envFluct_0p48ms,1,7,'NuFl',3);
 
subplot(5,4,9); 
autoCorrStatsAndPlot(statsmat_oneperunit.reliability_R_0p24ms,1,10,'Reliability');
subplot(5,4,10); 
autoCorrStatsAndPlot(statsmat_oneperunit.reliability_R_0p24ms,1,9,'Reliability',1);
subplot(5,4,11); 
autoCorrStatsAndPlot(statsmat_oneperunit.reliability_R_0p24ms,1,8,'Reliability',2);
subplot(5,4,12); 
autoCorrStatsAndPlot(statsmat_oneperunit.reliability_R_0p24ms,1,7,'Reliability',3);
 
subplot(5,4,13); 
autoCorrStatsAndPlot(statsmat_oneperunit.SACpeak1lag_p0001,1,10,'SAC peak 1');
subplot(5,4,14); 
autoCorrStatsAndPlot(statsmat_oneperunit.SACpeak1lag_p0001,1,9,'SAC peak 1',2);
subplot(5,4,15); 
autoCorrStatsAndPlot(statsmat_oneperunit.SACpeak1lag_p0001,1,8,'SAC peak 1',3);
subplot(5,4,16); 
autoCorrStatsAndPlot(statsmat_oneperunit.SACpeak1lag_p0001,1,7,'SAC peaks 1',4);
 
subplot(5,4,17); 
autoCorrStatsAndPlot(statsmat_oneperunit.SACpeak2lag_p0001,1,10,'SAC peak 2');
subplot(5,4,18); 
autoCorrStatsAndPlot(statsmat_oneperunit.SACpeak2lag_p0001,1,9,'SAC peak 2',2);
subplot(5,4,19); 
autoCorrStatsAndPlot(statsmat_oneperunit.SACpeak2lag_p0001,1,8,'SAC peak 2',3);
subplot(5,4,20); 
autoCorrStatsAndPlot(statsmat_oneperunit.SACpeak2lag_p0001,1,7,'SAC peak 2',4);


