% -------------------- More statistical models figure ------------------------

% For the supplemental material this was split into two figures via
% cropping in Word. 

xSize = 28; ySize = 18; xLeft = (21-xSize)/2;yTop = (30-xSize)/2;
figh=figure;set(figh,'PaperUnits','centimeters');
set(figh,'paperposition',[xLeft yTop xSize ySize],'units','centimeters','position',[2 4 xSize ySize],'DefaultAxesFontSize',8);%,'DefaultAxesColorOrder',colorOrder

% a. An expanded 
ah(1) = axes('position',[.05 .64 .45 .31]);hold on;
% Need to run the models before this. 
RegressionAllModelsExactPanel;
text(min(xlim),max(ylim)*1.05,'a. Statistical models of envelope classification in all data','fontweight','bold');

ah(2) = axes('position',[.57 .64 .4 .31]);hold on;
RegressionSACModelsPanel;
text(min(xlim),max(ylim)*1.05,'b. Adding combinations of SAC peak-trough values','fontweight','bold');


% print -dtiff -r300 FigureS5_S7_moremodels.tiff;

