function [periodlag, SACvalues] = plotSACwithStats(unit,unitVSoutputs,statsmat)
%UNTITLED2 Summary of this function goes here

% Extract all the information required for display from unitVSoutputs and
% statsmat
unit.ModPeriod  = 1e3/unit.ModF;
unit.Inds = find([statsmat.unitid{:}] == unit.unitID); % All the datasets for this unit. 
unit.Ind = unit.Inds(unit.dataSetNum);                 % The dataset we want.  

unit.amFreqInd = find([unitVSoutputs(unit.Ind).SACpeaks.amfreq] == unit.ModF);
unit.lagLimits = [-1.1 1.1].*unit.ModPeriod;
unit.lagAxisDt =  diff(unitVSoutputs(unit.Ind).SACpeaks(unit.amFreqInd).lagaxis(1:2));
unit.lagaxis = unitVSoutputs(unit.Ind).SACpeaks(unit.amFreqInd).lagaxis;
unit.SAC = unitVSoutputs(unit.Ind).SACpeaks(unit.amFreqInd).smoothsac;
unit.SAC = unit.SAC(unit.lagaxis>=min(unit.lagLimits*1.2) & unit.lagaxis<=max(unit.lagLimits*1.2));
unit.lagaxis = unit.lagaxis(unit.lagaxis>=min(unit.lagLimits*1.2) & unit.lagaxis<=max(unit.lagLimits*1.2));

% In some units the peaks are wrong (?) in SACpeaks - corrected in unitVSoutputs.
unit.peakLags = unitVSoutputs(unit.Ind).SACpeaks(unit.amFreqInd).setofsigpeaklags_p0001;

unit.trofLags =  unitVSoutputs(unit.Ind).SACpeaks(unit.amFreqInd).lagaxis( ...
    unitVSoutputs(unit.Ind).SACpeaks(unit.amFreqInd).setofsigtrofs_p0001 );

% Select out the peaks and troughs of interest.
% These are the ones we use for analysis.  
unit.peakLags = unit.peakLags(unit.peakLags>=min(unit.lagLimits) & unit.peakLags<=max(unit.lagLimits));
unit.trofLags = unit.trofLags(unit.trofLags>=min(unit.lagLimits) & unit.trofLags<=max(unit.lagLimits));

% Plot the SAC
plot(unit.lagaxis/unit.ModPeriod, unit.SAC,unit.lineStyle);
hold on; xlim([-1.3 1.3]);box off;

periodlags = unit.lagaxis/unit.ModPeriod;
SACvalues = unit.SAC;

% Plot the Peaks.
[~, unit.peakInds] = ismember( unit.peakLags, unit.lagaxis);
%plot((unit.peakLags+ unit.lagAxisDt)/unit.ModPeriod,unit.SAC(unit.peakInds+1),unit.peakStyle);
plot((unit.peakLags)/unit.ModPeriod,unit.SAC(unit.peakInds),unit.peakStyle);

% Plot vertical lines at the periods and zero. 
line([1 1],[0 6],'linestyle','--','color','k');
line(-[1 1],[0 6],'linestyle','--','color','k');
line([0 0],[0 6],'linestyle','--','color','k');
xlabel('Autocorrelation lag (modulation periods)','fontweight','bold'); 
ylabel('Coincidence count','fontweight','bold');

% Add in trough and salience measures. 
[~, unit.trofInds] = ismember( unit.trofLags, unit.lagaxis);
plot(unit.trofLags/unit.ModPeriod,unit.SAC(unit.trofInds),unit.trofStyle);

% % And salience measures - finding them again for double goodness. 
unit.peak0_lagind = find(unit.lagaxis == statsmat.SACpeak0lag_p0001{unit.Ind}(unit.amFreqInd))
plotSalienceBar(unit.saliencePositions(1),unit.SAC,unit.peak0_lagind, ...
    statsmat.SACpeak0sal_p0001{unit.Ind}(unit.amFreqInd),unit.salienceLineColour,0);
 
unit.peak1_lagind = find(unit.lagaxis == statsmat.SACpeak1lag_p0001{unit.Ind}(unit.amFreqInd))
plotSalienceBar(unit.saliencePositions(2),unit.SAC,unit.peak1_lagind, ...
    statsmat.SACpeak1sal_p0001{unit.Ind}(unit.amFreqInd),unit.salienceLineColour,1);

unit.peak2_lagind = find(unit.lagaxis == statsmat.SACpeak2lag_p0001{unit.Ind}(unit.amFreqInd))
plotSalienceBar(unit.saliencePositions(3),unit.SAC,unit.peak2_lagind, ...
    statsmat.SACpeak2sal_p0001{unit.Ind}(unit.amFreqInd),unit.salienceLineColour,2);

unit.peak3_lagind = find(unit.lagaxis == statsmat.SACpeak3lag_p0001{unit.Ind}(unit.amFreqInd))
plotSalienceBar(unit.saliencePositions(4),unit.SAC,unit.peak3_lagind, ...
    statsmat.SACpeak3sal_p0001{unit.Ind}(unit.amFreqInd),unit.salienceLineColour,3);

unit.peak4_lagind = find(unit.lagaxis == statsmat.SACpeak4lag_p0001{unit.Ind}(unit.amFreqInd))
plotSalienceBar(unit.saliencePositions(5),unit.SAC,unit.peak4_lagind, ...
    statsmat.SACpeak4sal_p0001{unit.Ind}(unit.amFreqInd),unit.salienceLineColour,4);

unit.peak5_lagind = find(unit.lagaxis == statsmat.SACpeak5lag_p0001{unit.Ind}(unit.amFreqInd))
plotSalienceBar(unit.saliencePositions(6),unit.SAC,unit.peak5_lagind, ...
    statsmat.SACpeak5sal_p0001{unit.Ind}(unit.amFreqInd),unit.salienceLineColour,5);
text(1.3,5.8,'Peak-to-trough sizes');

% Number the peaks in the main plot
% n.b. We take them from the extracted statistics in stats mat so we 
% gauranteed the correct order. 
poslags = [ statsmat.SACpeak0lag_p0001{unit.Ind}(unit.amFreqInd) ...
            statsmat.SACpeak1lag_p0001{unit.Ind}(unit.amFreqInd) ...
            statsmat.SACpeak2lag_p0001{unit.Ind}(unit.amFreqInd) ...
            statsmat.SACpeak3lag_p0001{unit.Ind}(unit.amFreqInd) ...
            statsmat.SACpeak4lag_p0001{unit.Ind}(unit.amFreqInd) ...
            statsmat.SACpeak5lag_p0001{unit.Ind}(unit.amFreqInd) ...
            ];
%poslags = unit.peakLags(unit.peakLags>=0.05*min(unit.lagLimits) & unit.peakLags<=max(unit.lagLimits))
[~, posinds] = ismember( poslags, unit.lagaxis);
peaknums = [0:length(poslags)-1];
if ~isempty(peaknums)
    arrayfun(@(x,y,n) text(x,y,num2str(n),'color',unit.salienceLineColour), ...
            (poslags + unit.lagAxisDt)/unit.ModPeriod , ...
            unit.SAC(posinds+1)+0.25, peaknums);
end;

% Add information about the unit. 
txt = sprintf('%d: %s (CF:%d), @%d Hz', unit.unitID,  statsmat.type{unit.Ind}, ...
    statsmat.cf{unit.Ind},unit.ModF);
text(-1,5.5,txt);

xlim([-1.5 max(unit.saliencePositions)+.1]);
set(gca,'xtick',[-1:0.5:1]);
text(min(xlim),max(ylim)*1.1,'d. Examples of shuffled autocorrelograms','fontweight','bold');



function plotSalienceBar(xpos,SAC,peakInd,salience,linecolour,linenum)

if ~isempty(peakInd)
    line([xpos xpos], ...
        [SAC(peakInd+1) SAC(peakInd+1) - salience], ...
        'linewidth',2,'color',linecolour);
    
    text(xpos-0.05,SAC(peakInd+1)+.2,num2str(linenum),'color',linecolour);
    
end;
