

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear all dimension variables.
clear   panelGap leftMargin bottomMargin  rastGap  bottomSideMargin rastX rastY rastAc ...
        phX  phY phAc  phGap isiX isiY  phY isiAc isiGap rastAx toneAx phAx isihAx ...
        raphLeftMargin thisUp graphSizeX graphSizeY graphYGap1  graphYGap2  imageSizeX  ...
        graphSizeX graphYGap3 dprimeX  dprimeY dprimeAc psthX psthY psthAc tmtfX tmtfY ...
        tmtfAc discUp discAc discX discY tmtfUp dprimeUp psthUp;

unit2runstr = num2str(unit2run);
load('rawdata\allUnitEgs',['unit' unit2runstr 'list'],['unit' unit2runstr 'spikes']);
thislist = eval(['unit' unit2runstr 'list']);
thesespikes = eval(['unit' unit2runstr 'spikes']);

% The unit number associated with these data.
this_eg_ind = find(unitids == unit2run); 

thisModLevel = unique(thislist(levCondInds,strcmp(EXPLOGLIST,'modLevel')));
if length(thisModLevel)>1
    error('Specified level is not unique');
end;

% Check that the model level for "this_eg_ind" is the one specified.
this_eg_modLevel = arrayfun(@(x) x.modLevel,  [unitoutputs(this_eg_ind).unitinfo]);

% this_eg_ind may not be unique. Most likely because there are more than
% one level.z
this_eg_ind = this_eg_ind(find(this_eg_modLevel == thisModLevel));

% A final specifier that Will make sure we get the right one - hand coded.
this_eg_ind = this_eg_ind(unitDataSetIndex);

% Check that the modulation frequency is 100%.
depthModChk = sum(thislist(:,strcmp(EXPLOGLIST,'depthMod'))~=1);

thislist = thislist(levCondInds,:);
thesespikes = thesespikes(levCondInds,:);


periods = 1000./freqConds;
analWind = [20 100];
% scatAxEnd = 15;
% scatAxLim = [0 scatAxEnd];

carrF = thislist(1,15);

% first collect everything that will need to be displayed
% if unit2run > 90000000
%     thisDir = bngFn('rhode_modelocking\fromrhode\cn91data');
% else
%     thisDir = bngFn('rhode_modelocking\fromrhode\cn88data');
% end
% thisDir = currentDir;
thisDir = 'rawdata'
expNo = unit2runstr(1:5);unitNo = unit2runstr(end-1:end);
eval(['load( ''' thisDir '\' num2str(expNo) 'unitStructs\Exp' num2str(expNo) 'U' num2str(unitNo) '.mat'' )']);
% eval(['load ' thisDir '\Exp' num2str(expNo) 'U' num2str(unitNo) '.mat']);
eval(['thisStruct = Exp' num2str(expNo) 'U' num2str(unitNo) 'struct;']);
eval(['thisModStruct = Exp' num2str(expNo) 'U' num2str(unitNo) 'ModStruct;']);
modStruct = thisModStruct(modStrChoice);
% figure
% for ii = 1:9
%     subplot(3,3,ii);
%     plot(modStruct.modPSTHtimes,squeeze(modStruct.modPSTHs(1,ii,:)));
% end

% PSTHs - normalised per second
psthLevs = thisStruct.psth.levels(psthChoice);
psthVals = thisStruct.psth.psthTenthms.values(psthChoice,:);
psthTimes = thisStruct.psth.psthTenthms.times;
psthTimes = psthTimes*1000; % in ms

% VS MTFs
sigSI = thislist(:,30);
allRay = thislist(:,32);

% rMTFs
scountVals = thislist(:,18);
scountVals(scountVals==0) = NaN;

modFreqs = thislist(:,13);
allModFreqs = modStruct.freqs';
levels = thislist(:,12);
allModLevs = unique(levels);
zScores = thislist(:,66);
overRef1 = thislist(:,70);
overRef1(isnan(overRef1)) = 0;
% CI vals
ciVals = thislist(:,60);
ciVals(ciVals==0) = NaN;

%% ------------- RERUN CLASSIFIER -------------

% This is really just to get the distance matrix. 
if runclassifier
    % Find the classifier specs.
    bestInd = unitoutputs(this_eg_ind).bestClassifier.index;
    bestClassifier = unitoutputs(this_eg_ind).wroutput(bestInd).classifier;
    bestClassifier.store_distances = 1;
    % Find which frequencies to run in the classifier. 
    pars.classifier_amfreq_inds = ismember(allModFreqs,unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions);
    wr_rerun = WR_Classifier([thesespikes{ pars.classifier_amfreq_inds}], ...
        unitoutputs(this_eg_ind).wroutput(bestInd).stimulusset,bestClassifier);
end;

%% ----------- CONTROL MODEL ------------

if runmodel

    % Derive the VS again - because we want the phase info too. 

    % Construct the input parameters for the control model. 
    pars.amfreqs = unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions;
    pars.nsweeps =  2*unique(unitoutputs(this_eg_ind).wroutput(bestInd).sweepspercondition);
    pars.dur_s = 0.08;
    pars.classifier_amfreq_inds = ismember(allModFreqs,unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions);
    pars.dcrate = scountVals(pars.classifier_amfreq_inds )/pars.dur_s;
    
    % Compute the angle from the spike times
    phi = cellfun(@(x,y)   2*pi*rem([x{:}],y)/y , ...
        thesespikes(pars.classifier_amfreq_inds),arrayfun(@(x) 1000/x,pars.amfreqs,'uni',false),'uni',false);
    cplx = cellfun(@(x)  cos( x ) + i*sin( x ),phi,'uni',false);
    ang = cellfun(@(x) angle(mean(x)), cplx);
    vs =  cellfun(@(x) abs(mean(x)), cplx);
     
    pars.phaselocking = [vs ang];
   
    % Checking. 
    % ind = 10
    % figure; hist(phi{ind},100); %[0:pi/32:2*pi])
    % hold on;
    % plot(ang(ind),15,'+'); ylim([0 20])

    % Find the classifier specs.
    bestInd = unitoutputs(this_eg_ind).bestClassifier.index;

    % Works well
    ctlModel = VSControlModel(pars,unitoutputs(this_eg_ind).wroutput(bestInd).classifier);

end;


%% ------------ the figure ---------------

colorgrey = [.4 .4 .4];
lightgrey = [.7 .7 .7];
colorctlmodel =  [0.8 0.5 0.5];

if newfig
    % drawing the figure
    g=figure;set(g,'PaperUnits','centimeters');
    if length(freqConds)==4
        xSize = 18;ySize = 14;
    else
         xSize = 18;ySize = 28;
    end;
    % centre on A4 paper
    xLeft = (21-xSize)/2;yTop = (30-ySize)/2;
    set(g,'paperposition',[xLeft yTop xSize ySize]);
    set(g,'units','centimeters')
    set(g,'position',[8 0 xSize ySize]);
end;

% ---------------------- sizes ---------------------------
panelGap = 1.3% .8; 

% Define some general values. 
if length(freqConds)==4
     leftMargin = 8; %7.5;
     bottomMargin = 4;
     rastGap = 0.6 %0.5;
     bottomSideMargin = 4;
     rastX = 4.5; %5.5; 
     rastY = 1.1; 
     rastAc = leftMargin; 
     phY = 1.7;   
     graphSizeY = 1.6; %1.8;
     graphYGap2 = .6;
     graphYGap3 = 1.1;% 1.2;
     graphLeftMargin =  panelGap +  0.6; % 1.1;
else   
    leftMargin = 12;
    bottomMargin = 7.5-rowOffset;
    rastGap =0.6;
    bottomSideMargin = 4-rowOffset;
    rastX =5.5; 
    rastY = 1.1; rastAc = leftMargin;  
    phY = 1.5;   
    graphSizeY = 1.35;
    graphYGap2 = .2;
    graphYGap3 = 1;
    graphLeftMargin =  panelGap +  0.2; % 1.1;
end;

% PH, ISIH and raster dimensions.
phX = 1.5; 
phAc = rastX+rastAc+panelGap; phGap = rastGap;
isiX = 1.7; isiY = phY; isiAc = phX+phAc+panelGap;isiGap = rastGap;
for jj=1:length(freqConds)
    thisUp = bottomMargin + (jj-1) * (phY+rastGap);
    rastAx(jj) = axes('position',[rastAc/xSize thisUp/ySize rastX/xSize rastY/ySize]);
    toneAx(jj) = axes('position',[rastAc/xSize (thisUp+rastY)/ySize rastX/xSize (phY-rastY-.05)/ySize]);
    if plotIHPH
        phAx(jj) = axes('position',[phAc/xSize thisUp/ySize phX/xSize phY/ySize]);hold on;
        if plotIH
            isihAx(jj) = axes('position',[isiAc/xSize thisUp/ySize isiX/xSize isiY/ySize]);hold on;
        end;
    end;
end


% rastAx = axes('position',[leftMargin/xSize bottomMargin/ySize 7/xSize 11/ySize]);hold on;
% phAx = axes('position',[8.5/xSize bottomMargin/ySize 2.5/xSize 11/ySize]);hold on;
% isihAx = axes('position',[11.5/xSize bottomMargin/ySize 2.5/xSize 11/ySize]);hold on;


graphSizeX = 3.8;
graphYGap1 = 1.5;
imageSizeX = graphSizeX;

dprimeX = graphSizeX;dprimeY = graphSizeY;dprimeAc = graphLeftMargin;
psthX = graphSizeX;psthY = graphSizeY;psthAc = graphLeftMargin;
tmtfX = graphSizeX;tmtfY = graphSizeY;tmtfAc = graphLeftMargin;

tmtfUp = bottomMargin;
dprimeUp = bottomMargin;

if length(freqConds)==4
    discAc = graphLeftMargin;
    discUp = bottomSideMargin+dprimeY+graphYGap2;
    discX = imageSizeX; discY = imageSizeX*.95;
else
    discUp = tmtfUp;
    discAc = graphLeftMargin+dprimeX+panelGap*1.5;
    discX = imageSizeX; discY = imageSizeX
end;

% tmtfUp = bottomSideMargin + (length(freqConds)==2)*discY;
% dprimeUp = bottomSideMargin + (length(freqConds)==2)*discY;
% psthUp = bottomSideMargin+dprimeY+(length(freqConds)==4)*(discY+graphYGap2) + graphYGap3;

psthUp = bottomMargin + dprimeY  +(length(freqConds)==4)*(discY + graphYGap2) + graphYGap3;



% PSTH AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
axes('position',[psthAc/xSize psthUp/ySize psthX/xSize psthY/ySize]);
bar(psthTimes,psthVals,'k');
maxY = round2(max(psthVals),100,'ceil');
set(gca,'xlim',[0 60],'ytick',[0 maxY],'yticklabel',[0 maxY],'ylim',[0 maxY],'FontSize',8);
%set(gca,'units','centimeters','box','off');
innPos = get(gca,'position');outPos = get(gca,'outerposition');furthLeft = outPos(1) - innPos(1);
xlabel('Time (ms)','fontweight','bold');yh=ylabel('Spikes/s','fontweight','bold');
box off;
%set(yh,'units','centimeters');
% if length(freqConds)==4
%     oldPos = get(yh,'position');
%     set(yh,'position',[furthLeft-.02 oldPos(2:3)])
% end;

% dprime vertical axis DDDDDDDDDDDDDDDDDDDDDDDDDDDD
axes('position',[tmtfAc/xSize tmtfUp/ySize tmtfX/xSize tmtfY/ySize]);
nonnanFreqs = allModFreqs(~isnan(sigSI));
topPanelXtick = 0:freqgap:nonnanFreqs(end);
set(gca,'ylim',[0 7],'ytick',0:2:6,'xtick',topPanelXtick, ...
    'xlim',[nonnanFreqs(1)-10 nonnanFreqs(end)+10],'FontSize',8);
ylabel('d prime','fontweight','bold','color','b');
xlabel('Modulation frequency (Hz)','fontweight','bold');

% tMTF DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
axes('position',[tmtfAc/xSize tmtfUp/ySize tmtfX/xSize tmtfY/ySize]);
plot(nonnanFreqs,sigSI(~isnan(sigSI)),'k','linewidth',2)
set(gca,'ylim',[0 1],'ytick',[0 .5 1],'xtick',topPanelXtick,'xlim',[nonnanFreqs(1)-10 nonnanFreqs(end)+10],...
    'xticklabel','','FontSize',8, ...
    'Yaxislocation','right');
ylabel('VS','fontweight','bold');


% dprime FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
% Plots d' on same axis/6.
hold on;
nonnanFreqs = allModFreqs(~isnan(sigSI));
topPanelXtick = 0:freqgap:nonnanFreqs(end);
bestInd = unitoutputs(this_eg_ind).bestClassifier.index;
% plot(bestwroutputs(this_eg_ind).stimulusconditions, ...
%     bestwroutputs(this_eg_ind).dprime/6,'b','linewidth',2)
% dprime for control model 
if ~isempty(who('ctlModel'))
    plot(ctlModel.wr.stimulusconditions, ctlModel.wr.dprime/7, ...
        'color',colorctlmodel,'linewidth',2);
end;

plot(unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions, ...
    unitoutputs(this_eg_ind).wroutput(bestInd).dprime/7,'b','linewidth',2)



isiGoto = 20;
rastStart = 20;
rastGoto = 80;
isibins = 0:.5:60;
shuffisibins = 0:.5:60;
phbins = 0:1/32:(1-1/32);
sampF = 1/100000;stimulustime =  [rastStart/1000:sampF:(rastGoto/1000-sampF)];

for jj=1:length(freqConds)
    thisCond = find(modFreqs==freqConds(jj));
    % spiketimes
    theseSpikeTimes = thesespikes{thisCond};
    tmptrain = sTimesCellTo2D([theseSpikeTimes{:}],length(theseSpikeTimes));
    statsout=getModelockCriteria_scholes(tmptrain, freqConds(jj), analWind, 0);
    tmptrain(tmptrain<analWind(1) | tmptrain>analWind(2)) = NaN;
    
    if plotIHPH && plotIH
        % ISIH %%%%% FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        scatterVals = diff(tmptrain); allScatVals = scatterVals(~isnan(scatterVals));
        isiHist = hist(allScatVals,isibins);
        ISIphsurVals = hist(statsout.PHsur.ISI,shuffisibins);

        allisi2plot = isiHist;
        % SHOW 90% OF ISIS
        %     cumISI = cumsum(allisi2plot);
        %     cumBins = isibins(cumISI<=cumISI(end)*.95);

        % ...OR SHOW 3.5 PERIODS
        cumBins = periods(jj)*3.5;

        % ...OR MANUAL
        if freqConds(jj)<=100
            cumBins = periods(jj)*1.5;
        elseif freqConds(jj)<=200
            cumBins = periods(jj)*2.5;
        elseif freqConds(jj)<=300
            cumBins = periods(jj)*4.5;
        elseif freqConds(jj)<=400
            cumBins = periods(jj)*6.5;
        else
            cumBins = periods(jj)*7.5;
        end

        for kk = 1:floor(cumBins/periods(jj))
            axes(isihAx(jj));
            line([kk*periods(jj) kk*periods(jj)],[0 2000],'color',lightgrey,'linestyle','--');
        end
        bar(isihAx(jj),isibins,allisi2plot,'facecolor','k');
        %     plot(isihAx(jj),isibins,allisi2plot,'k','linewidth',1);
        plot(isihAx(jj),shuffisibins,ISIphsurVals,'r','linewidth',1);

        set(isihAx(jj),'yticklabel','','xlim',[0 cumBins(end)],'ytick',[0 isiMax],'fontsize',8);
        maxISIs(jj) = ceil(max([allisi2plot ISIphsurVals]));
    end;
    
    % RASTER %%%%%%%%%% DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    spikeMat = cell(1,theseNSweeps);
    for kk=1:theseNSweeps
        % Get this spike train.
        tmptrain2 = thesespikes{thisCond}{kk};tmptrain2 = tmptrain2(tmptrain2<110);
        spikeMat{kk} = tmptrain2;
    end
    
    newspikeMat = sTimesCellTo2D(spikeMat,theseNSweeps);
    
    thisTrialArr = 1:theseNSweeps;
    triallen = 110;fs=1000;trialgap=1;
    newspikeMat(isnan(newspikeMat))=0;
    trials=repmat(thisTrialArr,size(newspikeMat,1),1);
    reltimes=mod(newspikeMat,triallen);
    reltimes(~reltimes)=NaN;
    numspikes=prod(size(newspikeMat));
    xx=ones(3*numspikes,1)*nan;
    yy=ones(3*numspikes,1)*nan;
    trialNums = repmat((1:theseNSweeps),3*size(newspikeMat,1),1);
    trialNums = trialNums(:);
    
    yy(1:3:3*numspikes)=(trials(:)-1)*trialgap;
    yy(2:3:3*numspikes)=yy(1:3:3*numspikes)+1;
    
    %scale the time axis to ms
    xx(1:3:3*numspikes)=reltimes(:)*1000/fs;
    xx(2:3:3*numspikes)=reltimes(:)*1000/fs;
    xlim=[rastStart,rastGoto];
    
    plot(rastAx(jj),xx, yy, 'k', 'linewidth',1.5);
    axis (rastAx(jj),[xlim,0,(theseNSweeps*trialgap)]);
    set(rastAx(jj),'xtick',0:20:100,'ytick',[],'xticklabel','','yticklabel','');
    
    stimuluswaveform = cos(stimulustime*2*pi*carrF) ...
    .*(1- 1*cos(2*pi*freqConds(jj)*stimulustime));

    plot(toneAx(jj),stimulustime,stimuluswaveform,'color',colorgrey);
    axis(toneAx(jj),[rastStart/1000 rastGoto/1000 -2 2],'off')
    
    if plotIHPH    
        % PH %%%%%%% EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        stInWin1D = tmptrain(~isnan(tmptrain));
        ph = (stInWin1D - floor(stInWin1D/periods(jj))*periods(jj))/periods(jj);
        phvals = hist(ph,phbins);
        phIsurVals = hist(statsout.Isur.ph,phbins);
        % centre PHs using last min value
        newshift = rem(find(phvals==min(phvals),1,'last'),length(phvals));
        newshift = 0;
        ph2plot=[phvals(newshift+1:end) phvals(1:newshift)];
        %     max(ph2plot)
        
        if ~isempty(who('ctlModel'))
           modelCond = find(ctlModel.wr.stimulusconditions==freqConds(jj));
           plot(phAx(jj),phbins,ctlModel.ph{modelCond}/2,'color',colorctlmodel, ...
               'linewidth',2);
        end;

        
        bar(phAx(jj),phbins,ph2plot,'facecolor','k');
        %plot(phAx(jj),phbins,phIsurVals,'r','linewidth',1);
        set(phAx(jj),'xtick',0:.5:1,'xlim',[0 1],'xticklabel','','yticklabel','','ylim',[0 phMax],'ytick',[0 phMax]);
        axes(phAx(jj));
        text(0.05,.9*phMax,['VS = ' num2str(round2(sigSI(thisCond),.01))],'fontweight','bold','FontSize',8);
    
    end;
end
set(rastAx(1),'xticklabel',0:20:100,'fontsize',8)
xlabel(rastAx(1),'Time (ms)','fontsize',8,'fontweight','bold');

if plotIHPH
    set(phAx(1),'xticklabel',0:.5:1,'yticklabel',[0 phMax],'fontsize',8)
    xlabel(phAx(1),'Phase','fontsize',8,'fontweight','bold');
    if plotIH
        set(isihAx(:),'ylim',[0 max(maxISIs)],'ytick',[0 max(maxISIs)]);
        set(isihAx(1),'yticklabel',[0 max(maxISIs)],'fontsize',8)
        xlabel(isihAx(1),'Interval (ms)','fontsize',8);
    end;
end;

% labels
% psthX = graphSizeX/xSize;psthY = graphSizeY/ySize;psthAc = graphLeftMargin/xSize;psthUp = graphPosY/ySize;
% tmtfX = graphSizeX/xSize;tmtfY = graphSizeY/ySize;tmtfAc = (graphLeftMargin+graphSizeX+graphGap)/xSize;tmtfUp = graphPosY/ySize;
% sppX = graphSizeX/xSize;sppY = graphSizeY/ySize;sppAc = (graphLeftMargin+graphSizeX*2+graphGap*2)/xSize;sppUp = graphPosY/ySize;
%
% rastX = 2.4/xSize;rastY = 1.5/ySize;rastAc = leftMargin/xSize;rastUp = 8.3/ySize;rastGap = 0.5;
% phX = 2.4/xSize;phY = 1.5/ySize;phAc = leftMargin/xSize;phUp = 5.9/ySize;phGap = 0.5;
% sacX = 2.4/xSize;sacY = 1.4/ySize;sacAc = leftMargin/xSize;sacUp = 1/ySize;sacGap = 0.5;
% scatXY = 2.4;scatAc = leftMargin/xSize;scatUp = 2.7/ySize;scatGap = 0.5;scatBoxGap = 0.05;

axes('position',[0 0 1 1],'visible','off');
% LETTERS
if length(freqConds)==4    
    abc_x = .02;
    text(.005,.95,unitType,'fontweight','bold','FontSize',16);
    text(abc_x,.9,'a','fontweight','bold','FontSize',14);
    text(abc_x,.7,'b','fontweight','bold','FontSize',14);
    text(abc_x,.42,'c','fontweight','bold','FontSize',14);
    text(.39,.9,'d','fontweight','bold','FontSize',14);
    text(.72,.9,'e','fontweight','bold','FontSize',14);
    if plotIH
        text(.87,.9,'f','fontweight','bold','FontSize',14);
    end;
    if plotConfMat && plotIH
        text(.7,.47,'f','fontweight','bold','FontSize',14)
    end;
else
    text(.005,(psthUp+psthY+1.3)/ySize,unitType,'fontweight','bold','FontSize',16)
    text(.01,(psthUp+psthY+0.8)/ySize,'a','fontweight','bold','FontSize',14)
    text(.01,(tmtfUp+tmtfY+0.4)/ySize,'b','fontweight','bold','FontSize',14)
    text(.39,(psthUp+psthY+0.8)/ySize,'c','fontweight','bold','FontSize',14)
    text(.63,(psthUp+psthY+0.8)/ySize,'d','fontweight','bold','FontSize',14)
%     text(.005,.95-rowOffset/ySize,unitType,'fontweight','bold','FontSize',16)
%     text(.035,.9-rowOffset/ySize,'a','fontweight','bold','FontSize',14)
%     text(.035,.7-rowOffset/ySize,'b','fontweight','bold','FontSize',14)
%     text(.39,.9-rowOffset/ySize,'c','fontweight','bold','FontSize',14)
%     text(.76,.9-rowOffset/ySize,'d','fontweight','bold','FontSize',14)
end;   
    
% FREQUENCIES LABEL
phX = 3; phY = 1.7; phAc = rastX+rastAc+panelGap; phGap = rastGap;
for ii =1:length(freqConds)
    if length(freqConds)==4
       thisUp = (bottomMargin + (ii-1) * (phY+rastGap) + phY/2) / ySize;
       fx = 0.38; %0.35
    else
        thisUp = (0.75 + bottomMargin + (ii-1) * (phY+rastGap) + phY/2) / ySize;
        fx = 0.67
    end;
    text(fx,thisUp,[num2str(freqConds(ii)) 'Hz'],'FontSize',10,'fontweight','bold');
end

if plotConfMat 
    % discrim performance EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    axes('position',[discAc/xSize discUp/ySize discX/xSize discY/ySize]);
    %imagesc(bestwroutputs(this_eg_ind).stimulusconditions,bestwroutputs(this_eg_ind).stimulusconditions,bestwroutputs(this_eg_ind).confusionmatrix);
    imagesc(unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions, ...
        unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions, ...
        unitoutputs(this_eg_ind).wroutput(bestInd).confusionmatrix);

    if length(freqConds)==4
        set(gca,'FontSize',8,'xtick',topPanelXtick,'ytick',topPanelXtick,'xticklabel','');
    else
        set(gca,'FontSize',8,'xtick',topPanelXtick,'ytick',topPanelXtick);
        xlabel('Modulation frequency (Hz)','fontweight','bold');
    end;        
    
    %set(gca,'units','centimeters');
    innPos = get(gca,'position');outPos = get(gca,'outerposition');furthLeft = outPos(1) - innPos(1);
    yh=ylabel('Classifier choice (Hz)','fontweight','bold');
%     if length(freqConds)==4
%          set(yh,'units','centimeters');oldPos = get(yh,'position');
%          set(yh,'position',[furthLeft oldPos(2:3)]);
%     end;
end;

