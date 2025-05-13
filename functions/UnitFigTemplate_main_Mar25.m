

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear all dimension variables.
clear   panelGap leftMargin bottomMargin  rastGap  bottomSideMargin rastX rastY rastAc ...
        phX  phY phAc  phGap isiX isiY  phY isiAc isiGap rastAx toneAx phAx isihAx ...
        raphLeftMargin thisUp graphSizeX graphSizeY graphYGap1  graphYGap2  imageSizeX  ...
        graphSizeX graphYGap3 dprimeX  dprimeY dprimeAc psthX psthY psthAc tmtfX tmtfY ...
        tmtfAc discUp discAc discX discY tmtfUp dprimeUp psthUp;

unit2runstr = num2str(unit2run);

% ------------- Load up the spike times -----------

% All the spike times are stored in a special structure for that neuron - for convenience.  
load('rawdata\allUnitEgs',['unit' unit2runstr 'list'],['unit' unit2runstr 'spikes']);
thislist = eval(['unit' unit2runstr 'list']);         % This is an old set of stats and conditions for a very early analysis for the  .
thesespikes = eval(['unit' unit2runstr 'spikes']);    % This is the set of all spike times for the neuron - all conditions 

% ----- Look up the neuron in unit outputs and do some checks ------

% Find the unit number associated with these data within unitoutputs
unitoutputs_notempty = find(arrayfun(@(x) ~isempty(x.unitinfo.unitid),  unitoutputs));
this_eg_ind = unitoutputs_notempty(find(arrayfun(@(x) (x.unitinfo.unitid == unit2run),  unitoutputs(unitoutputs_notempty))));

% THIS WAS NOT DOING IT RIGHT.
%this_eg_ind = find(unitids == unit2run); 

% Check the unit number is correct.
this_eg_in_unitoutputs_chk = arrayfun(@(x) x.unitid,  [unitoutputs(this_eg_ind).unitinfo]);
fprintf('Checking unit ID is correct... should be %d.\nFound: ', unit2run);
fprintf('%d, ', this_eg_in_unitoutputs_chk); fprintf('\n');

% Retrieve the modulation level stored in the old table of stats.
thisModLevel = unique(thislist(levCondInds,strcmp(EXPLOGLIST,'modLevel')));
if length(thisModLevel)>1
    error('Specified level is not unique');
end;
fprintf('Checking sound level. Should be %d.\nFound: ',thisModLevel);

% Look up the sound levels in unitoutputs  
this_eg_modLevel = arrayfun(@(x) x.modLevel,  [unitoutputs(this_eg_ind).unitinfo]);
fprintf('%d, ', this_eg_modLevel); fprintf('\n');

% this_eg_ind may not be unique. Most likely because there are more than
% one level.
this_eg_ind = this_eg_ind(find(this_eg_modLevel == thisModLevel));
fprintf('Position of data of %d dB SPL in unitoutputs:',thisModLevel);
fprintf('%d,', this_eg_ind); fprintf('\n');

% A final specifier in case there are multiple records at the same
% modulation level.
this_eg_ind = this_eg_ind(unitDataSetIndex);
if unitDataSetIndex!=1
    fprintf('Taking option %d of these.\n',unitDataSetIndex);   
end;
    
% Check that the modulation frequency is 100%.
depthModChk = sum(thislist(:,strcmp(EXPLOGLIST,'depthMod'))~=1);
if depthModChk == 1
    error('Modulation depth is not 100%');
else
    fprintf('Modulation depth is confirmed to be 100%%.\n');
end;

% With all the checks out of the way we can select the right spike times 
% and the stats list specific to them. 
thislist = thislist(levCondInds,:);
thesespikes = thesespikes(levCondInds,:);

% -------------- Load up the type PSTHs and prepare  --------------------

thisDir = 'rawdata';
expNo = unit2runstr(1:5);unitNo = unit2runstr(end-1:end);
eval(['load( ''' thisDir '\' num2str(expNo) 'unitStructs\Exp' num2str(expNo) 'U' num2str(unitNo) '.mat'' )']);
eval(['thisStruct = Exp' num2str(expNo) 'U' num2str(unitNo) 'struct;']);
eval(['thisModStruct = Exp' num2str(expNo) 'U' num2str(unitNo) 'ModStruct;']);
modStruct = thisModStruct(modStrChoice);

% PSTHs - normalised per second
psthLevs = thisStruct.psth.levels(psthChoice);
psthVals = thisStruct.psth.psthTenthms.values(psthChoice,:);
psthTimes = thisStruct.psth.psthTenthms.times;
psthTimes = psthTimes*1000; % in ms

% ------------- Initialise some important variables ---------------

periods = 1000./freqConds;
analWind = [20 100];
carrF = thislist(1,15);

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
    
    wsmooth = 0.0001;  % This is a very weak constraint. 
    wnoise = 0.0001;   % This amount of noise works well to stabilise the fit.
    niter = 10;        % Noise is added to initial parameters.  
    Zmax = 8;         % Final limit on Z becuase it makes no sense to have very high values.  
    plims = [0 .99];
    
    % Find the classifier specs.
    bestInd = unitoutputs(this_eg_ind).bestClassifier.index;
    bestClassifier = unitoutputs(this_eg_ind).wroutput(bestInd).classifier;
    %bestClassifier = unitoutputs(this_eg_ind).wroutput(1).classifier;
    %bestClassifier = unitoutputs(this_eg_ind).wroutput(2).classifier;
    bestClassifier.store_distances = 1;
    % Find which frequencies to run in the classifier (run them all). 
    modFreqs = allModFreqs(allModFreqs<=2300 & ~isnan(scountVals));
    pars.classifier_amfreq_inds = 1:length(modFreqs); % ismember(allModFreqs,unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions);
    stimulusset = ones(length(thesespikes{1}),1)*modFreqs';
    wr_rerun = WR_Classifier([thesespikes{ pars.classifier_amfreq_inds}], ...
        stimulusset,bestClassifier);
    [softmax_Z softmax_B softmax_err softmax_fit] = fitSoftMax2ConfMat(wr_rerun.confusionmatrix, ...
        wsmooth,wnoise,niter,plims,Zmax);
    wr_rerun.softmax_Z = softmax_Z;
    wr_rerun.softmax_B = softmax_B;
    wr_rerun.softmax_err = softmax_err;
    wr_rerun.softmax_fit = softmax_fit;
%wr_rerun = WR_Classifier([thesespikes{ pars.classifier_amfreq_inds}], ...
    %    unitoutputs(this_eg_ind).wroutput(bestInd).stimulusset,bestClassifier);
end;

%% ----------- CONTROL MODEL ------------

if runmodel

    % Derive the VS again - because we want the phase info too. 

    % Construct the input parameters for the control model. 
    pars.dur_s = 0.08;
    pars.dcrate = scountVals(pars.classifier_amfreq_inds )/pars.dur_s;
    pars.amfreqs = wr_rerun.stimulusconditions(~isnan(scountVals(pars.classifier_amfreq_inds ))) ; % wr_rerun.stimulusconditions(~isnan(scountVals)) ; %  unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions;
    pars.nsweeps = 2*unique(wr_rerun.sweepspercondition); 
    pars.classifier_amfreq_inds = ismember(allModFreqs,pars.amfreqs); 
    
    % Compute the angle from the spike times
    phi = cellfun(@(x,y)   2*pi*rem([x{:}],y)/y , ...
        thesespikes(pars.classifier_amfreq_inds),arrayfun(@(x) 1000/x,pars.amfreqs,'uni',false),'uni',false);
    cplx = cellfun(@(x)  cos( x ) + i*sin( x ),phi,'uni',false);
    ang = cellfun(@(x) angle(mean(x)), cplx);
    vs =  cellfun(@(x) abs(mean(x)), cplx);
     
    pars.phaselocking = [vs ang];
   
    % Find the classifier specs.
    bestInd = unitoutputs(this_eg_ind).bestClassifier.index;

    % Works well
    ctlModel = VSControlModel(pars,unitoutputs(this_eg_ind).wroutput(bestInd).classifier);

    % Softmax fit to model.
    [softmax_Z softmax_B softmax_err softmax_fit] = fitSoftMax2ConfMat(ctlModel.wr.confusionmatrix, ...
            wsmooth,wnoise,niter,plims,Zmax);
    ctlModel.wr.softmax_Z = softmax_Z;
    ctlModel.wr.softmax_B = softmax_B;
    ctlModel.wr.softmax_err = softmax_err;
    ctlModel.wr.softmax_fit = softmax_fit;

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
    discX = imageSizeX*1.38; discY = imageSizeX*.95;
else
    discUp = tmtfUp;
    discAc = graphLeftMargin+dprimeX+panelGap*1.5;
    discX = imageSizeX; discY = imageSizeX
end;
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

% dprime vertical axis DDDDDDDDDDDDDDDDDDDDDDDDDDDD
axes('position',[tmtfAc/xSize tmtfUp/ySize tmtfX/xSize tmtfY/ySize]);
nonnanFreqs = allModFreqs(~isnan(sigSI));
topPanelXtick = 0:freqgap:nonnanFreqs(end);
set(gca,'ylim',[-0.1*metric_scalefactor metric_scalefactor*1.1],'ytick',[0:4:8],'xtick',topPanelXtick, ...
    'xlim',[nonnanFreqs(1)-10 nonnanFreqs(end)+10],'FontSize',8);
ylabel('c'' metric','fontweight','bold','color','b');
xlabel('Modulation frequency (Hz)','fontweight','bold');

% tMTF DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
axes('position',[tmtfAc/xSize tmtfUp/ySize tmtfX/xSize tmtfY/ySize]);
plot(nonnanFreqs,sigSI(~isnan(sigSI)),'k','linewidth',2)
set(gca,'ylim',[-0.1 1.1],'ytick',[0 .5 1],'xtick',topPanelXtick,'xlim',[nonnanFreqs(1)-10 nonnanFreqs(end)+10],...
    'xticklabel','','FontSize',8, ...
    'Yaxislocation','right');
ylabel('VS','fontweight','bold');


% dprime FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
% Plots d' on same axis/6.
hold on;
nonnanFreqs = allModFreqs(~isnan(sigSI));
topPanelXtick = 0:freqgap:nonnanFreqs(end);
bestInd = unitoutputs(this_eg_ind).bestClassifier.index;
% dprime for control model 
if ~isempty(who('ctlModel'))
    plot(ctlModel.wr.stimulusconditions, ctlModel.wr.softmax_Z/metric_scalefactor, ...
        'color',colorctlmodel,'linewidth',2);
end;

% plot(unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions, ...
%     unitoutputs(this_eg_ind).wroutput(bestInd).dprime/7,'b','linewidth',2)
plot(wr_rerun.stimulusconditions, ...
    wr_rerun.softmax_Z/metric_scalefactor,'b','linewidth',2);

% % Also plot B.
% plot(wr_rerun.stimulusconditions, ...
%     wr_rerun.softmax_B/7,'b:','linewidth',2);

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
    statsout=getModelockCriteria(tmptrain, freqConds(jj), analWind, 0);
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

    % Convert to proportion of choices.
    nchoice = sum(wr_rerun.confusionmatrix,1);    
    pchoice = wr_rerun.confusionmatrix./(nchoice(1)*ones(length(wr_rerun.stimulusconditions)))
    imagesc(wr_rerun.stimulusconditions, ...
        wr_rerun.stimulusconditions, ...
        pchoice);

    if length(freqConds)==4
        set(gca,'FontSize',8,'xtick',topPanelXtick,'ytick',topPanelXtick,'xticklabel','');
    else
        set(gca,'FontSize',8,'xtick',topPanelXtick,'ytick',topPanelXtick);
        xlabel('Modulation frequency (Hz)','fontweight','bold');
    end;        
    
    ch = colorbar('TickLabels',{});  
    ch_pos = get(ch,'Position');
    text(max(pars.amfreqs)*1.075,max(pars.amfreqs)*1.05,num2str(rnd2ndp(min(pchoice(:)),1)));
    text(max(pars.amfreqs)*1.05,-30,num2str(rnd2ndp(max(pchoice(:)),1)));
    innPos = get(gca,'position');outPos = get(gca,'outerposition');furthLeft = outPos(1) - innPos(1);
    yh=ylabel('Classifier choice (Hz)','fontweight','bold');
end;

