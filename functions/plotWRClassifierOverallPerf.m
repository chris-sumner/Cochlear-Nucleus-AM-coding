function  [filename] = plotWRClassifierOverallPerf(unitoutput,thisdatastartind,options)

wroutput = unitoutput.wroutput;
unitinfo = unitoutput.unitinfo;

% Plot the results for each tau
figure; set(gcf,'position',[100 100 800 800]);
for i=1:length(unitinfo.classifierTauList)
    ah(1)=  subplot(4,4,4+2*(i-1)+1);
    ah(2)=  subplot(4,4,4+2*(i-1)+2);
    plotOutput_WRClassifier(wroutput(i),ah);
    subplot(ah(1)); ylabel(['Choice tau:' num2str(unitinfo.classifierTauList(i))]);
end;

% Plot the VS
subplot(4,3,1); set(gca,'fontsize',8);
plot(unitinfo.usedamfreqs,unitinfo.VS(unitinfo.conditionswithspikes));
xlabel('f_m_o_d (Hz)'); ylabel('VS'); xlim([unitinfo.usedamfreqs(1) unitinfo.usedamfreqs(end)]);
% Also plot CV here.
text(max(xlim)*.7,max(ylim)*.8,['CV:' num2str(rnd2ndp(unitinfo.CV,2))]);

% Plot the Z and spikes per cycle.
subplot(4,3,2); set(gca,'fontsize',8);
plot(unitinfo.usedamfreqs,unitinfo.Z(unitinfo.conditionswithspikes));
hold on; plot(unitinfo.allamfreqs(unitinfo.MLPoissTest==1),unitinfo.Z(unitinfo.MLPoissTest==1),'.r','markersize',10);

xlim([unitinfo.usedamfreqs(1) unitinfo.usedamfreqs(end)]);
plot(unitinfo.usedamfreqs,unitinfo.spikesperperiod(unitinfo.conditionswithspikes),'r');
xlabel('f_m_o_d (Hz)'); ylabel('Z (blue) and S/P (red)');

specstr = sprintf('#%d, UNIT:%d (%.3gkHz,%gdBSPL), %s, %ddB SPL, fC:%gHz, modD:%g', ...
    thisdatastartind,  unitinfo.unitid, round(unitinfo.cf)/1e3, round(unitinfo.cf_thr), ...
    unitinfo.TypeName, unitinfo.modLevel,unitinfo.carrFreq, unitinfo.depthMod);
title(specstr,'fontsize',12,'fontweight','bold');

% Plot the dprimes for all taus.
subplot(4,3,3);
dprimes = reshape([wroutput(:).dprime],length(wroutput(1).dprime),length(unitinfo.classifierTauList));
plot(unitinfo.usedamfreqs,dprimes);
xlim([unitinfo.usedamfreqs(1) unitinfo.usedamfreqs(end)]);
lh = legend(num2str(unitinfo.classifierTauList'));   set(lh,'fontsize',6);
xlabel('f_m_o_d (Hz)'); ylabel('d''');

if nargin>2 && strcmp(options,'save');
    filename = ['N' num2str(thisdatastartind) '_U' num2str(unitinfo.unitid) '.tiff' ];
    print('-dtiff','-r200',['WR_figsAllFreqs\' filename]);
else
    filename = [];
end;