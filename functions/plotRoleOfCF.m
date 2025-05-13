% --- Role of CF in discrimination individually and in population ----

xSize = 20; ySize = 12;
figh=figure;set(figh,'PaperUnits','centimeters');
set(figh,'paperposition',[1 1 xSize ySize],'units','centimeters','position',[2 4 xSize ySize],'DefaultAxesFontSize',8);%,'DefaultAxesColorOrder',colorOrder

typeColorMap;
subplot(1,2,1);
for ti = 1:6
    hold on;
    inds = strmatch( typeorder{ti},stats_oneperunit_stacked.rationalisedType);
    semilogx(stats_oneperunit_stacked.cf(inds),stats_oneperunit_stacked.dprimes(inds), ...
        '.','color',typecolormap(ti,:));
end;
legend(typeorder);
xlabel('CF (Hz)'); ylabel(['Classifier ' measure]);
title('Discrimination as a function of CF (all f_m_o_d)');
ylim([-1 10]);

% Population analysis.
%load datafiles\popnAnalByTypeLvl;

subplot(1,2,2);
allcfs = [];
for ti = 1:6
    cflist = cellfun(@(x) x.unitinfo.carrFreq,growingPopnByTypeLvl(1,ti).best1st.eachUnit,'Uni',false);
    semilogy([cflist{:}],'color',typecolormap(ti,:));
    allcfs = [allcfs [cflist{:}]];
    hold on;
end;
legend(typeorder);
ylabel('CF (Hz)'); xlabel('Order added to cluster');
ylim([2.4e3 35e3]);
set(gca,'ytick',[2.5e3, 5e3, 10e3, 20e3, 30e3]); box off;
title('The CF of neuron added to a cluster, best Z 1st');

fprintf('CF statistics:\n');
fprintf('\tMean: %g kHz, s.d.: %g kHz\n',mean(allcfs)/1e3,std(allcfs/1e3));
fprintf('\tAbove 5kHz: %g %%\n',100*sum(allcfs>5e3)/length(allcfs));

% print -dtiff -r300 effectOfCF.tiff