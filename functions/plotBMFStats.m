function plotBMFStats(stats_selected_BMF_stacked,stats_selected_BMFdp_stacked,coreoptions,options)


% Using the BMF stacked stats to make tables and histograms of the
% BMF for both statistics. 
rowandcol =  {'rationalisedType','usedamfreqs_reBMF'};
stats.BMFtable = makePopnTable(stats_selected_BMF_stacked, rowandcol{:} ,'BMF',coreoptions{:}, ...
    options{:},'hist',[50:100:2050],'fn','Rayleigh>13.8');
rowandcol =  {'rationalisedType','usedamfreqs_reBMFdp'};
stats.BMFdptable = makePopnTable(stats_selected_BMFdp_stacked, rowandcol{:} ,'BMFdp',coreoptions{:}, ...
    options{:},'hist',[50:100:2050]);

% This shows the radically different range of BMFs.
figure('position',[100 20 1000 1000],'paperposition',[.5 .5 20 20]); 

subplot(3,2,1);
plot([50:100:2050],squeeze(stats.BMFtable.hist)','o-');
xlabel('BMF-VS (Hz)'); ylabel('# units'); %legend(stats.BMFtable.rowValues);
ylim([0 26]); xlim([0 1500]); box off;
text(0, 28,'A. BMF-VS distribution','fontweight','bold');

subplot(3,2,2);
plot([50:100:2050],squeeze(stats.BMFdptable.hist)','o-');
xlabel('BMF-d'' (Hz)'); ylabel('# units'); lh = legend(stats.BMFdptable.rowValues);
ylim([0 26]);  xlim([0 1500]); box off;
set(lh,'box','off');
text(0, 28,'B. BMF-d'' distribution','fontweight','bold');

% % Box plots. 
% a1 =subplot(2,1,1);
% boxplot(stats_selected_BMF_stacked.BMF,stats_selected_BMF_stacked.rationalisedType);
% set(a1,'yscale','log','ytick',[25 100 500 2000]); ylabel('BMF-d'' (Hz)');
% ylim([20 4e3]); box off;  text(.75,5e3,'A. Distribution of best-VS','fontweight','bold');
% 
% a2 = subplot(2,1,2);
% boxplot(stats_selected_BMFdp_stacked.BMFdp,stats_selected_BMFdp_stacked.rationalisedType);
% set(a2,'yscale','log','ytick',[25 100 500 2000]); xlabel('Unit type'); ylabel('BMF-VS (Hz)');
% ylim([20 4e3]); box off; text(.75,5e3,'B. Distribution of best d''','fontweight','bold');

%  Box plots across level.
% This shows that BMF gradually shifts up with level for VS, but not really
% for d'. Not really very informative. 

inds = (stats_selected_BMF_stacked.modLevel==30 | stats_selected_BMF_stacked.modLevel==50 | ...
        stats_selected_BMF_stacked.modLevel==70 )  &  stats_selected_BMF_stacked.depthMod==1;
a1 =subplot(3,1,2);
boxplot(stats_selected_BMF_stacked.BMF(inds),{stats_selected_BMF_stacked.rationalisedType(inds) stats_selected_BMF_stacked.modLevel(inds)});
set(a1,'yscale','log','ytick',[25 100 500 2000]); ylabel('BMF-VS (Hz)');
ylim([20 4e3]); box off;  text(.75,5e3,'C. Distribution of best-VS by unit type and level','fontweight','bold');

inds = (stats_selected_BMFdp_stacked.modLevel==30 | stats_selected_BMFdp_stacked.modLevel==50 | ...
        stats_selected_BMFdp_stacked.modLevel==70) &  stats_selected_BMFdp_stacked.depthMod==1;

a2 = subplot(3,1,3);
boxplot(stats_selected_BMFdp_stacked.BMFdp(inds),{stats_selected_BMFdp_stacked.rationalisedType(inds) stats_selected_BMFdp_stacked.modLevel(inds)});
set(a2,'yscale','log','ytick',[25 100 500 2000]); xlabel('Unit type'); ylabel('BMF-d'' (Hz)');
ylim([20 4e3]); box off; text(.75,5e3,'D. Distribution of best d'' by unit type and level','fontweight','bold');


