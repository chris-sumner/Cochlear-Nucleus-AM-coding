function plotMTFStats(stats_selected_BMF_stacked,stats_selected_BMFdp_stacked,coreoptions,options,plotopt);

% Using the BMF stacked stats to make tables and histograms of the
% BMF for both statistics. 
rowandcol =  {'rationalisedType','usedamfreqs_reBMF'};
stats.VS_MTFCutOffBW = makePopnTable(stats_selected_BMF_stacked, rowandcol{:} ,'VS_MTFCutOffBW',coreoptions{:}, ...
    options{:},'hist',[50:200:2050],'fn','Rayleigh>13.8');
stats.VS_LoCornerHz = makePopnTable(stats_selected_BMF_stacked, rowandcol{:} ,'VS_LoCornerHz',coreoptions{:}, ...
    options{:},'hist',[50:200:2050],'fn','Rayleigh>13.8');
stats.VS_HiCornerHz = makePopnTable(stats_selected_BMF_stacked, rowandcol{:} ,'VS_HiCornerHz',coreoptions{:}, ...
    options{:},'hist',[50:200:2050],'fn','Rayleigh>13.8');
% Threshold statistics.
stats.VS_MTFAbsThrBW = makePopnTable(stats_selected_BMF_stacked, rowandcol{:} ,'VS_MTFAbsThrBW',coreoptions{:}, ...
    options{:},'hist',[50:200:2050],'fn','Rayleigh>13.8');
stats.VS_LoFAbsThr = makePopnTable(stats_selected_BMF_stacked, rowandcol{:} ,'VS_LoFAbsThr',coreoptions{:}, ...
    options{:},'hist',[50:200:2050],'fn','Rayleigh>13.8');
stats.VS_HiFAbsThr = makePopnTable(stats_selected_BMF_stacked, rowandcol{:} ,'VS_HiFAbsThr',coreoptions{:}, ...
    options{:},'hist',[50:200:2050],'fn','Rayleigh>13.8');


rowandcol =  {'rationalisedType','usedamfreqs_reBMFdp'};
stats.dprime_MTFCutOffBW = makePopnTable(stats_selected_BMFdp_stacked, rowandcol{:} ,'dprime_MTFCutOffBW',coreoptions{:}, ...
    options{:},'hist',[50:200:2050]);
stats.dprime_LoCornerHz = makePopnTable(stats_selected_BMFdp_stacked, rowandcol{:} ,'dprime_LoCornerHz',coreoptions{:}, ...
    options{:},'hist',[50:200:2050]);
stats.dprime_HiCornerHz = makePopnTable(stats_selected_BMFdp_stacked, rowandcol{:} ,'dprime_HiCornerHz',coreoptions{:}, ...
    options{:},'hist',[50:200:2050]);
% Threshold statistics.
stats.dprime_MTFAbsThrBW = makePopnTable(stats_selected_BMFdp_stacked, rowandcol{:} ,'dprime_MTFAbsThrBW',coreoptions{:}, ...
    options{:},'hist',[50:200:2050]);
stats.dprime_LoFAbsThr = makePopnTable(stats_selected_BMFdp_stacked, rowandcol{:} ,'dprime_LoFAbsThr',coreoptions{:}, ...
    options{:},'hist',[50:200:2050]);
stats.dprime_HiFAbsThr = makePopnTable(stats_selected_BMFdp_stacked, rowandcol{:} ,'dprime_HiFAbsThr',coreoptions{:}, ...
    options{:},'hist',[50:200:2050]);



%%  Looking at the -3dB statistics.

% Bandwidths are a bit different.
figure('position',[100 100 800 800],'paperposition',[.5 .5 16 16]); 

subplot(3,2,1);
plot([50:200:2050],squeeze(stats.VS_MTFCutOffBW.hist)','o-');
xlabel('30% Bandwidth-VS (Hz)'); ylabel('# units'); %legend(stats.VS_MTFCutOffBW.rowValues);
ylim([0 20]); xlim([0 2100]); box off; title('MTF-VS shape statistics');


subplot(3,2,2);
plot([50:200:2050],squeeze(stats.dprime_MTFCutOffBW.hist)','o-');
xlabel('30% Bandwidth-d'' (Hz)'); ylabel('# units'); %legend(stats.dprime_MTFCutOffBW.rowValues);
ylim([0 20]);xlim([0 2100]);box off; title('MTF-d'' shape statistics');


% Low corner frequencies are quite a lot lower and less numerous.
subplot(3,2,3);
plot([50:200:2050],squeeze(stats.VS_LoCornerHz.hist)','o-');
xlabel('70% Lower corner frequency-VS (Hz)'); ylabel('# units'); %legend(stats.VS_LoCornerHz.rowValues);
ylim([0 20]);xlim([0 2100]);box off; 

subplot(3,2,4);
plot([50:200:2050],squeeze(stats.dprime_LoCornerHz.hist)','o-');
xlabel('70% Lower corner frequency-d'' (Hz)'); ylabel('# units'); %legend(stats.dprime_LoCornerHz.rowValues);
ylim([0 20]);xlim([0 2100]);box off; 

subplot(3,2,5);
plot([50:200:2050],squeeze(stats.VS_HiCornerHz.hist)','o-');
xlabel('70% Upper corner frequency-VS (Hz)'); ylabel('# units'); 
ylim([0 22]);xlim([0 2100]);box off; 

subplot(3,2,6);
plot([50:200:2050],squeeze(stats.dprime_HiCornerHz.hist)','o-');
xlabel('70% Upper corner frequency-d'' (Hz)'); ylabel('# units'); 
ylim([0 22]);xlim([0 2100]);box off; 
lh = legend(stats.VS_HiCornerHz.rowValues); set(lh,'box','off')

if nargin<5
    print -dtiff -r150 mtffigures\MTF_CornerFreqStatDists.tiff
end;

%% Looking at the "threshold" statistics 


figure('position',[100 100 800 800],'paperposition',[.5 .5 16 16]); 

subplot(3,2,1);
plot([50:200:2050],squeeze(stats.VS_MTFAbsThrBW.hist)','o-');
xlabel('Threshold (1)-VS bandwidth (Hz)'); ylabel('# units'); %legend(stats.VS_LoCornerHz.rowValues);
ylim([0 20]); xlim([0 2100]); box off; title('MTF-VS threshold statistics');

subplot(3,2,2);
plot([50:200:2050],squeeze(stats.dprime_MTFAbsThrBW.hist)','o-');
xlabel('Threshold (1)-d'' bandwidth (Hz)'); ylabel('# units'); %legend(stats.VS_LoCornerHz.rowValues);
ylim([0 20]);xlim([0 2100]);box off; title('MTF-d'' threshold statistics');

% Low corner frequencies are quite a lot lower and less numerous.
subplot(3,2,3);
plot([50:200:2050],squeeze(stats.VS_LoFAbsThr.hist)','o-');
xlabel('Threshold (1)-VS lower frequency (Hz)'); ylabel('# units'); %legend(stats.VS_LoCornerHz.rowValues);
ylim([0 20]);xlim([0 2100]);box off; 

subplot(3,2,4);
plot([50:200:2050],squeeze(stats.dprime_LoFAbsThr.hist)','o-');
xlabel('Threshold (1)-d'' lower frequency (Hz)'); ylabel('# units'); %legend(stats.VS_LoCornerHz.rowValues);
ylim([0 20]);xlim([0 2100]);box off; 

subplot(3,2,5);
plot([50:200:2050],squeeze(stats.VS_HiFAbsThr.hist)','o-');
xlabel('Threshold (0.2)-VS upper frequency (Hz)'); ylabel('# units'); 
ylim([0 22]);xlim([0 2100]);box off; 

subplot(3,2,6);
plot([50:200:2050],squeeze(stats.dprime_HiFAbsThr.hist)','o-');
xlabel(' Threshold (1)-d'' upper frequency (Hz)'); ylabel('# units'); 
ylim([0 22]);xlim([0 2100]);box off; 
lh = legend(stats.dprime_HiFAbsThr.rowValues); set(lh,'box','off')

if nargin<5
    print -dtiff -r150 mtffigures\MTF_CutOffStatDists.tiff
end;


%% How statistics change with level. 

rowandcol =  {'rationalisedType','modLvl_rationalised'};

stats.dprimeLvl_MTFCutOffBW = makePopnTable(stats_selected_BMFdp_stacked, rowandcol{:} ,'dprime_MTFCutOffBW',coreoptions{:}, ...
    options{1:2},'hist',[50:200:2050]);
stats.VSLvl_MTFCutOffBW = makePopnTable(stats_selected_BMF_stacked, rowandcol{:} ,'VS_MTFCutOffBW',coreoptions{:}, ...
    options{1:2},'hist',[50:200:2050]);

figure('position',[100 100 800 800],'paperposition',[.5 .5 16 16]); 
subplot(2,1,1);
boxplot( stats_selected_BMFdp_stacked.dprime_MTFCutOffBW, ...`
    {stats_selected_BMFdp_stacked.rationalisedType stats_selected_BMFdp_stacked.modLvl_rationalised}, ...
    'factorgap',[5 1]);

% Delete some of the messy type labels. 
kids =  get(gca,'children');
%delete(kids.Children([19:3:35 21:3:36]));
box off;
xlabel('Neuron type/sound level');
ylabel('d'' cut-off (-10dB) frequency (Hz)');
text(1,2100,'A. d'' cutoff as a function of sound level and unit type','fontweight','bold');

subplot(2,1,2);
boxplot( stats_selected_BMF_stacked.VS_MTFCutOffBW, ...`
    {stats_selected_BMF_stacked.rationalisedType stats_selected_BMF_stacked.modLvl_rationalised}, ...
    'factorgap',[5 1]);

% Delete some of the messy type labels. 
kids =  get(gca,'children');
delete(kids.Children([19:3:35 21:3:36]));
box off;
xlabel('Neuron type/sound level');
ylabel('VS cut-off (-10dB) frequency (Hz)');
text(1,2100,'B. VS cutoff as a function of sound level and unit type','fontweight','bold');

if nargin<5
    print -dtiff -r150 mtffigures\MTF_CutOffStatsByLvl.tiff
end;





