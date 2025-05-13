% Model-locking related statistics figure for the suplemental material. 
figh = figure('position',[100 -100 800 700],'paperposition',[.5 .5 16 14]);
typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;

plotopts = {'figh',figh, ...
            'rowtitle','a. Frequency dependence of Z_I_S_I', ...
            'roworder', [5 6 1 2 3 4], ...
            'colorlist',typecolormap, ...
            'axeslist',{'Z'}, ...
            'axes_y',[0.6 0.3], ...
            'axes_x',[0.5 0.4], ...
            'offset_x',.05};        
%            'offset_x',.1};        
po1 = plot_VS_dp_SPP5(stats_oneperunit_stacked,'rationalisedType',{},coreoptions,plotopts);

subplot('position',[0.55 0.6 0.3 0.3]);
hold on; 
line([0 7],[2 2],'linestyle','--','color',[.7 .7 .7]);
bh = boxplot(stats_oneperunit_stacked.Z_meanBelow1k,stats_oneperunit_stacked.rationalisedType, ...
    'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');
ylim([0 20]);xlim([0.5 6.5]);
text(min(xlim),max(ylim)*1.1,'b. Mean Z_I_S_Is for each dataset (f_m_o_d<1kHz)','fontweight','bold');
ylabel('Z_I_S_I','fontweight','bold'); xlabel('Neuron type'); box off;
%set(bh,'linewidth',1);

subplot('position',[0.05 0.15 0.3 0.3]);
hold on; 
line([0 7],[2 2],'linestyle','--','color',[.7 .7 .7]);
bh2 = boxplot(stats_oneperunit_stacked.VS_Isur_meanBelow1k,stats_oneperunit_stacked.rationalisedType,...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');
ylabel('VS','fontweight','bold'); xlabel('Neuron type');
ylim([0 1]); xlim([0.5 6.5]); box off;
text(min(xlim),max(ylim)*1.1,'c. Phase-locking after interval-shuffling','fontweight','bold');

subplot('position',[0.55 0.15 0.3 0.3]);
hold on; 
line([0 7],[13.8 13.8],'linestyle','--','color',[.7 .7 .7]);
bh = boxplot(stats_oneperunit_stacked.Rayleigh_Isur_meanBelow1k,stats_oneperunit_stacked.rationalisedType, ...
    'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');
ylim([0 30]);xlim([0.5 6.5]);
ylabel('Rayeigh','fontweight','bold'); xlabel('Neuron type'); box off;
text(min(xlim),max(ylim)*1.1,'d. Rayleigh statistic after interval-shuffling','fontweight','bold');

%print -dtiff -r300 FigureS4.tiff

