function VS_DP_lowLvl_paperfig_v2(stats_oneperunit_stacked, stats_oneperunit_BMFdp_stacked, stats_oneperunit_BMF_stacked, coreoptions)
% VS and DP combined into one figure, for only the  low sound level
% dataset.

% ------------ Draft plot of vector strength for the paper ---------------

figh = figure('position',[100 -100 800 700],'paperposition',[.5 .5 16 14]);

typeColorMap;

typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;

% typeorder = {'ChS','ChT','PBU','PL','PLN','On'};
% typecolorlist = {'m','g','r','k','b','c'};
% typecolormap = [1 0 1; 0 1 0; 1 0 0; 0 0 0; 0 0 1; 0 1 1];
% chScolorlist = { [.75 0 .75],[1 0 1 ],[1 .5 1]};
% PLcolorlist = { [0 0 0 ],[.4  .4 .4 ],[.7 .7 .7]};
            

%%%%%%%%%%%%%%%%%%%%%%%%%   D PRIME    %%%%%%%%%%%%%%%%%%%%%%%%%%

%% Summary plot for all unit types - one dataset per unit.  This doesn't need any extra options.
plotopts = {'figh',figh, ...
            'rowtitle','b. MTF-d'' ', ...
            'roworder', [5 6 1 2 3 4], ...
            'colorlist',typecolormap, ...
            'axeslist',{'DP'}, ...
            'axes_y',[0.6 0.3], ...
            'axes_x',[0.5 0.4], ...
            'offset_x',.5};        
%            'offset_x',.1};        
po1 = plot_VS_dp_SPP5(stats_oneperunit_stacked,'rationalisedType',{},coreoptions,plotopts);
lh = legend(gca);
%set(lh,'position',[.35,.7,0.08,0.15]);
set(lh,'position',[.8,.75,0.08,0.15]);
set(gca,'xtick',[-1000:500:2000]); xlim([-500 2000]);
ylabel('Mean d''','fontweight','bold');


%% MTF shape.


mtftypes = unique(stats_oneperunit_BMFdp_stacked.dprime_MTFType);
unittypes = typeorder; %{'ChS','ChT','PBU','PL','PLN','On'};

%DP table ad VS table (no need to do again).
coli = 1; vs_nos = []; dp_nos = [];
for ui = 1:length(unittypes)
        li = 1;
        
        vsmtftype_unittype = stats_oneperunit_BMF_stacked.VS_MTFType( ...
            strcmp(stats_oneperunit_BMF_stacked.rationalisedType,unittypes{ui}));
 
       vs_nos(coli,:) =  ...
           [ sum(strcmp(vsmtftype_unittype,'lowpass')) ...
            sum(strcmp(vsmtftype_unittype,'bandpass')) ...
            sum(~strcmp(vsmtftype_unittype,'lowpass') & ~strcmp(vsmtftype_unittype,'bandpass'))];
        
      dpmtftype_unittype = stats_oneperunit_BMFdp_stacked.dprime_MTFType( ...
          strcmp(stats_oneperunit_BMFdp_stacked.rationalisedType,unittypes{ui}));
 
       dp_nos(coli,:) =  ...
           [ sum(strcmp(dpmtftype_unittype,'lowpass')) ...
            sum(strcmp(dpmtftype_unittype,'bandpass')) ...
            sum(~strcmp(dpmtftype_unittype,'lowpass') & ~strcmp(dpmtftype_unittype,'bandpass'))];

    coli=coli+1;
end;


% Doesn't really make the point so not included.
% maxv = max([sum(dp_nos,2);  sum(dp_nos,2)]);
% subplot('position', [0.4 0.6 0.23 0.3]); hold on;
% for li = 1:size(dp_nos,1)
%     bh = bar([li li+1],[dp_nos(li,:); 0 0 0],'stacked'); 
%     typei = li;
%     lvli = 1;
%     thiscolor = typecolormap(typei,:);
%     thiscolor(thiscolor==0) = 0.25*lvli;
%     set(bh(1),'facecolor',thiscolor,'edgecolor',thiscolor,'linewidth',1);
%     set(bh(2),'facecolor','w','edgecolor',thiscolor,'linewidth',1);
%     set(bh(3),'facecolor',[.7 .7 .7],'edgecolor',thiscolor,'linewidth',1);
% end;
% box off;
% ylabel('# neurons','fontweight','bold');
% set(gca,'xticklabel',unittypes,'xtick',[1:6]);
% text(0,maxv*1.1*1.1,'b. MTF-d'' classification','fontweight','bold');
% xlabel('Neuron type','fontweight','bold');
% 
% text(1.9,8,'lowpass','rotation',90,'fontsize',8);
% text(1.9,32,'b/pass','rotation',90,'fontsize',8);
% ylim([0 maxv*1.1])
% xlim([.25 6.75]);



%% Box plots of BMF distribution for the low level ones. 
% [0.65 0.1 0.3 0.3]
subplot('position',[0.1 0.6 0.3 0.3]);
hold on; 
line([0 7],[50 50],'linestyle','--','color',[.7 .7 .7]);
bh = boxplot(stats_oneperunit_BMFdp_stacked.BMFdp,stats_oneperunit_BMFdp_stacked.rationalisedType, ...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');
set(bh,'linewidth',1);
bh = boxplot(stats_oneperunit_BMFdp_stacked.BMFdp,stats_oneperunit_BMFdp_stacked.rationalisedType, ...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');
ylabel('BMF-d'' (Hz)','fontweight','bold');
ylim([0 1e3]); box off;  text(min(xlim),1.1e3,'a. BMF-d''','fontweight','bold');
xlabel('Neuron type','fontweight','bold');
set(gca,'ytick',[100 200:200:1000]); 








%%%%%%%%%%%%%%%%%%%%%%%%%%      VS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Summary plot for all unit types - one dataset per unit.  This doesn't need any extra options.
plotopts = {'figh',figh, ...
            'rowtitle','d. MTF-VS', ...
            'roworder', [5 6 1 2 3 4], ...
            'colorlist',typecolorlist, ...
            'axeslist',{'VS'}, ...
            'axes_y',[0.1 0.3], ...
            'axes_x',[0.5 0.4], ...
            'offset_x',.5};        
po1 = plot_VS_dp_SPP5(stats_oneperunit_stacked,'rationalisedType',{},coreoptions,plotopts);
legend off;
ylabel('Mean VS','fontweight','bold');

%% VS-MTF type - didn't look so good so I didn't use. 

% maxv = max([sum(vs_nos,2);  sum(vs_nos,2)]);
% subplot('position', [0.4 0.1 0.23 0.3]); hold on;
% for li = 1:size(vs_nos,1)
%     bh = bar([li li+1],[vs_nos(li,:); 0 0 0],'stacked'); 
%     typei = li;
%     lvli = 1;
%     thiscolor = typecolormap(typei,:);
%     thiscolor(thiscolor==0) = 0.25*lvli;
%     set(bh(1),'facecolor',thiscolor,'edgecolor',thiscolor,'linewidth',1);
%     set(bh(2),'facecolor','w','edgecolor',thiscolor,'linewidth',1);
%     set(bh(3),'facecolor',[.7 .7 .7],'edgecolor',thiscolor,'linewidth',1);
% end;
% box off;
% ylabel('# neurons','fontweight','bold');
% set(gca,'xticklabel',unittypes,'xtick',[1:6]);
% text(0,maxv*1.1*1.1,'e. MTF-VS'' classification','fontweight','bold');
% xlabel('Neuron type','fontweight','bold');
% 
% text(1.9,8,'lowpass','rotation',90,'fontsize',8);
% text(1.9,36,'b/pass','rotation',90,'fontsize',8);
% ylim([0 maxv*1.1])
% xlim([.25 6.75]);



%% Box plots of BMF distribution for the low level ones. 
%stats_oneperunit_BMF_stacked.rationalisedType
%[0.65 0.1 0.3 0.3]
subplot('position',[0.1 0.1 0.3 0.3]);
hold on; 
line([0 7],[50 50],'linestyle','--','color',[.7 .7 .7]);
bh = boxplot(stats_oneperunit_BMF_stacked.BMF,stats_oneperunit_BMF_stacked.rationalisedType, ...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');
set(bh,'linewidth',1);
bh = boxplot(stats_oneperunit_BMF_stacked.BMF,stats_oneperunit_BMF_stacked.rationalisedType, ...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');

ylabel('BMF-VS (Hz)','fontweight','bold');
ylim([0 1e3]); box off;  text(min(xlim),1.1e3,'c. BMF-VS','fontweight','bold');
xlabel('Neuron type','fontweight','bold');
set(gca,'ytick',[100 200:200:1000]); 






