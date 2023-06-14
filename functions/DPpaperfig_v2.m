function VSpaperfig_v2(stats_oneperunit_stacked,stats_stacked, stats_oneperunit_BMFdp_stacked, stats_selected_BMFdp_stacked, coreoptions)
% VS paper fig v2.

% ------------ Draft plot of vector strength for the paper ---------------

figh = figure('position',[100 100 1200 700],'paperposition',[.5 .5 24 14]);

typeorder = {'ChS','ChT','PBU','PL','PLN','On'};
typecolorlist = {'m','g','r','k','b','c'};
typecolormap = [1 0 1; 0 1 0; 1 0 0; 0 0 0; 0 0 1; 0 1 1];
chScolorlist = { [.75 0 .75],[1 0 1 ],[1 .5 1]};
PLcolorlist = { [0 0 0 ],[.4  .4 .4 ],[.7 .7 .7]};

%% Summary plot for all unit types - one dataset per unit.  This doesn't need any extra options.
plotopts = {'figh',figh, ...
            'rowtitle','A. MTF-d'' at low sound levels', ...
            'roworder', [1 2 4 5 6 3], ...
            'colorlist',{'m','g','r','k','b','c'}, ...
            'axeslist',{'DP'}, ...
            'axes_y',[0.55 0.3], ...
            'axes_x',[0.2 0.18], ...
            'offset_x',.05};        
po1 = plot_VS_dp_SPP4(stats_oneperunit_stacked,'rationalisedType',{},coreoptions,plotopts);
lh = legend(gca);
set(lh,'position',[.23,.55,0.08,0.15])
set(gca,'xtick',[-1000:500:2000]); xlim([-500 2000]);
ylabel('Mean d''','fontweight','bold');


%% Inset plot showing the normalised function. 
plotopts([3 4 9:18]) = {  'rowtitle','', ...
                      'axeslist', {'DPn'}, ...
                      'axes_y',[0.75 0.1], ...
                      'axes_x',[0.15 0.08], ...
                      'offset_x',   0.23, ...
                      'sem','nosem'};
po1B = plot_VS_dp_SPP4(stats_oneperunit_stacked,'rationalisedType',{},coreoptions,plotopts);
axes(po1B.ahs); 
legend off;
xlabel(''); 
ylabel('d'' / d''_m_a_x','fontweight','bold','fontsize',8);
set(gca,'xticklabel',[],'ytick',[0.5 1],'yticklabel',[0.5 1]);

%[0.05 0.1 0.4 0.3]

%% ChS neuron type and looking across level instead.
options = {'fn','depthMod==1','fn','strcmp(rationalisedType,''ChS'')'};   
plotopts([3:6 7:18]) = {  'rowtitle','d. ChS neurons', ...
                      'colorlist',chScolorlist, ...
                      'roworder',  [1:3], ...
                      'axeslist', {'DP'}, ...
                      'axes_y',[0.1 0.3], ...
                      'axes_x',[0.2 0.18], ...
                      'offset_x',   0.05, ...
                      'sem','sem'};        
po2 = plot_VS_dp_SPP4(stats_stacked,'modLvl_rationalised',options,coreoptions,plotopts);
ylabel('Mean d''','fontweight','bold');


%% PL neuron type and looking across level.
options = {'fn','depthMod==1','fn','strcmp(rationalisedType,''PL'')'};     
plotopts([3:4 5 6 15:16]) = {  'rowtitle','e. PL neurons', ...
                           'colorlist',PLcolorlist, ...
                           'offset_x', 0.3};               
po3 = plot_VS_dp_SPP4(stats_stacked,'modLvl_rationalised',options,coreoptions,plotopts);
ylabel('Mean d''','fontweight','bold');


%% Box plots of BMF distribution for the low level ones. 

subplot('position',[0.37 0.55 0.11 0.3]);
hold on; 
line([0 7],[50 50],'linestyle','--','color',[.7 .7 .7]);
bh = boxplot(stats_oneperunit_BMFdp_stacked.BMFdp,stats_oneperunit_BMFdp_stacked.rationalisedType, ...
     'colors',[typecolorlist{:}],'grouporder',typeorder,'boxstyle','filled');
ylabel('BMF-d'' (Hz)','fontweight','bold');
ylim([20 1e3]); box off;  text(.75,1.1e3,'b. BMF-d''','fontweight','bold');
xlabel('Neuron type','fontweight','bold');
set(gca,'ytick',[100 200:200:1000],'xticklabel',[]); 


%% MTF shape.


% ----------------- SPLIT BY SOUND LEVEL ---------------

mtftypes = unique(stats_selected_BMFdp_stacked.dprime_MTFType);
unittypes = {'ChS','ChT','PBU','PL','PLN','On'};
levels = [30 50 70];

% VS table.
coli = 1; vs_nos = []; dp_nos = [];
for ui = 1:length(unittypes)
    for li = 1:length(levels)
       vsmtftype_unittype = stats_selected_BMFdp_stacked.VS_MTFType( ...
            strcmp(stats_selected_BMFdp_stacked.rationalisedType,unittypes{ui}) & ...
            stats_selected_BMFdp_stacked.modLevel>=(levels(li)-10) & stats_selected_BMFdp_stacked.modLevel<(levels(li)+10));
 
       vs_nos(coli,:) =  ...
           [ sum(strcmp(vsmtftype_unittype,'lowpass')) ...
            sum(strcmp(vsmtftype_unittype,'bandpass')) ...
            sum(~strcmp(vsmtftype_unittype,'lowpass') & ~strcmp(vsmtftype_unittype,'bandpass'))];
        
      dpmtftype_unittype = stats_selected_BMFdp_stacked.dprime_MTFType( ...
          strcmp(stats_selected_BMFdp_stacked.rationalisedType,unittypes{ui}) & ...
          stats_selected_BMFdp_stacked.modLevel>=(levels(li)-10) & stats_selected_BMFdp_stacked.modLevel<(levels(li)+10));
 
       dp_nos(coli,:) =  ...
           [ sum(strcmp(dpmtftype_unittype,'lowpass')) ...
            sum(strcmp(dpmtftype_unittype,'bandpass')) ...
            sum(~strcmp(dpmtftype_unittype,'lowpass') & ~strcmp(dpmtftype_unittype,'bandpass'))];

        coli = coli+1;
    end;
    % create a gap between types.
    coli=coli+1;
end;

%                       'axes_y',[0.55 0.3], ...
%                       'axes_x',[0.2 0.15], ...
%                       'offset_x',   0.6, ...


maxv = max([sum(dp_nos,2);  sum(dp_nos,2)]);
subplot('position', [0.55 0.55 0.41 0.3]); hold on;
for li = 1:size(dp_nos,1)
    bh = bar([li li+1],[dp_nos(li,:); 0 0 0],'stacked'); 
    typei = 1+floor(li/4);
    lvli = rem(li,4);
    thiscolor = typecolormap(typei,:);
    thiscolor(thiscolor==0) = 0.25*lvli;
    set(bh(1),'facecolor',thiscolor);
    set(bh(2),'facecolor','w');
    set(bh(3),'facecolor',[.7 .7 .7]);
end;
box off;
ylabel('# neurons','fontweight','bold');
set(gca,'xticklabel',unittypes,'xtick',[2 6 10 14 18 22]);
text(0,maxv*1.1*1.1,'c. MTF-d'' classification','fontweight','bold');
xlabel('Neuron type/Sound level','fontweight','bold');


text(4.5,25,'30','fontsize',7); text(5.5,37,'50','fontsize',7); 
text(6.5,41,'70','fontsize',7); text(5,47,'dB SPL','fontsize',7);
text(2.9,8,'lowpass','rotation',90,'fontsize',7);
ylim([0 maxv*1.1])

%% Peak VS against level.
grouporder = {'ChS,30','ChS,50','ChS,70','ChT,30','ChT,50','ChT,70', ...
      'PBU,30','PBU,50','PBU,70','PL,30','PL,50','PL,70', ...
      'PLN,30','PLN,50','PLN,70', 'On,30','On,50','On,70'};

subplot('position', [0.58 0.1 0.38 0.3]);
boxplot( stats_selected_BMFdp_stacked.dprimes, ...`
    {stats_selected_BMFdp_stacked.rationalisedType stats_selected_BMFdp_stacked.modLvl_rationalised}, ...
    'factorgap',[5 1]);

bh = boxplot( stats_selected_BMFdp_stacked.dprimes, ...`
    {stats_selected_BMFdp_stacked.rationalisedType stats_selected_BMFdp_stacked.modLvl_rationalised}, ...
     'grouporder',grouporder, 'colors',char([typecolorlist{:}]'*[1 1 1])','factorgap',[5 1], ...
     'boxstyle','filled');


% Delete some of the messy type labels. 
kids =  get(gca,'children');
delete(kids.Children([19:3:35 21:3:36]));
box off;
xlabel('Neuron type/sound level','fontweight','bold');
ylabel('d'' at BMF (Hz)','fontweight','bold');
text(0,max(ylim)*1.1,'f. Maximum d''','fontweight','bold');












