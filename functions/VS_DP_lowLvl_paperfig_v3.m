function VS_DP_lowLvl_paperfig_v3(stats_oneperunit_stacked, stats_oneperunit_BMFdp_stacked, stats_oneperunit_BMF_stacked, coreoptions,dplim,measure,xlsxnames)
% VS and DP combined into one figure, for only the  low sound level
% dataset.

% ------------ Draft plot of vector strength for the paper ---------------

figh = figure('position',[100 -100 800 700],'paperposition',[.5 .5 16 14]);

typeColorMap;

typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;
dpmax = max(dplim);
dpmin = min(dplim);

%%%%%%%%%%%%%%%%%%%%%%%%%   D PRIME    %%%%%%%%%%%%%%%%%%%%%%%%%%

%% Summary plot for all unit types - one dataset per unit.  This doesn't need any extra options.
plotopts = {'figh',figh, ...
            'rowtitle','', ...
            'roworder', [5 6 1 2 3 4], ...
            'colorlist',typecolormap, ...
            'axeslist',{'DP'}, ...
            'axes_y',[0.6 0.3], ...
            'axes_x',[0.5 0.4], ...
            'offset_x',.5};        
po1 = plot_VS_dp_SPP5(stats_oneperunit_stacked,'rationalisedType',{},coreoptions,plotopts);
lh = legend(gca);
ylim([dpmin dpmax]);
set(lh,'position',[.8,.75,0.08,0.15]);
set(gca,'xtick',[-500:500:2000]); xlim([-500 2000]);
ylabel(['Mean ' measure],'fontweight','bold');
text(min(xlim),dpmax*1.1,'b. MTF from classifier','fontweight','bold');


%% MTF shape.

mtftypes = unique(stats_oneperunit_BMFdp_stacked.dprime_MTFType);
unittypes = typeorder;

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


%% Box plots of BMF distribution for the low level ones. 
% [0.65 0.1 0.3 0.3]
subplot('position',[0.1 0.6 0.3 0.3]);
hold on; 
line([0 7],[50 50],'linestyle','--','color',[.7 .7 .7]);
bh = boxplot(stats_oneperunit_BMFdp_stacked.BMFdp,stats_oneperunit_BMFdp_stacked.rationalisedType, ...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');
ylabel('BMF from classifier (Hz)','fontweight','bold');
ylim([0 1e3]); box off;  text(min(xlim),1.1e3,'a. BMF from classifier','fontweight','bold');
xlabel('Neuron type','fontweight','bold');
set(gca,'ytick',[100 200:200:1000]); 

% Option to write raw values to file.
if nargin>6
    xlsxfilename = xlsxnames{1};
    xlsxsheet = xlsxnames{2};

    % Output individual data values to an Excel file (PLOSBiol requirement).
    typed = cellfun(@(x) ~isempty(x), stats_oneperunit_BMFdp_stacked.rationalisedType);
    xlswrite(xlsxfilename, ...
        [ num2cell(stats_oneperunit_BMFdp_stacked.BMFdp(typed)') ...
         stats_oneperunit_BMFdp_stacked.rationalisedType(typed)'], ...
        [xlsxsheet 'A']);

end;

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

%% Box plots of BMF distribution for the low level ones. 
%stats_oneperunit_BMF_stacked.rationalisedType
%[0.65 0.1 0.3 0.3]
subplot('position',[0.1 0.1 0.3 0.3]);
hold on; 
line([0 7],[50 50],'linestyle','--','color',[.7 .7 .7]);
bh = boxplot(stats_oneperunit_BMF_stacked.BMF,stats_oneperunit_BMF_stacked.rationalisedType, ...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');

ylabel('BMF-VS (Hz)','fontweight','bold');
ylim([0 1e3]); box off;  text(min(xlim),1.1e3,'c. BMF-VS','fontweight','bold');
xlabel('Neuron type','fontweight','bold');
set(gca,'ytick',[100 200:200:1000]); 

% Option to write raw values to file.
if nargin>6
    % Output individual data values to an Excel file (PLOSBiol requirement).
    typed = cellfun(@(x) ~isempty(x), stats_oneperunit_BMF_stacked.rationalisedType);
    xlswrite(xlsxfilename, ...
        [ num2cell(stats_oneperunit_BMF_stacked.BMFdp(typed)') ...
         stats_oneperunit_BMF_stacked.rationalisedType(typed)'], ...
        [xlsxsheet 'C']);

end;



