            
%% ------------------------------------------------------------------------
% Looking at data where frequency steps are narrower.

% How many of what steps there are. 
figure; 
freqsteps = unique(cellfun(@(x) x(2)-x(1),statsmat_selected.allamfreqs));
hist(cellfun(@(x) x(2)-x(1),statsmat_selected.allamfreqs),freqsteps)

statsmat_selected.amfreqstep = cellfun(@(x) x(2)-x(1),statsmat_selected.allamfreqs,'uni',false)

% 3 Ch-T datasets with 5Hz steps at very low modulation frequencies.
statsmat_selected_5Hzsteps = select_Datasets(statsmat_selected,'fn','amfreqstep==5');
% 21 datasets with 25Hz steps.
statsmat_selected_25Hzsteps = select_Datasets(statsmat_selected,'fn','amfreqstep==25');
% 196 datasets with 50Hz steps. 
statsmat_selected_50Hzsteps = select_Datasets(statsmat_selected,'fn','amfreqstep==50');


figh = figure('position',[100 -100 1000 800],'paperposition',[.5 .5 20 16]);
typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;


% 8 Ch-S and PBU neurons at 3 different signal levels when delta-f is 25 Hz. 
subplot(2,2,1);
% 7 ChS, one PBU.
unique([statsmat_selected_25Hzsteps.unitid{:}])
tmp.dprimes = statsmat_selected_25Hzsteps.dprimes;
tmp.amfreqs = statsmat_selected_25Hzsteps.usedamfreqs_dprime;
tmp.type = statsmat_selected_25Hzsteps.rationalisedType;
tmp.level = statsmat_selected_25Hzsteps.modLvl_rationalised;
uorder = [1:3 15:17 4:14 18:21]; % Order to display in - enables neat legend
legstr = {};
levelstyles = {'-',':','--'};

for ui = uorder
    hold on;        
    typecol = strcmp(tmp.type{ui},typeorder);
    levelstyle = levelstyles{find([30 50 70] == tmp.level{ui})};
    plot(tmp.amfreqs{ui},tmp.dprimes{ui},'linewidth',1,'linestyle',levelstyle,'color',typecolormap(typecol,:));
    legstr = {legstr{:} [tmp.type{ui} ' ' num2str(tmp.level{ui}) 'dB SPL']};
end;
lh = legend(legstr(1:6));
xlabel('Modulation frequency (Hz)','fontweight','bold'); ylabel(['Classifier ' measure],'fontweight','bold');
set(lh,'box','off');
xlim([0 1000]); ylim([-0.5 maxscore]);
text(0,maxscore*1.1,'a. Examples of discrimination with 25Hz steps','fontweight','bold');  
clear tmp;

% A single ChT neuron - run at 3 different signal levels. 
subplot(2,2,2);
hold on;
tmp.dprimes = statsmat_selected_5Hzsteps.dprimes;
tmp.amfreqs = statsmat_selected_5Hzsteps.usedamfreqs_dprime;
tmp.dpover1 = cellfun(@(x) max(x>=1),tmp.dprimes)
cellfun(@(x,y) plot(x,y,'color',[0.7 0.7 0.7]),tmp.amfreqs(~tmp.dpover1),tmp.dprimes(~tmp.dpover1));
cellfun(@(x,y) plot(x,y,'linewidth',1),tmp.amfreqs(tmp.dpover1),tmp.dprimes(tmp.dpover1));
lh = legend(cellfun(@(x) num2str(x), statsmat_selected_5Hzsteps.modLevel, 'UniformOutput', false));
set(lh,'box','off');
xlim([0 50]); ylim([-0.5 maxscore]);
text(0,maxscore*1.1,'b. Slow modulations in a single Ch-T neuron','fontweight','bold');  
xlabel('Modulation frequency (Hz)','fontweight','bold'); ylabel(['Classifier ' measure],'fontweight','bold');
clear tmp;

% print('-r300','-dtiff','figures\FigureS3_smallamsteps.tif');
% saveas(gcf,'figures\FigureS3_smallamsteps','fig')
