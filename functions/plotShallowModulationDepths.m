            
%% ------------------------------------------------------------------------
% Looking at shallow modulation depths.

% There are relatively few examples at low modulation depths. 
% We just plot them all.

figh = figure('position',[100 -100 1000 800],'paperposition',[.5 .5 20 16]);
typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;
ymax  = 6;

subplot(2,2,1);
hold on;
tmp.dprimes = statsmat_modless1.dprimes(strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.5);
tmp.amfreqs = statsmat_modless1.usedamfreqs_dprime(strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.5);
tmp.dpover1 = cellfun(@(x) max(x>=1),tmp.dprimes)
cellfun(@(x,y) plot(x,y,'color',[0.7 0.7 0.7]),tmp.amfreqs(~tmp.dpover1),tmp.dprimes(~tmp.dpover1));
cellfun(@(x,y) plot(x,y,'linewidth',1),tmp.amfreqs(tmp.dpover1),tmp.dprimes(tmp.dpover1));
xlim([25 800]); ylim([-0.5 ymax]);
text(0,ymax*1.1,'a. ChS when modulation depth=50%','fontweight','bold');  
xlabel('Modulation frequency (Hz)','fontweight','bold'); ylabel(['Classifier ' measure],'fontweight','bold');
clear tmp;

subplot(2,2,2);
hold on;
tmp.dprimes = statsmat_modless1.dprimes(strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.2);
tmp.amfreqs = statsmat_modless1.usedamfreqs_dprime(strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.2);
tmp.dpover1 = cellfun(@(x) max(x>=1),tmp.dprimes)
cellfun(@(x,y) plot(x,y,'color',[0.7 0.7 0.7]),tmp.amfreqs(~tmp.dpover1),tmp.dprimes(~tmp.dpover1));
cellfun(@(x,y) plot(x,y,'linewidth',1),tmp.amfreqs(tmp.dpover1),tmp.dprimes(tmp.dpover1));
clear tmp;
text(0,ymax*1.1,'b. ChS when modulation depth=20%','fontweight','bold');  
xlabel('Modulation frequency (Hz)','fontweight','bold'); ylabel(['Classifier ' measure],'fontweight','bold');
xlim([25 800]); ylim([-0.2 ymax]);

subplot(2,2,3);
hold on;
tmp.dprimes = statsmat_modless1.dprimes(~strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.5);
tmp.amfreqs = statsmat_modless1.usedamfreqs_dprime(~strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.5);
tmp.type = statsmat_modless1.rationalisedType(~strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.5);
tmp.dpover1 = cellfun(@(x) max(x>=1),tmp.dprimes)
cellfun(@(x,y) plot(x,y,'color',[0.7 0.7 0.7]),tmp.amfreqs(~tmp.dpover1),tmp.dprimes(~tmp.dpover1));
cellfun(@(x,y) plot(x,y,'linewidth',1),tmp.amfreqs(tmp.dpover1),tmp.dprimes(tmp.dpover1));
lh = legend([tmp.type(~tmp.dpover1) tmp.type(tmp.dpover1)])
set(lh,'box','off','fontsize',8);
xlim([25 800]); ylim([-0.2 ymax]);
text(0,ymax*1.1,'c. Other types when modulation depth=50%','fontweight','bold');  
xlabel('Modulation frequency (Hz)','fontweight','bold'); ylabel(['Classifier ' measure],'fontweight','bold');
clear tmp;


subplot(2,2,4);
hold on;
tmp.dprimes = statsmat_modless1.dprimes(~strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.2);
tmp.amfreqs = statsmat_modless1.usedamfreqs_dprime(~strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.2);
tmp.type = statsmat_modless1.rationalisedType(~strcmp(statsmat_modless1.rationalisedType,'ChS') & [statsmat_modless1.depthMod{:}]==0.2);
tmp.dpover1 = cellfun(@(x) max(x>=1),tmp.dprimes)
cellfun(@(x,y) plot(x,y,'color',[0.7 0.7 0.7]),tmp.amfreqs(~tmp.dpover1),tmp.dprimes(~tmp.dpover1));
cellfun(@(x,y) plot(x,y,'linewidth',1),tmp.amfreqs(tmp.dpover1),tmp.dprimes(tmp.dpover1));
lh = legend([tmp.type(~tmp.dpover1) tmp.type(tmp.dpover1)])
set(lh,'box','off','fontsize',8);
xlim([25 800]); ylim([-0.2 ymax]);
text(0,ymax*1.1,'d. Other types when modulation depth=20%','fontweight','bold');  
xlabel('Modulation frequency (Hz)','fontweight','bold'); ylabel(['Classifier ' measure],'fontweight','bold');
clear tmp;

% print('-r600','-dtiff','figures\FigureS2_shallowmods.tif');
% saveas(gcf,'figures\FigureS2_shallowmods','fig')
% 
