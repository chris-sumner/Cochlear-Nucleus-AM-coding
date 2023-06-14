

%% -------- Statistical modelling figure for the paper  -----------

xSize = 27; ySize = 7; xLeft = (21-xSize)/2;yTop = (30-xSize)/2;
g=figure;set(g,'PaperUnits','centimeters');
set(g,'paperposition',[xLeft yTop xSize ySize],'units','centimeters','position',[2 4 xSize ySize],'DefaultAxesFontSize',8);%,'DefaultAxesColorOrder',colorOrder
axesw = 0.175;

typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;


% --------------- A. mean unnormalised MTF-d' functions --------

jj=2; ah(jj) = axes('position',[.05 .2 axesw .67]);hold on;
[~, colororder]= ismember(dPrime_reBMFdp.rowValues,typeorder)
[dh lh] = plotTable3(dPrime_reBMFdp,'sem','colorlist',typecolormap(colororder,:)); 
lstrs = get(lh,'String'); 
axis([0 1550 0 4.5]);
xlabel('f_m_o_d (Hz)','fontweight','bold'); 
ylabel('d''','fontweight','bold');
text(0, 4.5*1.1, 'a. Mean MTF-d''s','fontweight','bold');



% --------------- B. Peak d' box plot vs. type. ----------------

subplot('position', [0.27 0.2 axesw 0.67]);

bh = boxplot( stats_oneperunit_BMFdp_stacked.dprimes, ...`
    {stats_oneperunit_BMFdp_stacked.rationalisedType}, ...
     'grouporder',typeorder, 'colors',typecolormap, ...
     'boxstyle','filled');

box off; ylim([0 7]);
xlabel('Neuron type','fontweight','bold');
ylabel('d'' at BMF (Hz)','fontweight','bold');
text(0,max(ylim)*1.1,'b. Maximum d''','fontweight','bold');

% ------------ C. mean normalised mtf-d' functions -------------

[~, colororder]= ismember(dPrime_reBMFdp_reBestDP.rowValues,typeorder)
jj=4; ah(jj) = axes('position',[.51 .2 axesw .67]);hold on;
plot(flist,meanfModFn,'k:','linewidth',3);
[dh lh] = plotTable3(dPrime_reBMFdp_reBestDP,'sem','colorlist',typecolormap(colororder,:)); 
plot(flist,meanfModFn,'k:','linewidth',3);
lstrs = get(lh,'String'); 
lh= legend({'mean',lstrs{1:end-1}}); set(lh,'box','off');
axis([0 1550 0 1]);
xlabel('f_m_o_d (Hz)','fontweight','bold'); 
ylabel('Scaling factor','fontweight','bold');
text(0, 1.1, 'c. Normalised MTF-d''','fontweight','bold');


% --------- D. scatter plot showing the fit to all MTFs --------

jj=3; ah(jj) = axes('position',[.77 .2 axesw .67]);hold on;

for ti = 1:length(typeorder)
    inds = find(strcmp(stats_stacked_4models.rationalisedType,typeorder{ti}));
    plot(nltables.stim_upperbound.dprimes(inds),nlm.stim_theoreticalmax.fitted(inds),'.','color',typecolormap(ti,:));   
end;

line([-1 8],[-1 8],'linestyle','--','color','k','linewidth',1);
xlim([-.5 7]); ylim([-.5 7]);
xlabel('d'' (measured)','fontweight','bold'); 
ylabel('d'' (modeled)','fontweight','bold');
text(-.5, 7.7, 'd. Predicting all d''s','fontweight','bold');
set(gca,'xtick',[0:7],'ytick',[0:7]);
text(0,6,sprintf('R^2: %.3g',(nlm.stim_theoreticalmax.overallR)^2));