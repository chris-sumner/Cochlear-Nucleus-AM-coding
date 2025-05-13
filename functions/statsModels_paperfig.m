

%% -------- Statistical modelling figure for the paper  -----------

xSize = 27; ySize = 7; xLeft = (21-xSize)/2;yTop = (30-xSize)/2;
g=figure;set(g,'PaperUnits','centimeters');
set(g,'paperposition',[xLeft yTop xSize ySize],'units','centimeters','position',[2 4 xSize ySize],'DefaultAxesFontSize',8);%,'DefaultAxesColorOrder',colorOrder
axesw = 0.175;

typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;

% extract these from the theoretical model
flist = nlm.stim_theoreticalmax.flist;
meanfModFn =  nlm.stim_theoreticalmax.normalised_fitted_fn;

% N.B. the fitted population MTF is normalised to 1 but here we scale it 
% to fit over the functions to emphasise the fit to shape. This change in 
% scale is of no functional importance.
meanfModFn = meanfModFn*mean(dPrime_reBMFdp_reBestDP.mean(:,1))

% --------------- A. mean unnormalised MTF-d' functions --------


jj=2; ah(jj) = axes('position',[.05 .2 axesw .67]);hold on;
[~, colororder]= ismember(dPrime_reBMFdp.rowValues,typeorder)
[dh lh] = plotTable3(dPrime_reBMFdp,'sem','colorlist',typecolormap(colororder,:)); 
lstrs = get(lh,'String'); 
axis([0 1550 0 maxscore]);
xlabel('f_m_o_d (Hz)','fontweight','bold'); 
ylabel(['Classifier ' measure],'fontweight','bold');
text(0, maxscore*1.1, 'a. Mean classifier MTFs','fontweight','bold');

% --------------- B. Peak d' box plot vs. type. ----------------

subplot('position', [0.27 0.2 axesw 0.67]);

bh = boxplot( stats_oneperunit_BMFdp_stacked.dprimes, ...`
    {stats_oneperunit_BMFdp_stacked.rationalisedType}, ...
     'grouporder',typeorder, 'colors',typecolormap, ...
     'boxstyle','filled');

% Output individual data values to an Excel file (PLOSBiol requirement).
typed = cellfun(@(x) ~isempty(x), stats_oneperunit_BMFdp_stacked.rationalisedType);
xlswrite(xlsxfilename, ...
    [ num2cell(stats_oneperunit_BMFdp_stacked.dprimes(typed)') ...
    stats_oneperunit_BMFdp_stacked.rationalisedType(typed)'], ...
    xlsxsheet);
 
box off; ylim([0 maxscore]);
xlabel('Neuron type','fontweight','bold');
ylabel(['Classifier ' measure ' at BMF (Hz)'],'fontweight','bold');
text(0,max(ylim)*1.1,['b. Maximum classifier ' measure],'fontweight','bold');



% ------------ C. mean normalised mtf-d' functions -------------


% N.B. This is relative to the peak d' but not frequency normalised to BMF.
[~, colororder]= ismember(dPrime_reBMFdp_reBestDP.rowValues,typeorder)
jj=4; ah(jj) = axes('position',[.51 .2 axesw .67]);hold on;
plot(flist,meanfModFn,'k:','linewidth',3);
[dh lh] = plotTable3(dPrime_reBMFdp_reBestDP,'sem','colorlist',typecolormap(colororder,:)); 
plot(flist,meanfModFn,'k:','linewidth',3);

ymax = 1.02;
lstrs = get(lh,'String'); 
lh= legend({'best-fit',lstrs{1:end-1}}); set(lh,'box','off');
axis([0 1550 0 ymax]);
xlabel('f_m_o_d (Hz)','fontweight','bold'); 
ylabel(['Relative ' measure ' (unitless)'] ,'fontweight','bold');
text(0, 1.1*ymax, 'c. Normalised classifier MTFs','fontweight','bold');


% --------- D. scatter plot showing the fit to all MTFs --------

jj=3; ah(jj) = axes('position',[.77 .2 axesw .67]);hold on;

for ti = 1:length(typeorder)
    inds0 = find(~isinf(stats_stacked_4models.dprimes));
    inds = find(strcmp(stats_stacked_4models.rationalisedType(inds0),typeorder{ti}));
    inds = inds(inds<length(nlm.stim_theoreticalmax.measureV));
    plot(nlm.stim_theoreticalmax.measureV(inds),nlm.stim_theoreticalmax.predmeasureV(inds),'.','color',typecolormap(ti,:));   
end;

line([-1 8],[-1 8],'linestyle','--','color','k','linewidth',1);
xlabel(['Classifier ' measure ' (data)'],'fontweight','bold'); 
ylabel(['Classifier ' measure ' (model)'],'fontweight','bold');
xlim([mmin mmax]); ylim([pmin pmax]);
text(mmax*0.05,pmax*.9,sprintf('R^2: %.3g',(nlm.stim_theoreticalmax.overallR)^2));
text(mmin, pmax*1.1, 'd. Predicting classifier performance','fontweight','bold');


