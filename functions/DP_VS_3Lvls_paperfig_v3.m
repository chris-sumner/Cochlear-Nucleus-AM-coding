function DP_VS_3Lvls_paperfig_v1(stats_stacked,  stats_selected_BMFdp_stacked, stats_selected_BMF_stacked, coreoptions,dplim, measure)

dpmax = max(dplim);
dpmin = min(dplim);

% ------------ MTF shape stats wqith level for the paper ---------------

figh = figure('position',[100 100 1200 800],'paperposition',[.5 .5 24 16]);

typeColorMap;

typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;
roworder = [5 6 1 2 3 4];
chScolorlist = {typecolormap(3,:)*1.2 typecolormap(3,:)*.8 typecolormap(3,:)*0.5};
chScolorlist{1}(1) = 1;

PLcolorlist = {typecolormap(1,:)*1.2 typecolormap(1,:)*.8 typecolormap(1,:)*0.5};
PLcolorlist{1}(3) = 1;

% Look for the threshold number of neurons to plot for means.
ind = find(strcmp('minn',coreoptions));           % minimum number of units to include condition.
if ~isempty(ind)
    min_n = coreoptions{ind+1};
end;

plotopts = {'figh',figh, ...
            'rowtitle','a. MTF from classifier at low sound levels', ...
            'roworder', roworder, ...
            'colorlist',typecolormap, ...
            'axeslist',{'DP'}, ...
            'axes_y',[0.6 0.27], ...
            'axes_x',[0.2 0.18], ...
            'offset_x',.05, ...
            'xlim',[-500 1500], ...
            'sem','sem'};        

%%%%%%%%%%%%%%%%%%%%        D-PRIME     %%%%%%%%%%%%%%%%%%%%%%
 
y1 = 0.1; y2 = 0.63;
h1 = 0.27;


%% ChS neuron type and looking across level instead.
options = {'fn','strcmp(rationalisedType,''ChS'')'};   
plotopts([3:6 7:16]) = {  'rowtitle','', ...
                      'colorlist',chScolorlist, ...
                      'roworder',  [1:3], ...
                      'axeslist', {'DP'}, ...
                      'axes_y',[y2 h1], ...
                      'axes_x',[0.13 0.18], ...
                      'offset_x',   0.54};     
                  
po2 = plot_VS_dp_SPP4(stats_stacked,'modLvl_rationalised',options,coreoptions,plotopts);
ylabel(['Classifier ' measure],'fontweight','bold');
xlim([-500 1500]);
ylim([dpmin dpmax]);
text(min(xlim),dpmax*1.1,'c. Mean ChS classifier MTF','fontweight','bold');


%% PL neuron type and looking across level.
options = {'fn','strcmp(rationalisedType,''PL'')'};     
plotopts([3:4 5 6 15:18]) = {  'rowtitle','', ...
                           'colorlist',PLcolorlist, ...
                           'offset_x', 0.78, ...
                           'xlim',[-1000 1500]};               
po3 = plot_VS_dp_SPP4(stats_stacked,'modLvl_rationalised',options,coreoptions,plotopts);
ylabel(['Classifier ' measure],'fontweight','bold');
ylim([dpmin dpmax]);
text(min(xlim),dpmax*1.1,'d. Mean PL classifier MTF','fontweight','bold');



%% MTF shape.


% ----------------- SPLIT BY SOUND LEVEL ---------------

mtftypes = unique(stats_selected_BMFdp_stacked.dprime_MTFType);
unittypes = typeorder; %{'ChS','ChT','PBU','PL','PLN','On'};
levels = [30 50 70];

% VS table.
coli = 1; vs_nos = []; dp_nos = []; dp_propLP = [];
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
    
       if sum( dp_nos(coli,:))>=min_n 
            dp_propLP(ui,li) = sum(strcmp(dpmtftype_unittype,'lowpass'))/sum( dp_nos(coli,:));
       else
           dp_propLP(ui,li)  = nan;
       end;
       
       if sum( vs_nos(coli,:))>=min_n 
            vs_propLP(ui,li) = vs_nos(coli,1)/sum( vs_nos(coli,:));
       else
            vs_propLP(ui,li) = nan;
       end;
       
       coli = coli+1;
    end;

    % create a gap between types.
    coli=coli+1;
end;


% ------------------- plot the d' MTF shape proportions. 

subplot('position', [0.3 y2 0.17 h1]); hold on;

for li = 1:size(dp_propLP,1)   
    typei = li;
    thiscolor = typecolormap(typei,:);
    plot(levels+0.5*li-1,dp_propLP(li,:),'+-','color',thiscolor,'linewidth',1); hold on;
end;
box off;
ylabel('Prop. lowpass','fontweight','bold');
text(25,1.05*1.1,'b. Classification MTF','fontweight','bold');
xlabel('Sound level (dB SPL)','fontweight','bold');
xlim([25 75]);
ylim([-.1 1.05]);
set(gca,'xtick',[30 50 70]);

% ---------------- Tests describing the proportions -----------------

% Type by type comparison of MTF shape proporitions across VS vs. d'
fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs VS, for all neuron types:\n');
fprintf('30dB:'); chi2tests.VS_vs_dp.p_30 = chi2test([sum(vs_nos([1:4:end],1:2)); sum(dp_nos([1:4:end],1:2))]);
fprintf('50dB:'); chi2tests.VS_vs_dp.p_50 = chi2test([sum(vs_nos([2:4:end],1:2)); sum(dp_nos([2:4:end],1:2))]);
fprintf('70dB:'); chi2tests.VS_vs_dp.p_70 = chi2test([sum(vs_nos([3:4:end],1:2)); sum(dp_nos([3:4:end],1:2))]);
fprintf('All levels:'); chi2tests.VS_vs_dp.p_allLvls  = chi2test([sum(vs_nos(:,1:2),1); sum(dp_nos(:,1:2),1)]);

% Type by type comparison of MTF shape proporitions across VS vs. d'
for typei = 1:length(typeorder)
    thistype = typeorder{typei};
    typeoffset = 4*(find(strcmp(typeorder,thistype))-1)+1;
    fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs VS, for %s:\n',thistype);
    fprintf('30dB:'); chi2tests.VS_vs_dp.(['p_' thistype '_30']) = chi2test([vs_nos(typeoffset,1:2); dp_nos(typeoffset,1:2)]);
    fprintf('50dB:'); chi2tests.VS_vs_dp.(['p_' thistype '_50']) = chi2test([vs_nos(typeoffset+1,1:2); dp_nos(typeoffset+1,1:2)]);
    fprintf('70dB:'); chi2tests.VS_vs_dp.(['p_' thistype '_70']) = chi2test([vs_nos(typeoffset+2,1:2); dp_nos(typeoffset+2,1:2)]);
    fprintf('All levels:'); chi2tests.VS_vs_dp.(['p_' thistype '_allLvls'])  = chi2test([sum(vs_nos(typeoffset+(0:2),1:2),1); sum(dp_nos(typeoffset+(0:2),1:2),1)]);
end;

% Test describing the effect of level on the proportion of MTF types - c'. 
% Following https://uk.mathworks.com/matlabcentral/answers/606106-chi-square-test-non-parametric
fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs level for all neuron types:\n');
% There are a small proporton of "allpass" and "highpass" MTFs - excluding
% thgem does not change anything. 
%inds = stats_selected_BMFdp_stacked.dprime_MTFType=="lowpass" | stats_selected_BMFdp_stacked.dprime_MTFType == "bandpass";
%[cp.Table,cp.CHI2,cp.P,cp.LABELS] = crosstab(stats_selected_BMFdp_stacked.modLvl_rationalised(inds),stats_selected_BMFdp_stacked.dprime_MTFType(inds));
[cp.Table,cp.CHI2,cp.P,cp.LABELS] = crosstab(stats_selected_BMFdp_stacked.modLvl_rationalised,stats_selected_BMFdp_stacked.dprime_MTFType);
fprintf('\tProportion of low pass:\n'); disp([30,50,70]); disp((cp.Table(:,2)./(sum(cp.Table'))')');
fprintf('\tp=%g, chi^2 =%g, N=%d\n',cp.P, cp.CHI2, sum(cp.Table(:)));

% Test describing the effect of level on the proportion of MTF types - VS
fprintf('\nChi-squared  test of proportions MTF-shape: VS vs level for all neuron types:\n');
% There are a small proporton of "allpass" and "highpass" MTFs - excluding
[vs.Table,vs.CHI2,vs.P,vs.LABELS] = crosstab(stats_selected_BMF_stacked.modLvl_rationalised,stats_selected_BMF_stacked.VS_MTFType);
fprintf('\tProportion of low pass:\n'); disp([30,50,70]); disp((vs.Table(:,2)./(sum(vs.Table'))')');
fprintf('\tp=%g, chi^2 =%g, N=%d\n',vs.P, vs.CHI2, sum(vs.Table(:)));


%% Peak DP against level.

subplot('position', [0.07 y2 0.17 h1]);
%roworder =  [1 2 4 5 6 3];

% D-PRIME 
rowandcol =  {'rationalisedType','modLvl_rationalised'};
% stats.dPrime_reBMFdp = makePopnTable(stats_selected_BMFdp_stacked,rowandcol{:},'dprimes',coreoptions{:} ...
%     ,options{1:2},'fn','Rayleigh>13.8','roworder',roworder);
stats.dPrime_reBMFdp = makePopnTable(stats_selected_BMFdp_stacked,rowandcol{:},'dprimes',coreoptions{:} ...
    ,'fn','Rayleigh>13.8','roworder',roworder);
            
% Plot
for ri = roworder
    errorbar([30 50 70]+ri*.5-1.5,stats.dPrime_reBMFdp.mean(ri,:), ...
        stats.dPrime_reBMFdp.sd(ri,:)./sqrt(stats.dPrime_reBMFdp.n(ri,:)), ...
        'color',typecolorlist(ri,:),'linewidth',1); 
    hold on;
end;

% plotTable2(stats.dPrime_reBMFdp,'sem','colorlist',{'m','g','r','k','b','c'});
ylim([dpmin dpmax]); xlim([25 75]); box off;
xlabel('Sound level (dB SPL)','fontweight','bold');
ylabel(['Classifier ' measure ' at BMF (Hz)'],'fontweight','bold');
text(25,dpmax*1.1,['a. Maximum ' measure ' (at BMF)'],'fontweight','bold');
set(gca,'xtick',[30 50 70]);


%%%%%%%%%%%%%%%%%       VS          %%%%%%%%%%%%%%%%%%


%% ChS neuron type and looking across level instead.
%options = {'fn','depthMod==1','fn','strcmp(rationalisedType,''ChS'')'};   
options = {'fn','strcmp(rationalisedType,''ChS'')'};   
plotopts([3:6 7:18]) = {  'rowtitle','', ...
                      'colorlist',chScolorlist, ...
                      'roworder',  [1:3], ...
                      'axeslist', {'VS'}, ...
                      'axes_y',[y1 h1], ...
                      'axes_x',[0.13 0.18], ...
                      'offset_x',   0.54, ...
                      'xlim',[-750 1250]};        
po2 = plot_VS_dp_SPP4(stats_stacked,'modLvl_rationalised',options,coreoptions,plotopts);
ylabel('Mean VS','fontweight','bold');
text(-750,0.99,['g. Mean ChS MTF-VS'],'fontweight','bold');
xlim([-750 1250]);
text(25,dpmax*1.1,['a. Maximum ' measure ' (at BMF)'],'fontweight','bold');


%% PL neuron type and looking across level.
%options = {'fn','depthMod==1','fn','strcmp(rationalisedType,''PL'')'};     
options = {'fn','strcmp(rationalisedType,''PL'')'};     
plotopts([3:4 5 6 15:18]) = {  'rowtitle','', ...
                           'colorlist',PLcolorlist, ...
                           'offset_x', 0.78, ...
                            'xlim',[-750 1250]};               
po3 = plot_VS_dp_SPP4(stats_stacked,'modLvl_rationalised',options,coreoptions,plotopts);
ylabel('Mean VS','fontweight','bold');
xlim([-1000 1500]);
text(-1000,0.99,['g. Mean PL MTF-VS'],'fontweight','bold');


%% MTF shape 


maxv = max([sum(vs_nos,2);  sum(vs_nos,2)]);
subplot('position', [0.3 y1 0.17 h1]); hold on;

for li = 1:size(vs_propLP,1)   
    typei = li;;
    thiscolor = typecolormap(typei,:);
    plot(levels+0.5*li-1,vs_propLP(li,:),'+-','color',thiscolor,'linewidth',1); hold on;
end;
box off;
ylabel('Prop. lowpass','fontweight','bold');
%set(gca,'xticklabel',unittypes,'xtick',[2 6 10 14 18 22]);
text(25,1.05*1.1,'f. MTF-VS classification','fontweight','bold');
xlabel('Sound level (dB SPL)','fontweight','bold');
xlim([25 75]);
ylim([-.1 1.05]);
set(gca,'xtick',[30 50 70]);




%% Peak VS against level.
subplot('position', [0.07 y1 0.17 h1]);

% VS 
rowandcol =  {'rationalisedType','modLvl_rationalised'};
% stats.VS_reBMF = makePopnTable(stats_selected_BMF_stacked,rowandcol{:},'VS',coreoptions{:}, ...
%     options{1:2},'fn','Rayleigh>13.8','roworder',roworder);
stats.VS_reBMF = makePopnTable(stats_selected_BMF_stacked,rowandcol{:},'VS',coreoptions{:}, ...
   'fn','Rayleigh>13.8','roworder',roworder);

% Plot
for ri = roworder
    errorbar([30 50 70]+ri*.5-1.5,stats.VS_reBMF.mean(ri,:), ...
        stats.VS_reBMF.sd(ri,:)./sqrt(stats.VS_reBMF.n(ri,:)), ...
        'color',typecolorlist(ri,:),'linewidth',1); 
    hold on;
end;

%plotTable2(stats.VS_reBMF,'sem','colorlist',{'m','g','r','k','b','c'});
ylim([0 1]); xlim([25 75]); box off;
xlabel('Sound level (dB SPL)','fontweight','bold');
ylabel('VS at BMF (Hz)','fontweight','bold');

text(25,1*1.1,'e. Maximum VS (at BMF)','fontweight','bold');
set(gca,'xtick',[30 50 70]);

% -------- Plot the legend axis -----------

subplot('position', [0.35 0.48 0.3 0.05]);

r1y = 0.49; r2y = 0.51;
xvals = [0 0.1 0.2];
arrayfun(@(x,c) line([x x+0.03],[r1y r1y],'color',typecolormap(c,:),'linewidth',2),xvals,[1 3 5]);
arrayfun(@(x,t) text(x+0.05,r1y,typeorder{t}),xvals,[1 3 5]);

arrayfun(@(x,c) line([x x+0.03],[r2y r2y],'color',typecolormap(c,:),'linewidth',2),xvals,[2 4 6]);
arrayfun(@(x,t) text(x+0.05,r2y,typeorder{t}),xvals,[2 4 6]);

axis off;




function [p chistat] = chi2test(vs_nos)

% Get the numbers of observed of each type.
nLowpass = vs_nos(:,1);
nBandpass = vs_nos(:,2);
nEachType = nLowpass + nBandpass;

% pooled proportion of lowpass
p0 = sum(nLowpass)./sum( nEachType );
% expected numbers for each unit type.
nLowpass_expected = p0 * nEachType;
nBandpass_expected = (1-p0) * nEachType;

% Compile all expected and observed numbers.
observed = [nLowpass; nBandpass];
expected = [nLowpass_expected; nBandpass_expected];

% calculate the chi2 statistic.
chistat = sum( (observed - expected).^2 ./expected );

% Calculate the p value for 1 degree of freedom (2 MTF types).
p = 1- chi2cdf(chistat,1);

% Report this. 
fprintf('Chi-squared statistic: %g, p-value = %g, N= %d\n', chistat, p,min(nEachType) );
% N.B. Given that this is comparing two different meausures from each data
% set we are quoting the number of minimum number of datasets from each group. 

