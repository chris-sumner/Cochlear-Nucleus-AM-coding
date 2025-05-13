
%% stats Figure 

% This version of this panel shows a wider range of statistics.
% Too many for the main paper!

% Plot a bar chart showing the R values for the different models. 
Rfn = @(x) ( x.Rsquared.Ordinary );    % Quote R^2 to facilitate discussion in terms of variance explained. 

% Measures for when including envelope fluctuation
exactRs =  [nan  nan    nan     Rfn(nlm.VS)   Rfn(nlm.CI)   Rfn(nlm.numSACpeaks) ... 
             Rfn(nlm.SACpeaks1Sal_p0001) Rfn(nlm.SACpeaks2345Sal_p0001) ...
             Rfn(nlm.Z)   Rfn(nlm.SACpeaks12345Sal_p0001) Rfn(nlm.reliability) ...
             Rfn(nlm.envFluct) ];
         
R150s = [nan nan nan Rfn(nlm.VS_150ish) Rfn(nlm.CI_150ish) Rfn(nlm.numSACpeaks_150ish) Rfn(nlm.SACpeak1Sal_150ish_p0001) Rfn(nlm.SACpeak2345Sal_150ish_p0001)  ...
        Rfn(nlm.Z_150ish)   Rfn(nlm.SACpeak12345Sal_150ish_p0001) ...
        Rfn(nlm.reliability_150ish) Rfn(nlm.envFluct_150ish) ]

R150_2pars = [Rfn(nlm.TwoVars_PeaksPlusZ_150ish)  ...
              Rfn(nlm.TwoVars_PeaksPlusReliability_150ish)  ...
              Rfn(nlm.TwoVars_ZPlusReliability_150ish) ...
              Rfn(nlm.TwoVars_NuFlPlusReliability_150ish) ...
              Rfn(nlm.TwoVars_PeaksPlusNuFl_150ish) ...
              Rfn(nlm.TwoVars_ZPlusNuFl_150ish) ...
              ];

R150_3pars = [ Rfn(nlm.AllSensibleVars_minusNuFl_150ish) ...
    Rfn(nlm.AllSensibleVars_minusZ_150ish) ...
    Rfn(nlm.AllSensibleVars_minusSACPeaks_150ish) ...
    Rfn(nlm.AllSensibleVars_minusReliability_150ish) ... 
    ];
R150_4pars =  Rfn(nlm.AllSensibleVars_150ish);

stimR = [Rfn(nlm.stimonly) nan(1,11)];
nonSpecRs = [nan Rfn(nlm.type) Rfn(nlm.CV) nan(1,9)];

% Set up axes
hold on;
ylim([0.4 0.85]); xlim([0 20]);

% Plots a load of them.  
bh = bar([stimR;nonSpecRs;R150s]',2); 
bh_exact = bar([1:12], exactRs,0.3);
bh_2pars = bar([13.5:0.5:16], R150_2pars,0.8);
bh_3pars = bar([17:.5:18.5], R150_3pars,0.8);
bh_max = bar(19.5, R150_4pars,0.5);
set(bh(1),'facecolor','w'); set(bh(2),'facecolor','k'); %set(bh(3),'facecolor',[.7 .7 .7]);

% Lines for maximum and baseline for most models.  
line([0 22],[Rfn(nlm.stimonly) Rfn(nlm.stimonly)],'linestyle','--','color',[.5 .5 .5]);
line([0 22],[nlm.stim_theoreticalmax.overallR nlm.stim_theoreticalmax.overallR].^2,'linestyle','--','color',[.7 .7 .7]);

% Finish off the axes
th = text([1:8 9.25:1:12.25 14 14.6 15.2 17.6 19.5],.385*ones(1,17), ...
    {'stimuli','+neuron type','+CV','+VS','+CI','+#SAC peaks','+SAC peak 1 (1/f_m)', '+SAC peaks 2-5','+Z_I_S_I', ...
     '+SAC peaks 1-5','+Neural reliability', '+Neural fluctuation', ...
     '  best statistics ','  of 2 of the 4','+combinations', '+3 of the best', '+ the 4 best'});
set(th,'rotation',-90)
set(gca,'xticklabel',[]);
ylabel('Variance explained (%)','fontweight','bold');
box off;
text(1,.55,'stimuli ');
text(1,.52,' only');
line([1.5 1.5],[.48 .43]+0.02,'color','k');
line([1.4 1.5],[.45 .43]+0.02,'color','k');
line([1.6 1.5],[.45 .43]+0.02,'color','k');
text(15, 0.79,'model maximum');
lh = legend('Included in all models','+ neuron statistic','+ f_m_o_d=150Hz statistic', ...
    '+ f_m_o_d specific statistic','+2  f_m_o_d=150Hz statistics','+3  f_m_o_d=150Hz statistics', ...
    '+4 f_m_o_d=150Hz statistics','location','northwest'); 
set(lh,'box','off');
