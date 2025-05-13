
%% stats Figure 

% Plot a bar chart showing the R values for the different models. 
Rfn = @(x) ( x.Rsquared.Ordinary );    % Quote R^2 to facilitate discussion in terms of variance explained. 

R150s = [nan nan nan Rfn(nlm.VS_150ish) Rfn(nlm.SACpeaks01Sal_150ish_p0001) Rfn(nlm.SACpeak2345Sal_150ish_p0001)  ...
        Rfn(nlm.Z_150ish)   Rfn(nlm.SACpeak12345Sal_150ish_p0001) ...
        Rfn(nlm.reliability_150ish) ]

R150s_noChS = [  Rfn(nlm.Z_150ish_noChS) Rfn(nlm.SACpeak12345Sal_150ish_p0001_noChS)   Rfn(nlm.reliability_150ish_noChS)  ];    
    
R150_2pars = [Rfn(nlm.TwoVars_PeaksPlusZ_150ish)  ...
              Rfn(nlm.TwoVars_PeaksPlusReliability_150ish)  ...
              Rfn(nlm.TwoVars_ZPlusReliability_150ish)];

R150_3pars = [ Rfn(nlm.AllSensibleVars_minusNuFl_150ish)];


regressionpanel_ymin = Rfn(nlm.stimonly) - 0.05; % Sets the low limit of the vertical axis.

% Plots a load of them.  
hold on;
ylim([regressionpanel_ymin 0.85]); xlim([0 14]);
stimR = [Rfn(nlm.stimonly) nan(1,8)];
nonSpecRs = [nan Rfn(nlm.type) Rfn(nlm.CV) nan(1,6)];
bh = bar([stimR;nonSpecRs;R150s]',2); 
bh_noChS = bar([7.6 8.6 9.6], R150s_noChS,0.3);

bh_2pars = bar([10.5:0.5:11.5], R150_2pars,1);
bh_3pars = bar([13], R150_3pars,1);

set(bh(1),'facecolor','w'); set(bh(2),'facecolor','k'); %set(bh(3),'facecolor',[.7 .7 .7]);
  
line([0 18],[Rfn(nlm.stimonly) Rfn(nlm.stimonly)],'linestyle','--','color',[.5 .5 .5]);
line([0 18],[nlm.stim_theoreticalmax.overallR nlm.stim_theoreticalmax.overallR].^2,'linestyle','--','color',[.7 .7 .7]);

th = text([1:6 7.25:1:9.25 10.4 11 11.6 13],regressionpanel_ymin*.96*ones(1,13), ...
    {'stimuli','+neuron type','+CV','+VS','+SAC peak 1 (1/f_m)', '+SAC peaks 2-5','+Z_I_S_I', ...
     '+SAC peaks 1-5','+Neural reliability', ...
     '  best statistics ','  of 2 of the 3','+combinations', '+3 best'});
set(th,'rotation',-90)
set(gca,'xticklabel',[]);

ylabel('Variance explained (%)','fontweight','bold');
box off;
%text(0.5,0.88,'b. Predictive power of different models','fontweight','bold');
text(1,Rfn(nlm.stimonly)+0.15,'stimuli ');
text(1,Rfn(nlm.stimonly)+0.12,' only');
line([1.5 1.5],Rfn(nlm.stimonly)+[.1 .01],'color','k');
line([1.4 1.5],Rfn(nlm.stimonly)+[.03 .01],'color','k');
line([1.6 1.5],Rfn(nlm.stimonly)+[.03 .01],'color','k');
text(10, nlm.stim_theoreticalmax.overallR.^2+0.03,'model maximum');
lh = legend('Included in all models','+ neuron statistic','+ f_m_o_d=150Hz statistic','without ChS neurons','+2 statistics','+3 statistics','location','northwest'); 
set(lh,'box','off');
