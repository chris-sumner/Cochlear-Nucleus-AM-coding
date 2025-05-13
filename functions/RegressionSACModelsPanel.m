
%% stats Figure 

% This version of this panel shows a wider range of statistics.
% Too many for the main paper!

% Plot a bar chart showing the R values for the different models. 
Rfn = @(x) ( x.Rsquared.Ordinary );    % Quote R^2 to facilitate discussion in terms of variance explained. 

% Measures for when including envelope fluctuation
R150s = [ nan nan ...                     
        Rfn(nlm.SACpeak0Sal_150ish_p0001 ) Rfn(nlm.SACpeaks01Sal_150ish_p0001 ) Rfn(nlm.SACpeaks012Sal_150ish_p0001 ) ...
        Rfn(nlm.SACpeaks0123Sal_150ish_p0001 )  Rfn(nlm.SACpeaks01234Sal_150ish_p0001 ) Rfn(nlm.SACpeakAllSals_150ish )  ... 
        nan(1,16); nan(1,9) ...       
        Rfn(nlm.SACpeak1Sal_150ish_p0001) Rfn(nlm.SACpeak12Sal_150ish_p0001) Rfn(nlm.SACpeak123Sal_150ish_p0001) ...
        Rfn(nlm.SACpeak1234Sal_150ish_p0001) Rfn(nlm.SACpeak12345Sal_150ish_p0001) ...
        nan(1,10); nan(1,15) ...
        Rfn(nlm.SACpeak0Sal_150ish_p0001 ) Rfn(nlm.SACpeak0Plus2Sal_150ish_p0001) Rfn(nlm.SACpeak0Plus23Sal_150ish_p0001) ...
        Rfn(nlm.SACpeak0Plus234Sal_150ish_p0001) Rfn(nlm.SACpeak0Plus2345Sal_150ish_p0001) nan(1,4); nan(1,21) ...        
        Rfn(nlm.SACpeak23Sal_150ish_p0001) Rfn(nlm.SACpeak234Sal_150ish_p0001) Rfn(nlm.SACpeak2345Sal_150ish_p0001) ...
    ]

% Set up axes
hold on;
ylim([0.4 0.85]); xlim([0 25]);
set(bh(1),'facecolor','w'); set(bh(2),'facecolor','k'); %set(bh(3),'facecolor',[.7 .7 .7]);

% Plots a load of them.  
bh = bar([R150s]',4); 

% Lines for maximum and baseline for most models.  
line([0 26],[Rfn(nlm.stimonly) Rfn(nlm.stimonly)],'linestyle','--','color',[.5 .5 .5]);
line([0 26],[nlm.stim_theoreticalmax.overallR nlm.stim_theoreticalmax.overallR].^2,'linestyle','--','color',[.7 .7 .7]);

% Finish off the axes
th = text([3:8 10:14 16:20 22:24],.385*ones(1,19), ...
     {'+SAC peak 0 (~CI)', '+0-1','+0-2','+0-3','+0-4','+0-5', ...
    '+1','+1-2','+1-3','+1-4','+1-5', ...
     '+0','+0,2','+0,2-3','+0,2-4','+0,2-5', ...
     '+2-3','+2-4','+2-5' ...
     });
set(th,'rotation',-90)
set(gca,'xticklabel',[]);
ylabel('Variance explained (%)','fontweight','bold');
box off;
text(0.5,.55,'stimuli ');
text(0.3,.52,' only');
line([1.5 1.5],[.48 .43]+0.02,'color','k');
line([1.4 1.5],[.45 .43]+0.02,'color','k');
line([1.6 1.5],[.45 .43]+0.02,'color','k');
text(18, 0.79,'model maximum');
lh = legend('incrementally adding all SAC peaks f_m_o_d=150Hz','incrementally adding peaks 1-5', ...
            'incrementally adding all peaks 0 & 2-5','incrementally adding all peaks 2-5','location','northwest'); 
set(lh,'box','off');
