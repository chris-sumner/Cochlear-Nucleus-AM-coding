% --------------- Single panel on current axes =-----------------
% Create axes then call.
% Plots mean Z-score below 1kHz, by neuron type.

line([0 7],[2 2],'linestyle','--','color',[.7 .7 .7]);
bh2 = boxplot(stats_oneperunit_stacked.Z_NHPP_meanBelow1k,stats_oneperunit_stacked.rationalisedType,...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled');
ylabel('Z-score','fontweight','bold'); xlabel('Neuron type');
ylim([0 11]); xlim([0.5 6.5]); box off;
text(min(xlim),max(ylim)*1.1,'b. ISIs vs Poisson+DT','fontweight','bold');



