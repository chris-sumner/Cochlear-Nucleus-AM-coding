% -------------- Plot isolated panel --------
% # of peaks different neuron types at low modn frequencies. 

% Find all the data at 150Hz (most units) or 100Hz (minority of units).
bh3 = boxplot(stats_oneperunit_stacked.numSACpeaks_p0001(stackinds), ...
              stats_oneperunit_stacked.rationalisedType(stackinds),...
     'colors',[typecolorlist],'grouporder',typeorder,'boxstyle','filled', ...
     'outliersize',6);
ylabel('# SAC peaks','fontweight','bold'); xlabel('Neuron type');
ylim([-0.2 13]); xlim([0.5 6.5]); box off;
text(min(xlim),max(ylim)*1.1,'d. Interval structure @ f_m_o_d = 150Hz','fontweight','bold');
line([0 7],[2.9 2.9],'linestyle','--','color',[.7 .7 .7]);

