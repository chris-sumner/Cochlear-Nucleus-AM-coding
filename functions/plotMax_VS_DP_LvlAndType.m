function plotMax_VS_DP_LvlAndType(stats_selected_BMF_stacked,stats_selected_BMFdp_stacked)


figure('position',[100 100 800 800],'paperposition',[.5 .5 16 16]); 
subplot(2,1,1);
boxplot( stats_selected_BMFdp_stacked.dprimes, ...`
    {stats_selected_BMFdp_stacked.rationalisedType stats_selected_BMFdp_stacked.modLvl_rationalised}, ...
    'factorgap',[5 1]);

% Delete some of the messy type labels. 
kids =  get(gca,'children');
delete(kids.Children([19:3:35 21:3:36]));
box off;
xlabel('Neuron type/sound level');
ylabel('d'' at BMFdp');
text(1,7,'A. d'' maximum as a function of sound level and unit type','fontweight','bold');

subplot(2,1,2);
boxplot( stats_selected_BMF_stacked.VS, ...`
    {stats_selected_BMF_stacked.rationalisedType stats_selected_BMF_stacked.modLvl_rationalised}, ...
    'factorgap',[5 1]);

% Delete some of the messy type labels. 
kids =  get(gca,'children');
delete(kids.Children([19:3:35 21:3:36]));
box off;
xlabel('Neuron type/sound level');
ylabel('VS at BMF (Hz)');
text(1,1.1,'B. VS maximum as a function of sound level and unit type','fontweight','bold');


