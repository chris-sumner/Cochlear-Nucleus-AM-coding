%% ------------------------------------------------------------------------

% Temporal resolution of information in synapses?
% This plot shows how the "preferred" temporal resolution of the classifier
% depended on neuron type.

% Generally, neurons were better at small time windows. However, regular 
% firing neurons preferred longer windows. 

bins = [1 2 5 10 20 50];
typeorder = rationalisedTypeNames(2:7);
hcounts = zeros(length(typeorder),length(bins));
for ti = 1:6
    inds = strmatch( typeorder{ti},statsmat_selected.rationalisedType);
    hcounts(ti,:) = hist([statsmat_selected.classifier_tau{inds}],[1 2 5 10 20 50]);

    if ti == 1
        tauTable = table([statsmat_selected.type(inds)]', ...
            [statsmat_selected.classifier_tau{inds}]', ...
            'VariableNames',{'NeuronType','BestClassifierTau_ms'});
    else
        tauTable = vertcat(tauTable, ...
             table([statsmat_selected.type(inds)'], ...
            [statsmat_selected.classifier_tau{inds}]', ...
            'VariableNames',{'NeuronType','BestClassifierTau_ms'}));
    end;
end;
nhcounts = 100*hcounts./(sum(hcounts,2)*ones(1,6));

figh = figure('position',[100 100 500 400],'paperposition',[.5 .5 10 8]);
bh= bar(nhcounts');
lh = legend(typeorder);
set(lh,'box','off');
for i=1:length(bh)
    set(bh(i),'FaceColor',typecolorlist(i,:));
end;
set(gca,'xticklabel',bins);
box off; 
xlabel('Best classifier resolution (tau in ms)','fontweight','bold');
ylabel('Percentage of neuron type','fontweight','bold');
text(0,max(ylim)*1.2,'How classifier resolution varied with neuron type');


clear bh lh bins typeorder hcounts figh inds;

writetable(tauTable,'datavalues_FigS12.xlsx');
%print -dtiff -r300 FigureS8_bestclassifiertaus.tif
