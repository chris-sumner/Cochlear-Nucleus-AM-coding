function  npktable = SACPeaks_numByFmodAndType(stats_oneperunit_stacked,xlsname)
% Figure showing the number of SAC peaks, by modulation frequency for each
% neuron type.


figh = figure('position',[100 -100 1000 1000],'paperposition',[.5 .5 20 20]);
amflist = unique(stats_oneperunit_stacked.allamfreqs_rationalised);
amflist = amflist(amflist<1000);
peakedges = [0 2.5 3.5 5.5 7.5 9.5];
types = unique(stats_oneperunit_stacked.rationalisedType);
types = types(2:end);

for ti = 1: length(types)
    numpeakh = zeros(length(amflist),length(peakedges));
    for ami = 1:length(amflist)
        inds = find( stats_oneperunit_stacked.allamfreqs_rationalised==amflist(ami) & ...
                cellfun(@(x) strcmp(x,types{ti}),stats_oneperunit_stacked.rationalisedType) & ...
                ~isnan(stats_oneperunit_stacked.numSACpeaks_p0001));
        
        numpeakh(ami,:) =  hist(stats_oneperunit_stacked.numSACpeaks_p0001(inds), ...
            peakedges);
        numpeakh(ami,:) = numpeakh(ami,:)/sum(numpeakh(ami,:));
        
        % Add to the table of all values for the Excel file. 
        if ti == 1 && ami == 1
            npktable = table(  stats_oneperunit_stacked.unitid(inds)', ...
                               arrayfun(@(x) types{ti},[1:length(inds)]','uni',false), ...
                               stats_oneperunit_stacked.allamfreqs(inds)', ...         
                               stats_oneperunit_stacked.numSACpeaks_p0001(inds)', ...
                              'VariableNames', ...
                              {'NeuronID','NeuronType','AMfrequency_Hz','NumberOfSSACpeaks'});                          
        else
            npktable = vertcat(npktable, ...
                            table(  stats_oneperunit_stacked.unitid(inds)', ...
                              arrayfun(@(x) types{ti},[1:length(inds)]','uni',false), ...
                              stats_oneperunit_stacked.allamfreqs(inds)' , ...         
                              stats_oneperunit_stacked.numSACpeaks_p0001(inds)', ...
                              'VariableNames', ...
                              {'NeuronID', 'NeuronType','AMfrequency_Hz','NumberOfSSACpeaks'}));                          
        end;
        
    end;
    subplot(3,3,ti);
    bh = bar(amflist,numpeakh,'stacked');
    title(types{ti});
    %set(bh,'xtick',[2:2:10]);
    set(bh(1).Parent,'xtick',[150:200:950]);
    
    if ti==1 || ti==4
        ylabel('Proportion of datasets');
    end;
    if ti>=4
        xlabel('Modulation frequency (Hz)');
    end;
    box off;
    ylim([0 1]);
    
end;

npktable.description(npktable.NumberOfSSACpeaks<2.5)  = {'Zero-lag peak or period peaks missing'};
npktable.description(npktable.NumberOfSSACpeaks==3)  = {'Phase-locking'};
npktable.description(npktable.NumberOfSSACpeaks==4)  = {'1-2 peaks at < period-lag '};
npktable.description(npktable.NumberOfSSACpeaks==5)  = {'1-2 peaks at < period-lag '};
npktable.description(npktable.NumberOfSSACpeaks==6)  = {'3-4 additional peaks'};
npktable.description(npktable.NumberOfSSACpeaks==7)  = {'3-4 additional peaks'};
npktable.description(npktable.NumberOfSSACpeaks==8)  = {'5-6 additional peaks'};
npktable.description(npktable.NumberOfSSACpeaks==9)  = {'5-6 additional peaks'};
npktable.description(npktable.NumberOfSSACpeaks>9)  = {'7+ additional peaks'};

legend({'Zero-lag peak or period peaks missing','Phase-locking','1-2 peaks at < period-lag ','3-4 additional peaks', ...
        '5-6 additional peaks','7+ additional peaks'});
    
if nargin>1    
    writetable(   npktable,  xlsname );
end;
