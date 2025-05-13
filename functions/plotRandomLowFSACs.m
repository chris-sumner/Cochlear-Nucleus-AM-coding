% The aim here is a scattergun check/browse of the SACs and the statistics
% extracted from them.
% We are particularly interested in whether the statistics were extracted
% from the low modulation frequencies, so we focus there. 



figure;
ah = subplot('position',[.15 .2 .7 .7]); 

rSACplot = [];
rSACplot.ModF = [];

for saci=1:30
    
    % Generate the list of all units
    rSACplot.unitSetList = ([statsmat.unitid{:}]);
    rSACplot.unitList = unique(rSACplot.unitSetList);

    
    % It's possible that SACs were not run (not enough spikes or no
    % phase-locking).
    while isempty(rSACplot.ModF)
        % Select one at random
        rSACplot.unitInd = randi(length(rSACplot.unitList)); % Which index for the unit (not dataset);
        rSACplot.unitID = rSACplot.unitList(rSACplot.unitInd);

        % Pick a dataset
        rSACplot.numDataSets = length(find( rSACplot.unitSetList ==rSACplot.unitList(rSACplot.unitInd)));
        rSACplot.dataSetNum = randi(rSACplot.numDataSets);
        rSACplot.Inds = find(rSACplot.unitSetList == rSACplot.unitID);  % All the datasets for that unit
        rSACplot.Ind = rSACplot.Inds(rSACplot.dataSetNum);  % Overall index to that dataset. 

        % Pick a low modulation frequency.
        if ~isempty(unitVSoutputs(rSACplot.Ind).SACpeaks)
            rSACplot.setAMFreqs = [unitVSoutputs(rSACplot.Ind).SACpeaks.amfreq];
            rSACplot.setLowAMFreqs = rSACplot.setAMFreqs(rSACplot.setAMFreqs<400); 
            % If there are no low frequencies
            if length(rSACplot.setLowAMFreqs)<1
                continue;
            end;
            rSACplot.amFreqInd =  randi(length(rSACplot.setLowAMFreqs));
            rSACplot.ModF = rSACplot.setAMFreqs(rSACplot.amFreqInd);

            % A final check that there significant some peaks at that modulation frequency  
            if isnan(statsmat.numSACpeaks_p0001{rSACplot.Ind}(rSACplot.amFreqInd)) || ...
                    statsmat.cf{rSACplot.Ind}<3000                     % Also don't display low frequency neurons.

                rSACplot.ModF = [];
            end;
            
            
        end;

    end;

    % Specify the display details.
    rSACplot.lineStyle = 'r';
    rSACplot.peakStyle = 'k.';
    rSACplot.trofStyle = 'r.'; % Red dot makes the trough marker almost invisible.
    rSACplot.saliencePositions = [1.5:0.1:2];
    rSACplot.salienceLineColour = 'r';

    plotSACwithStats(rSACplot,unitVSoutputs,statsmat);
    pause;
    hold off;
    rSACplot = [];
    rSACplot.ModF = [];
end;

