
% ---- Panel. Make axes before calling-------
% Plots an example SAC for the example ChS and PL units. 

% Specify the unit
% This is the same neuron from Figure 2. 
ChS.unitID = 88299021;
ChS.dataSetNum = 7;     % This neuron was held for a long time and lots of conditions run.                    % This is at 30dB SPL for a lot of am frequencies.     
ChS.ModF = 125;       
ChS.lineStyle = 'r';
ChS.peakStyle = 'k.';
ChS.trofStyle = 'r.'; % Red dot makes the trough marker almost invisible.
ChS.saliencePositions = [1.5:0.1:2];
ChS.salienceLineColour = 'r';

plotSACwithStats(ChS,unitVSoutputs,statsmat);

% This is the same neurons from Figure 3.
PL.unitID = 88340053;
PL.dataSetNum = 2;         
PL.ModF = 150;   
PL.lineStyle = 'b';
PL.peakStyle = 'k.';
PL.trofStyle = 'b.'; % Blue dot makes the trough marker almost invisible.
PL.saliencePositions = [2.2:0.1:2.7];
PL.salienceLineColour = 'b';

plotSACwithStats(PL,unitVSoutputs,statsmat);

