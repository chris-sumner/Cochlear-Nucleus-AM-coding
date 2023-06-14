function allSpikeTimes88and91 = loadSpikeTimes
% loadSpikeTimes load and combine spike times from small(er) files.
%
% The whole dataset was too large to put on github.

load('rawdata\justSpikeTimesA.mat');
load('rawdata\justSpikeTimesB.mat');
load('rawdata\justSpikeTimesC.mat');
load('rawdata\justSpikeTimesD.mat');

allSpikeTimes88and91 = [spikeTimesA; spikeTimesB; spikeTimesC; spikeTimesD]