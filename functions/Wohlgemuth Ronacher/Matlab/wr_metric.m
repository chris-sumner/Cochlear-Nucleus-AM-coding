function D = wr_metric(spks, alpha, norm)
% WR_METRIC Compute Wohlgemuth-Ronacher spike distance metric.
%    WR_METRIC(SPKS, ALPHA) takes a cell array of N spike trains, SPK, and
%    returns the pairwise W-R spike distances in an N x N matrix. ALPHA
%    controls the decay of the filter (i.e., the spikes trains are smoothed
%    more for smaller values of ALPHA). SPKS should be a column, and the
%    spike times that appear in each cell should also be in columns.
%
%    WR_METRIC(SPKS, ALPHA, 'norm') normalises the output of the routine so
%    that, for a spike distance X:
%
%    * if ALPHA is small, sqrt(X) is the difference in spike counts
%    * if ALPHA is large, X is the number of spikes with no coincident
%      spike (obviously is ALPHA is very large no spikes coincide exactly)
%    * if ALPHA is scaled by 1/Q and the spike trains are dilated by a
%      factor of Q, X is unchanged
%
%    Example:
%       
%       % Create two one-second spike trains each containing 100 spikes
%       spkts1_s = 60*rand(100, 1);
%       spkts2_s = 60*rand(100, 1);
%
%       % Measure the distance between them using alpha = 1 (without
%       % normalisation)
%       D = wr_metric({spkts1_s; spkts2_s}, 1);
%
%    Note that the units of the spike times are unimportant, provided that
%    ALPHA is appropriately scaled.
%
%    For more implementation details, see wr_metric_mex.c, wr_metric.c and
%    wr_metric.h. For the original work, see
%       S. Wohlgemuth & B. Ronacher (2007), Auditory Discrimination of
%       Amplitude Modulations Based On Metric Distances of Spike Trains,
%       J. Neurophysiol. 97:3082-3092.
%
% Author:   Robert Mill <rob.mill.uk@gmail.com>
%           MRC IHR
% Date:     13 12 12
% Part of:  (W)ohlgemuth (R)onacher spike distance metric
% See also WR_EG1, WR_EG2.


% Use defaults?
if nargin < 3, norm = ''; end

% Run spike metric routine
D = wr_metric_mex(spks, alpha);

% Normalise output?
if strcmpi(norm, 'norm'), D = D * alpha^3 * 4; end

% end of file %
