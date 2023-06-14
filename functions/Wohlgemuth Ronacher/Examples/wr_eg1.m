% WR_EG1 Example plotting discrete spikes.
% 
% Author:   Robert Mill <rob.mill.uk@gmail.com>
%           MRC IHR
% Date:     14 12 12
% Part of:  (W)ohlgemuth (R)onacher spike distance metric
% See also WR_EG2.


% Set figure scale
scale = 1.5;

% Set spike times
S = [1, 2, 4, 7, 8];
T = [2, 3, 6, 7];

% Create a figure (for spikes)
figure('Units', 'inches', 'Position', [1 1 scale*[3 2]]);

% Add axes
hAxes = axes('Units', 'inches', 'Position', scale*[0.4 0.4 2.4 1.4]);
hold on;

% Plot stems
stem(S,  ones(size(S)), 'b', 'LineWidth', scale);
stem(T, -ones(size(T)), 'r', 'LineWidth', scale);

% Configure axes
set(hAxes, ...
    'XLim', [0 10], ...
    'YLim', [-1 1], ...
    'LineWidth', scale, ...
    'FontSize', 8*scale);

% Add labels
xlabel 'Time (s)'
ylabel 'Spikes'

% end of file %
