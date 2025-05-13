% WR_EG_2 Example of filtering spike trains and computing difference/area.
%
% Author:   Robert Mill <rob.mill.uk@gmail.com>
%           MRC IHR
% Date:     14 12 12
% Part of:  (W)ohlgemuth (R)onacher spike distance metric
% See also WR_EG1.


% Set figure scale
scale = 1.5;

% Set alpha
a = 1;

% Set spike times
S = [1, 2, 4, 7, 8];
T = [2, 3, 6, 7];

% Make time steps
dt   = 0.0001;
allT = 20;
ts   = 0:dt:allT;
N    = length(ts);

% Reserve space for spike train
yS = zeros(N, 1);
yT = zeros(N, 1);

% Make spike train
yS(round(S / dt)+1) = 1;
yT(round(T / dt)+1) = 1;

% Filter waveforms
as = [1 -1+a*dt];
as = conv(as, as);
yS = filter(dt, as, yS);
yT = filter(dt, as, yT);

% Create a figure (for filtered waveforms)
figure('Units', 'inches', 'Position', [1 1 scale*[3 2]]);

% Add axes
hAxes = axes('Units', 'inches', 'Position', scale*[0.4 0.4 2.4 1.4]);
hold on;

% Plot blue and pink areas
area(ts,  yS, 'LineWidth', scale, 'FaceColor', [0.5 0.5 0.9]);
area(ts, -yT, 'LineWidth', scale, 'FaceColor', [0.9 0.5 0.5]);

% Configure axes
set(hAxes, ...
    'XLim', [0 allT], ...
    'YLim', [-1 1], ...
    'LineWidth', scale, ...
    'FontSize', 8*scale);

xlabel 'Time (s)'
ylabel 'Filtered Spikes'

% Create a figure (for difference)
figure('Units', 'inches', 'Position', [1 1 scale*[3 2]]);

% Add axes
hAxes = axes('Units', 'inches', 'Position', scale*[0.4 0.4 2.4 1.4]);
hold on;

% Plot green area
area(ts, (yS - yT).^2, 'LineWidth', scale, 'FaceColor', [0.5 0.9 0.5]);

% Configure axes
set(hAxes, ...
    'XLim', [0 allT], ...
    'YLim', [0 0.5], ...
    'LineWidth', scale, ...
    'FontSize', 8*scale);

% Add labels
xlabel 'Time (s)'
ylabel 'Squared Difference'

% end of file %
