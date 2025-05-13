function [PeasronR] = autoCorrStatsAndPlot(statin,a,n,labelstr,gap)
% autoCorrStatsAndPlot How well a statistic correlates with itself.
%                      across modulation frequency.
%
% statin the statistic to plot. Hand in the field from statsmat. 
% a      index of the starting value e.g. 1 = lowest modulation frequency. 
% n      vector of how many more to plot. 
%        e.g a = 1, n = 1 will plot the first value against the second
%            a = 1, n = 6 will plot values 1-6 against 2-7
% labelstr name of the statistic for title label. 
% 
% Returns the Pearson correlation of all the values.

if nargin<5
    gap = 0;
end;

x_ind = a:a+n-1;
y_ind = a+1+gap:a+n+gap;

% Find the number of frequencies in each dataset. 
n_vals = cellfun(@(x) length(x),statin,'uni',false);

% Find all the valid values for the Pearson correlation.
x_vals = cellfun(@(x,y) x(x_ind(y_ind<y)),statin,n_vals,'uni',false); x_vals = [x_vals{:}]';
y_vals = cellfun(@(x,y) x(y_ind(y_ind<y)),statin,n_vals,'uni',false); y_vals = [y_vals{:}]';
% Compute Pearson correlation across all values. 
notnans = ~isnan(x_vals) &  ~isnan(y_vals);
R = corr(x_vals(notnans),y_vals(notnans));

% For axes limits. 
maxvals = nanmax([x_vals; y_vals]);
minvals = nanmin([x_vals; y_vals]);

% Plots all the pairs that exist. 
hold on;
cellfun(@(x,y) plot( x(x_ind(y_ind<y)),x(y_ind(y_ind<y)),'.'),statin,n_vals); 

% Label the axes. 
ylabel([labelstr ' ' num2str(a+1+gap) '-' num2str(a+n+gap)]); 
xlabel([labelstr ' ' num2str(a) '-' num2str(a+n-1)]);
axis([minvals maxvals  minvals maxvals]);
text(minvals+maxvals*0.1,maxvals*0.9,['L:' num2str(n) ', R:' num2str(R)]);

title([labelstr ' vs next ' labelstr ' across f_m_o_d']);

end

