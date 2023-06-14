function dataSummary(data, EXPLOGLIST, variable_list)
% dataSummary - summarise the ranges of values of each column in wholeList
%
% dataSummary(data, EXPLOGLIST, [variable_list])
%   data    - an numeric matrix of Scholes' specification.
%             or a cell array.
%   variable_list - cell array of names from EXPLOGLIST to display.
%                   'all' displays all variables.
%                   if ommitted it displays a default selection.
% This could obviously by modified to be more flexible than it is!

% A default list of basic stimulus parameters and variables. 
if nargin<3  
    variable_list = {'cf','cf_thr','modLevel','modFreq','amToneDur','carrFreq','depthMod'};
end;

% Option to display statistics for all vairables handed in.
if isnumeric(data) &&  isstr(variable_list) && strcmp(variable_list,'all')
    variable_list = EXPLOGLIST;
end;
    
if iscell(data) && isstr(data{1})
   uvals = unique(data);
   for i=1:length(uvals)
       nvals(i) = sum(strcmp(data,uvals{i}));
   end;
   fprintf('Strings:\n');
   disp(uvals');
   fprintf('Number:\n');
   disp(nvals);
   return;
end;

% Start by summarising the range of values of each:
for i=1:length(variable_list)
   fprintf('\n-------------------------------------\n');
   fprintf('Summary: %s\n', variable_list{i});
      
   vals = data(:,strcmp(EXPLOGLIST,variable_list{i}));
   uvals = unique(vals(~isnan(vals)));

   fprintf('NaNs: %d\n',sum(isnan(vals)));
   fprintf('Mean: %g\n',nanmean(vals(~isnan(vals))));
   fprintf('Maximum: %g\n',nanmax(vals(~isnan(vals))));
   fprintf('Minimum: %g\n',nanmin(vals(~isnan(vals))));

   if length(uvals)>20
      figure; hist(vals,20); xlabel(variable_list{i});
   else
      fprintf('Values:\n');
      disp(uvals');

      hvals = hist(vals,uvals);
      fprintf('Number:\n')
      disp(hvals);
   end;
             
end;

