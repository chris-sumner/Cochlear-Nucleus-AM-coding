function [dh lh] = plotTable3(tabstr,varargin)
% Version of plotTable for the paper.

% Plot style parameters
lw = 1;
transp = 0.2;
colorlist = {'m','g','r','k','b','c'};

% Colors
ind = find(strcmp(varargin,'colorlist'));
if ~isempty(ind)
    colorlist = varargin{ind +1};
end;

% Lines of means (to drive legend)
for i=1:length(tabstr.rowValues)
    plot(tabstr.colValues,tabstr.mean(i,:),'color',colorlist(i,:),'linewidth',lw);
    hold on;
end;

if sum(strcmp(varargin,'sem'))>0
    % Patch of standard errors.
    sem_upper = tabstr.mean + tabstr.sd./sqrt(tabstr.n);
    sem_lower = tabstr.mean - tabstr.sd./sqrt(tabstr.n);

    for i=1:length(tabstr.rowValues)
        nonnan = find(~isnan(tabstr.mean(i,:)));
        patch([tabstr.colValues(nonnan) fliplr(tabstr.colValues(nonnan))], ...
            [sem_upper(i,nonnan) fliplr(sem_lower(i,nonnan))],  ...
            colorlist(i,:),'edgecolor',colorlist(i,:),'facealpha',transp,'edgealpha',transp);
    end;

    % Redraw lines on top.
    % Lines of means (to drive legend)
    for i=1:length(tabstr.rowValues)
        dh(i) = plot(tabstr.colValues,tabstr.mean(i,:),'color',colorlist(i,:),'linewidth',lw);
    end;
end;

box off;

lvals = tabstr.rowValues';
if isnumeric(lvals)
    lvals = arrayfun(@(x) sprintf('%d',x),lvals,'uni',false);
end;
    
lh = legend(lvals);
set(lh,'box','off');

colName = tabstr.colName;
colName(colName=='_') = ' ';
xlabel(colName,'fontweight','bold');

valueName = tabstr.valueName;
valueName(valueName=='_') = ' ';
ylabel(valueName,'fontweight','bold');
