function [dh lh R optable] = scatterPlotStack(stackdata,xfield,yfield,groupfield,stackinds,xjit,logyR)

% xjit allows some scatter to be added to the x values to make it clearer. 
if nargin<6
    xjit = 0;
end;

if nargin<7
    logyR = false;
end;

typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;


groupstyles = {'s','d','o','^','+','.','.'};
groups = typeorder;

for i=1:length(groups)

    % Flag for all elements in a group.
    if iscell(groups) && isstr(groups{i})
        gflag = strcmp(stackdata.(groupfield),groups{i});
    elseif isnumeric(groups)
        gflag = stackdata.(groupfield)==groups(i);
    end;
    
    inds = find(gflag & stackinds);  % Find those units to include.
    if xjit~=0
        jitvals = -xjit*.5 + xjit*rand(size(inds));
    else 
        jitvals = 0;
    end;
    
    fprintf('%d of %s\n',length(inds),groups{i});
       
    dh{i} = plot(jitvals + stackdata.(xfield)(inds),stackdata.(yfield)(inds), ...
        groupstyles{i}, 'color',typecolorlist(i,:));
    hold on;

    % Add these to the table.
    if i == 1
         optable = table(stackdata.(groupfield)(inds)',stackdata.(xfield)(inds)',stackdata.(yfield)(inds)');
    else        
        optable = vertcat(optable, ...
            table(stackdata.(groupfield)(inds)',stackdata.(xfield)(inds)',stackdata.(yfield)(inds)'));
    end;
end;

optable = renamevars(optable,optable.Properties.VariableNames,{groupfield,xfield,yfield});

% Correlation coeffient.
xvals = stackdata.(xfield)(stackinds); yvals = stackdata.(yfield)(stackinds);
if logyR
    yvals =  real(log(yvals));
end;   
inds = ~isnan(xvals) & ~isnan(yvals) & ~isinf(xvals) & ~isinf(yvals);
R = corrcoef(xvals(inds),yvals(inds));

box off;
lh = legend(groups);
set(lh,'box','off');
xlabel(xfield,'fontweight','bold');
ylabel(yfield,'fontweight','bold');
