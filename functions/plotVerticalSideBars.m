function [dh lh R] = plotVerticalSideBars(stackdata,xfield,groupfield,stackinds)

typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;
groups = typeorder;

typeposition = [5 4 1 2 3 6]; % Which position to put each group.

%plotClrs = 'mgcrkb';
plotClrs = typecolorlist;
meanzbar_yvals = [-1:-1:-6]/3;      
zbar_hgt = 0.15; dpbar_wid = 0.15;
meandpbar_xvals = [-1:-1:-6]/2;

for i=1:length(groups)

    % Flag for all elements in a group.
    if iscell(groups) && isstr(groups{i})
        gflag = strcmp(stackdata.(groupfield),groups{i});
    elseif isnumeric(groups)
        gflag = stackdata.(groupfield)==groups(i);
    end;
    % Find those units to include.
    inds = find(gflag & stackinds);  
    fprintf('%d of %s\n',length(inds),groups{i});
    
    % Compute statistics.    
    
    % Get the values. Remove nans. 
    Zvals = stackdata.(xfield)(inds);
    Zvals = Zvals(~isnan(Zvals));
    
    Qdp = quantile(stackdata.(xfield)(inds),[.025 .25 .50 .75 .975]);
    mean_dp = Qdp(3);
    sd_dp = Qdp([2,4]);
    ci95_dp = Qdp([1,5]);

    
    % Where to plot this group.
    dpbar_xval = meandpbar_xvals(typeposition(i));

     patch([dpbar_xval-dpbar_wid dpbar_xval-dpbar_wid dpbar_xval+dpbar_wid dpbar_xval+dpbar_wid], ...
        [sd_dp fliplr(sd_dp)], ...
        plotClrs(i,:),'edgealpha',0,'facealpha',0.5); 
   patch([dpbar_xval-dpbar_wid dpbar_xval-dpbar_wid dpbar_xval+dpbar_wid dpbar_xval+dpbar_wid], ...
        [ci95_dp fliplr(ci95_dp)], ...
        plotClrs(i,:),'edgealpha',0,'facealpha',0.1); 
    
    line([dpbar_xval-0.2 dpbar_xval+0.2], [mean_dp mean_dp], ...
        'color',plotClrs(i,:),'linewidth',2);            
    
    hold on;
   
end;

set(gca,'XColor',[0.99 0.99 0.99 ],'xtick',[]);


