function [dh lh R] = plotHorizontalSideBars(stackdata,xfield,groupfield,stackinds)

typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;
groups = typeorder;

typeposition = [5 4 1 2 3 6]; % Which position to put each group.
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
    
    QZ = quantile(stackdata.(xfield)(inds),[.025 .25 .50 .75 .975]);
    mean_Z = QZ(3);
    sd_Z = QZ([2,4]);
    ci95_Z = QZ([1,5]);
    
    % Where to plot this group.
    zbar_yval = meanzbar_yvals(typeposition(i));

    
    patch([ci95_Z fliplr(ci95_Z) ], ...
        [zbar_yval-zbar_hgt zbar_yval-zbar_hgt zbar_yval+zbar_hgt zbar_yval+zbar_hgt], ...
        plotClrs(i,:),'edgealpha',0,'facealpha',0.1); 
    patch([sd_Z fliplr(sd_Z) ], ...
        [zbar_yval-zbar_hgt zbar_yval-zbar_hgt zbar_yval+zbar_hgt zbar_yval+zbar_hgt], ...
        plotClrs(i,:),'edgealpha',0,'facealpha',0.5); 
    line([mean_Z mean_Z],[zbar_yval-0.2 zbar_yval+0.2], ...
        'color',plotClrs(i,:),'linewidth',2); 
    
    hold on;
   
end;

set(gca,'YColor',[0.99 0.99 0.99 ],'ytick',[]);


