function [dh lh R] = plotHorizontalSideBars(stackdata,xfield,groupfield,stackinds)

groupstyles = {'sm','gd','co','r^','k+','.b','w.'};

groupstyles = {'m','d','g','c','k','s','b'};
groups = unique(stackdata.(groupfield));
groups = groups(cellfun(@(x) ~isempty(x),groups))

typeposition = [5 4 1 2 6 3]; % Which position to put each group.
%plotClrs = 'mgkbrc';
plotClrs = 'mgcrkb';
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
    mean_Z = nanmean(stackdata.(xfield)(inds));    
    sd_Z = nanstd(stackdata.(xfield)(inds));
    sem_Z = sd_Z/sqrt(length(inds)); 
    ci95_Z = 1.96*sem_Z;
    
    % Where to plot this group.
    zbar_yval = meanzbar_yvals(typeposition(i));

    
    patch([mean_Z-sd_Z mean_Z+sd_Z mean_Z+sd_Z mean_Z-sd_Z ], ...
        [zbar_yval-zbar_hgt zbar_yval-zbar_hgt zbar_yval+zbar_hgt zbar_yval+zbar_hgt], ...
        plotClrs(i),'edgealpha',0,'facealpha',0.1); 
    patch([mean_Z-ci95_Z mean_Z+ci95_Z mean_Z+ci95_Z mean_Z-ci95_Z ], ...
        [zbar_yval-zbar_hgt zbar_yval-zbar_hgt zbar_yval+zbar_hgt zbar_yval+zbar_hgt], ...
        plotClrs(i),'edgealpha',0,'facealpha',0.3); 
    line([mean_Z mean_Z],[zbar_yval-0.2 zbar_yval+0.2], ...
        'color',plotClrs(i),'linewidth',2); 
    
    hold on;
   
end;

set(gca,'YColor',[0.99 0.99 0.99 ],'ytick',[]);


