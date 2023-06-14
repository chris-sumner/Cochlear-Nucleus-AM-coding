function [dh lh R] = plotVerticalSideBars(stackdata,xfield,groupfield,stackinds)


groupstyles = {'m','d','g','c','k','s','b'};
groups = unique(stackdata.(groupfield));
groups = groups(cellfun(@(x) ~isempty(x),groups));

typeposition = [5 4 1 2 6 3]; % Which position to put each group.
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
    mean_dp = nanmean(stackdata.(xfield)(inds));    
    sd_dp = nanstd(stackdata.(xfield)(inds));
    sem_dp = sd_dp/sqrt(length(inds)); 
    ci95_dp = 1.96*sem_dp;
    
    % Where to plot this group.
    dpbar_xval = meandpbar_xvals(typeposition(i));
    
     patch([dpbar_xval-dpbar_wid dpbar_xval-dpbar_wid dpbar_xval+dpbar_wid dpbar_xval+dpbar_wid], ...
        [mean_dp-sd_dp mean_dp+sd_dp mean_dp+sd_dp mean_dp-sd_dp], ...
        plotClrs(i),'edgealpha',0,'facealpha',0.1); 
   patch([dpbar_xval-dpbar_wid dpbar_xval-dpbar_wid dpbar_xval+dpbar_wid dpbar_xval+dpbar_wid], ...
        [mean_dp-ci95_dp mean_dp+ci95_dp mean_dp+ci95_dp mean_dp-ci95_dp], ...
        plotClrs(i),'edgealpha',0,'facealpha',0.3); 
%     line([dpbar_xval dpbar_xval], [mean_dp-sd_dp mean_dp+sd_dp], ...
%         'color',plotClrs(i),'linewidth',2); 
    line([dpbar_xval-0.2 dpbar_xval+0.2], [mean_dp mean_dp], ...
        'color',plotClrs(i),'linewidth',2);            
    
    hold on;
   
end;

set(gca,'XColor',[0.99 0.99 0.99 ],'xtick',[]);
% axis([0 9 -2.5 0]); box off;
% xlabel('Z_I_S_I','fontweight','bold'); 




% 
% 
% % Plot the mean and s.d.s of the d prime and Z as bars outside the main
% % axes.
% axes('position',[x1_Zd+x2_Zd y1_Zd x3_Zd y2_Zd-xy_spc]);hold on;
% amfreq = 150; meanzbar_yvals = [-1:-1:-6]/3;
% zbar_hgt = 0.15; dpbar_wid = 0.15;
% meandpbar_xvals = [-1:-1:-6]/2;
% typeposition = [5 4 1 2 6 3];
% for i=1:length(typelist)
%     % Select the unit data 
%     clear tmp;
%     unittype= typelist{i};
%     if strcmp(typelist{i},'On')
%         unitinds = find(cellfun(@(x) ~isempty(strfind(x,unittype)),datavalues.type) & ...
%             testCriteria);
%     else
%         unitinds = find(strcmp(datavalues.type,unittype) & testCriteria);
%     end;
%     
%     datainds = find([datavalues.usedmodfs{unitinds}]==amfreq);
%     tmp.dprime =  [datavalues.dprimes{unitinds}];    
%     tmp.tmp = cellfun(@(x) x.mean_Z, datavalues.Zstats(unitinds),'Unif',false );
%     tmp.mean_Z = [tmp.tmp{:}];    
%     
%     mean_Z = nanmean(tmp.mean_Z(datainds));
%     sd_Z = nanstd(tmp.mean_Z(datainds));
%     sem_Z = sd_Z/sqrt(length(datainds)-1); 
%     ci95_Z = 1.96*sem_Z;
%     zbar_yval = meanzbar_yvals(typeposition(i));
%     
%     patch([mean_Z-sd_Z mean_Z+sd_Z mean_Z+sd_Z mean_Z-sd_Z ], ...
%         [zbar_yval-zbar_hgt zbar_yval-zbar_hgt zbar_yval+zbar_hgt zbar_yval+zbar_hgt], ...
%         plotClrs(i),'edgealpha',0,'facealpha',0.1); 
%     patch([mean_Z-ci95_Z mean_Z+ci95_Z mean_Z+ci95_Z mean_Z-ci95_Z ], ...
%         [zbar_yval-zbar_hgt zbar_yval-zbar_hgt zbar_yval+zbar_hgt zbar_yval+zbar_hgt], ...
%         plotClrs(i),'edgealpha',0,'facealpha',0.3); 
%     line([mean_Z mean_Z],[zbar_yval-0.2 zbar_yval+0.2], ...
%         'color',plotClrs(i),'linewidth',2); 
% %     line([mean_Z-sd_Z mean_Z+sd_Z],[zbar_yval zbar_yval], ...
% %         'color',plotClrs(i),'linewidth',2,'alpha',0.5); 
% end;
% set(gca,'YColor',[0.99 0.99 0.99 ],'ytick',[]);
% axis([0 9 -2.5 0]); box off;
% xlabel('Z_I_S_I','fontweight','bold'); 