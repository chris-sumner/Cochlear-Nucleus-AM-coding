function [dph vsh] = plotMTFsOneType(stats_oneperunit, type, uniti,oneplot)
% Crude plotting function to allow selewcton of MTFs for display.

if nargin<4
    oneplot = false;
end;
    
if islogical(oneplot) && ~oneplot
    figvs = figure('position',[100 100 800 700],'paperposition',[.5 .5 16 14]);
    figdp = figure('position',[100 100 800 700],'paperposition',[.5 .5 16 14]);
    oneplot = false;
elseif ishandle(oneplot)
    dph = oneplot(1);
    vsh = oneplot(2);
    oneplot = true;
elseif islogical(oneplot) && oneplot
    figh = figure('position',[100 100 800 500],'paperposition',[.5 .5 16 14]);
    dph = subplot(1,2,1);
    vsh = subplot(1,2,2);
end;    

typeColorMap;
typeorder = rationalisedTypeNames(2:7);
typecolorlist = typecolormap;
unittypes = typeorder; %{'ChS','ChT','PBU','PL','PLN','On'};

colorbase = typecolorlist(find(strcmp(typeorder,type)),:);

if nargin<3 || isempty(uniti)
    uniti = find(strcmp(stats_oneperunit.rationalisedType,type));
end;

colormap = min(1,([1:ceil(length(uniti)/2)]/floor(length(uniti)/2))'*colorbase);
colormap = [colormap; colormap];

%[1.2 0.8 0.5]

for ui=1:length(uniti)
    if oneplot
        subplot(dph);
    else
        figure(figdp);
        subplot(8,7,ui);
    end;
    if rem(ui,2)>0
        linestyle = '--';
    else
        linestyle = '-';
    end;
    plot( stats_oneperunit.usedamfreqs_dprime{uniti(ui)}, ...
        stats_oneperunit.dprimes{uniti(ui)}, ...
        'color',colormap(ui,:),'linewidth',1,'linestyle',linestyle); 
    set(gca,'fontsize',8);
    hold on;
    title(num2str(uniti(ui)));
    axis([0 2000 -0.2 8]);
end;
axis([0 2000 -0.2 8]);



for ui=1:length(uniti)
    if oneplot
        subplot(vsh);
    else
        figure(figvs);
        subplot(8,7,ui);
    end;
    set(gca,'fontsize',8);
        if rem(ui,2)>0
        linestyle = '--';
    else
        linestyle = '-';
    end;
    plot( stats_oneperunit.usedamfreqs{uniti(ui)}, ...
        stats_oneperunit.VS{uniti(ui)}, ...
        'color',colormap(ui,:),'linewidth',1,'linestyle',linestyle); 
    hold on;
    title(num2str(uniti(ui)));
    axis([0 2000 -0.2 1.2]);

end;
axis([0 2000 -0.2 1.2]);

