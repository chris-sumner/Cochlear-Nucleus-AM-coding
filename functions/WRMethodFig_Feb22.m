% Method figure to illustrate the function of the WR classifier.

amf1 = 75;
amf2 = 125;

fontsize1 = 14;
fontsize2 = 10;

% Number of spike distances to show
n_dists = 100;

newfig = true;
if newfig
    % drawing the figure
    g=figure;set(g,'PaperUnits','centimeters');
    xSize = 21; ySize = 11;
    % centre on A4 paper
    xLeft = (21-xSize)/2;yTop = (30-ySize)/2;
    set(g,'paperposition',[xLeft yTop xSize ySize]);
    set(g,'units','centimeters')
    set(g,'position',[1 0 xSize ySize]);
end;


% --------------- Plot the traces of spike distance -------------

% amfreq_inds = ismember(allModFreqs,unitoutputs(this_eg_ind).wroutput(bestInd).stimulusconditions);
amfreq_inds = find(allModFreqs<=2000 & ~isnan(scountVals));

a1 = subplot('position',[0.1 0.6 0.3 0.2]); axis off;
WREgFig(thesespikes(amfreq_inds),amf1,amf1,allModFreqs(amfreq_inds),1,a1);
text(120,0.15,'\Sigma','fontsize',fontsize1,'fontweight','bold')

a2 = subplot('position',[0.1 0.3 0.3 0.2]); axis off;
WREgFig(thesespikes(amfreq_inds),amf1,amf2,allModFreqs(amfreq_inds),1,a2);
text(120,0.6,'\Sigma','fontsize',fontsize1,'fontweight','bold')

% -------------- Left hand text ------------------

atxt = subplot('position',[0 0 0.1 1]); axis off;
texty = 0.35
text(0.15,texty+0.1,[num2str(amf1) 'Hz'],'fontweight','bold','fontsize',fontsize1)
text(0.15,texty+0.05,['vs'],'fontweight','bold','fontsize',fontsize1)
text(0.15,texty,[num2str(amf2) 'Hz'],'fontweight','bold','fontsize',fontsize1)

texty = 0.65
text(0.15,texty+0.1,[num2str(amf1) 'Hz'],'fontweight','bold','fontsize',fontsize1)
text(0.15,texty+0.05,['vs'],'fontweight','bold','fontsize',fontsize1)
text(0.15,texty,[num2str(amf1) 'Hz'],'fontweight','bold','fontsize',fontsize1)

% --------------- Plot distance matrix -----------

a4 = subplot('position',[0.5 0.1 0.3 0.7]); axis off;
imagesc(wr_rerun.wr_spktrdist(1:n_dists,1:n_dists))
ch = colorbar(a4,'SouthOutside');
text(n_dists*.025, n_dists*1.4,'distance between spike trains (arb.)','fontweight','bold', ...
 'fontsize',fontsize2);
box off; axis off;
hold on; 

% Draw boxes and lines linking the traces.
sq1x = 10*(amf1/25)-10+10*[0 0 1 1]; sq1y = 10*(amf1/25)-10+10*[0 1 1 0];
patch(sq1x,sq1y,'w','facealpha',0.1,'edgecolor','k','linewidth',2)

l1x = 10*(amf1/25)-10+10*[-5 0];l1y = 10*(amf1/25)-10+10*[0 0];
line(l1x,l1y,'color','k','clipping','off','linewidth',2);

sq1x = 10*(amf1/25)-10+10*[0 0 1 1]; sq1y = 10*(amf2/25)-10+10*[0 1 1 0];
patch(sq1x,sq1y,'w','facealpha',0.1,'edgecolor','k','linewidth',2)

l2x = 10*(amf1/25)-10+10*[-5 0];l2y = 10*(amf2/25)-10+10*[4 1];
line(l2x,l2y,'color','k','clipping','off','linewidth',2);

% Axes labels and arrows. 
text(0, n_dists*1.05,'stimulus spike train','fontweight','bold','fontsize',fontsize2);
text(n_dists*-0.05, n_dists*1.02,'comparison   spike train','fontweight','bold', ...
    'fontsize',fontsize2,'rotation',90);
% Stimulus spike train
% Comparison spike train. 

% -------------- illutrate the computation of distances -----
ylist = round([1:10:100] + 9*rand([1 10]));
amf3 = 200;

l1x = 10*(amf3/25)-10+10*[0 0]+3;l1y = [0 ylist(end)];
line(l1x,l1y,'linestyle',':','color',[0.3 0.3 0.3],'clipping','off','linewidth',2);
annotation('arrow',[0.718 0.718],[0.88 0.8])

l1x = 10*(amf3/25)-10+3 + [0 30];
arrayfun(@(x) ...
    line(l1x,[x x],'linestyle',':','color',[0.3 0.3 0.3], ...
    'clipping','off','linewidth',2),ylist);

line([l1x(2) l1x(2)]+15,[1 100],'color',[0.3 0.3 0.3],'linewidth',2,'clipping','off');
line([l1x(2)+10 l1x(2)+15],[1 1],'color',[0.3 0.3 0.3],'linewidth',2,'clipping','off');
line([l1x(2)+10 l1x(2)+15],[100 100],'color',[0.3 0.3 0.3],'linewidth',2,'clipping','off');
annotation('arrow',[0.85 0.88],[0.45 0.45]);

text(n_dists*1.3, n_dists*0.4,'choose','fontweight','bold', ...
    'fontsize',fontsize2);
text(n_dists*1.3, n_dists*0.5,'minimum','fontweight','bold', ...
    'fontsize',fontsize2);
text(n_dists*1.3, n_dists*0.6,'distance','fontweight','bold', ...
    'fontsize',fontsize2);

text(n_dists*0.5, n_dists*-0.33,'pick 1 spike train','fontweight','bold', ...
    'fontsize',fontsize2);
text(n_dists*0.53, n_dists*-0.23,'for f_m_o_d= \omega Hz','fontweight','bold', ...
    'fontsize',fontsize2);

text(n_dists*1.1, n_dists*0.85,'pick 1 train for each f_m_o_d','fontweight','bold', ...
    'fontsize',fontsize2,'rotation',90);

