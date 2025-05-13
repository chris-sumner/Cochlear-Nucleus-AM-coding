function plotPairWise(varnames,varargin)

xSize = 24; ySize = 22; xLeft = (21-xSize)/2;yTop = (30-xSize)/2;
figh=figure;set(figh,'PaperUnits','centimeters');
set(figh,'paperposition',[xLeft yTop xSize ySize],'units','centimeters','position',[2 4 xSize ySize],'DefaultAxesFontSize',8);%,'DefaultAxesColorOrder',colorOrder
paneldim = length(varargin);

% Plot the individual distribution along the main diagonal. Or the do
for pi = 1:paneldim
    subplot(paneldim,paneldim, paneldim*(pi-1) + pi);    
    hist(varargin{pi},50); box off;

    tx = -(diff(xlim)*(0.25+0.06*length(varnames{pi})) + min(xlim)); 
    ty = diff(ylim)*0.7 + min(ylim);    
    text(tx,ty,varnames{pi},'FontWeight','bold');
      xlabel(varnames{pi});
end;

% Loop through and plot each pairwise correlation.
for yi = 1:paneldim
    for xi = (yi+1):paneldim
        subplot(paneldim,paneldim, paneldim*(yi-1) + xi);    
        plot(varargin{xi},varargin{yi},'.'); box off;
        
        if xi == (yi+1)
            ylabel(varnames{yi});
        end;
            
        if yi==1
            title(varnames{xi},'FontWeight','bold');
        end;
        
        % Correlation coeffient.
        xvals = varargin{xi}; yvals = varargin{yi};
        inds = ~isnan(xvals) & ~isnan(yvals) & ~isinf(xvals) & ~isinf(yvals);
        [R,p] = corrcoef(xvals(inds),yvals(inds));

        tx = diff(xlim)*0.1 + min(xlim); 
        ty = diff(ylim)*0.8 + min(ylim);    

        if p(2)<.001
            text(tx,ty,['R^2=' sprintf('%.2g***',R(2).^2)],'FontWeight','bold');
        elseif p(2)<.01
            text(tx,ty,['R^2=' sprintf('%.2g**',R(2).^2)],'FontWeight','bold');
        elseif p(2)<.05
            text(tx,ty,['R^2=' sprintf('%.2g*',R(2).^2)],'FontWeight','bold');
        else
            text(tx,ty,['R^2=' sprintf('%.2g',R(2).^2)]);
        end;
    end;
end;



