function ahs = plotEgConfMatRow(unitoutput,ylims,prefix)

wroutput = unitoutput.wroutput( unitoutput.bestClassifier.index );
fit = wroutput.softmax_fit;
confmat = wroutput.confusionmatrix;

% Compute the probabilty function of the confusion matrix.
pconfmat = confmat./(ones(fit.m,1)*sum( confmat ));  % Convert confusion matrix to probability.

% Prepare confusion matrix axes labels
maxamf = max(wroutput.stimulusconditions);
minamf = min(wroutput.stimulusconditions);
amfrange = maxamf - minamf;
amstep = amfrange/4;
amflabels = 50*round((minamf + [1 2 3 4]*amstep)/50);
amflabelpos = fit.m*[0.25 0.5 0.75 1];


% plot the confusion matrix for the data.
ahs(1) = subplot('position',[0.25 ylims(1) 0.2 ylims(2)]);
cmax = max(pconfmat(:));
image(pconfmat,"CDataMapping","scaled");
set(gca,'CLim',[0 cmax]);
xlabel('Modulation frequency (Hz)'); ylabel('Classifier Choice (Hz)');
set(gca,'xtick',amflabelpos,'xticklabel',amflabels, ...
    'ytick',amflabelpos,'yticklabel',amflabels);

% Display some information. 
text(-1.1*fit.m,fit.m*0.1,[prefix ' ' unitoutput.unitinfo.TypeName ' neuron'],'fontweight','bold');
text(-1.1*fit.m,fit.m*0.3,sprintf('CF: %3.2g kHz',unitoutput.unitinfo.cf/1e3),'fontweight','bold');
text(-1.1*fit.m,fit.m*0.5,sprintf('RMS error: %2.2g',fit.rmserr_final),'fontweight','bold');

% Plot the confusion matrix for the fitted model.
fit.pmodel_final = softMaxEval(fit.Z_final,fit.B_final);
ahs(2) = subplot('position',[0.5 ylims(1) 0.2 ylims(2)]);
image(real(fit.pmodel_final),"CDataMapping","scaled"); 
set(gca,'CLim',[0 cmax],'xtick',amflabelpos,'xticklabel',amflabels);
xlabel('Modulation frequency (Hz)');

% Prep stuff for the fitted parameters
Z_errmax =  max(fit.Z);
Z_errmin =  min(fit.Z);
B_errmax =  max(fit.B);
B_errmin =  min(fit.B);


% Plot the fitted parameters. 
ahs(3) = subplot('position',[0.75 ylims(1) 0.18 ylims(2)]);
%plot(fit.Z'); hold on; 
plot(wroutput.stimulusconditions,fit.Z_final,'k','linewidth',1);hold on;
plot(wroutput.stimulusconditions,fit.B_final,'--b','linewidth',1);
plot(wroutput.stimulusconditions,wroutput.hitrate*8,'LineWidth',1);
plot(wroutput.stimulusconditions,Z_errmax,'k:'); hold on;
plot(wroutput.stimulusconditions,Z_errmin,'k:'); 
plot(wroutput.stimulusconditions,B_errmin,'b:'); 
plot(wroutput.stimulusconditions,B_errmax,'b:'); 
%plot(fit.B','--');
xlabel('Modulation frequency (Hz)'); ylabel('Parameter value');
plot(wroutput.stimulusconditions,fit.B_final,'--b','linewidth',1);
plot(wroutput.stimulusconditions,fit.Z_final,'k','linewidth',1);

% Add a scale for hit-rate.
if  max(Z_errmax) > 4
    hitticks = [0:0.25:1]; 
else
    hitticks = [0:0.1:1]; 
end;    
hitticky = hitticks*8;
hitticks = hitticks(hitticky < max(Z_errmax));
hitticky = hitticky(hitticky < max(Z_errmax));

line( maxamf*[1.02 1.05]'*ones(1,length(hitticky)), ...
        [1 1]'*hitticky,'color','k');
text(1.1*maxamf*ones(1,length(hitticky)),hitticky,num2str(hitticks'),'FontSize',8);
text(maxamf*1.3,max(Z_errmax)*0.2,'Hit-rate','rotation',90,'FontSize',8);


lh = legend({'Metric C','Bias B','hit-rate'},'Location','northeast');
ylim([min(Z_errmin)-.2,max([Z_errmax wroutput.hitrate*8])*1.3]);
xlim([minamf*0.95, maxamf*1.05]); 
set(lh,'box','off');
set(gca,'xtick',amflabels,'box','off');
