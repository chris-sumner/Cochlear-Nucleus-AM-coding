function plotSoftMaxFit2ConfMat(fit,confmat)
% plotSoftMaxFit2ConfMat Plots the fit from fitSoftMax2ConfMat
%

% Compute the probabilty function of the confusion matrix.
pconfmat = confmat./(ones(fit.m,1)*sum( confmat ));  % Convert confusion matrix to probability.

figure; 
subplot(2,2,1);
cmax = max(pconfmat(:));
image(pconfmat,"CDataMapping","scaled");
set(gca,'CLim',[0 cmax]);
xlabel('Class'); ylabel('Choice');title(sprintf('Data (max: %G)',cmax));
if ~isfield(fit,'pmodel_init')
    fit.pmodel_init = softMaxEval(fit.Z_init,fit.B_init);
end;
subplot(2,2,2); image(real(fit.pmodel_init),"CDataMapping","scaled");
set(gca,'CLim',[0 cmax]);
xlabel('Class');title(sprintf('Initial prediction (err: %g)',fit.rmserr_init));
%if ~isfield(fit,'pmodel_final')
    fit.pmodel_final = softMaxEval(fit.Z_final,fit.B_final);
%end;

subplot(2,2,3); image(real(fit.pmodel_final),"CDataMapping","scaled"); 
set(gca,'CLim',[0 cmax]);
xlabel('Class'); ylabel('Choice');title(sprintf('Final prediction (err: %g)',fit.rmserr_final));
subplot(2,2,4); 
plot(fit.Z'); hold on; 
plot(fit.B','--');
plot(fit.Z_init); plot(fit.B_init);
xlabel('Class'); ylabel('Metric (Z/B)');
plot(fit.B_final,'--k','linewidth',1);
plot(fit.Z_final,'k','linewidth',1);
title('Fitted parameters');