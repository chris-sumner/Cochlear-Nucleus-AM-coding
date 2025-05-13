% ----------------- Checking the softmax fits -------------------

% To start: load unitoutputs, and generate statsmat, as in mainAnalyses 
allZ =  real([statsmat.dprimes{:}]);

% How many are finite?
length(allZ)
sum(~isfinite(allZ ))  % All computed values are finite.  

% ---------------- How many did not converge? --------------

inds = find(arrayfun(@(x) isfield(x.wroutput(1),'confusionmatrix') , unitoutputs))
didnotconverge = find(arrayfun(@(x) sum(x.wroutput(x.bestClassifier.index).softmax_fit.fitstoconverge==-1) , unitoutputs(inds))==10)
% All of them converged. 

% -------------- How were errors distributed? ---------------
rmserrors = real(arrayfun(@(x) x.wroutput(x.bestClassifier.index).softmax_fit.rmserr_final , unitoutputs(inds)))

figure; 
x = [0:.001:0.12]
hist(rmserrors,x); hold on;
hist(rmserrors(didnotconverge),x);

cs = 100*cumsum(hist(rmserrors,x))/length(rmserrors)
plot(x,cs);
% line([1 1]*x(findnearest(cs, 95)),[0 200]);
% line([1 1]*x(findnearest(cs, 5)),[0 200]);
nanmean(rmserrors)
nanmedian(rmserrors)
nanstd(rmserrors)
xlim([0 0.12]);
line((quantile(rmserrors,[.05 .5 .95])'*[1 1])',([1;1;1]*[0 200])');

% The quantiles of the RMS errors. 
quantile(rmserrors,[.05 0.25 .5 .75 .95])
% 75% have an error of .02 or less (probability).
% 95% have an error of .056 or less (probability).


% ------- Compute the variance of the confusion matrix accounted for? -----
CM2pCM = @(x) x./(ones(length(x),1)*sum(x));
squash = @(x) x(:);

corr_coefs = arrayfun( @(x) corr( ...
    squash(CM2pCM( x.wroutput(1).confusionmatrix)), ...
    squash(softMaxEval([x.wroutput(1).softmax_fit.Z_final ...
                        x.wroutput(1).softmax_fit.B_final])) ),unitoutputs(inds) );    
rsquared = real(corr_coefs.^2);

figure;
x = [0:0.0025:1];
cs = 100*cumsum(hist(rsquared,x))/length(rsquared);
hist(rsquared,x); hold on;
plot(x,cs,'linewidth',2,'color','r');
line((quantile(rsquared,[.05 .25 .5 .75 .95])'*[1 1])',([1;1;1;1;1]*[0 180])');
xlim([0.5 1.01]);

% Quantiles of variance accounted for.
quantile(rsquared,[.05 .25 .5 .75 .95])
% 75% of fits have an R^2 of >0.93
% 50% have R^2 of >0.97

