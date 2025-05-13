function fModFactor_fit = fitTheoreticalMaxModel_logistic(stats_stacked_4models,meanfModFn)
% fitTheoreticalMaxModel_logistic does what it says. 
%

% List of all possible modulation frequencies (rationalised steps).
flist = unique(stats_stacked_4models.usedamfreqs_dprime_rationalised);

% Code specific to logistic model

    % We fit the individual MTFs to a mean function. Each MTF is related to the
    % mean by a change in gain. However, this function is related to the 
    % performance measure through a logistic function (link), which provides a bounded
    % scale (probability of being correct = chance -> 1).
    
    % This is the logistic function which maps an underlying variable to
    % proportion correct. 2nd function scales for chance performance.
    % x0 = 1;  % Somewhat arbitrary. Keeps the final function positive.
    % p0 = (pchance*(1+exp(x0))-1)/exp(x0)
    % logitfn = @(g) 1./(1 + exp(-(g*f_fn-x0)));
    % pcfn = @(p0,lt) p0 +(1-p0).*lt;
    %
    % Function to invert the logistic function and compute the underlying metric from
    % pecentage correct. It is not as good as fitting, because it cannot
    % account fully for the compression. 
    % c_metric_fn = @(p,p0,x0)  -log( (1 -p0)./(p-p0)  - 1 ) + x0;
    %
    % The approach is make a table with lots of gains, reducing some of the fitting
    % process to table lookup. 
    % gains = [0:0.1:10];
    % logittable = logitfn(gains');  % Generic table of functions for all neurons 
    % pctable = pcfn(p0,logittable); % This is effectively a lookup table for a neuron with a specific chance level. 
    % 
    % % They look like this. 
    % figure;
    % plot(f,pctable')
    
    % Some meta parameters for the fit.
    gains = [0:0.1:10];  % Total range of gains to consider.
    x0 = 1;              % Somewhat arbitrary. Keeps the final function positive.
                         % Means the weightings can be interpreted as "gains".
    
    % To speed up the fit we precompute as much as possible. 
    [fModFactor_fit, p0] = prepTheoreticalMaxFit_logistic(stats_stacked_4models,gains,x0, ...
        meanfModFn - min(meanfModFn));
    fModFactor_fit.p0_stacked = p0; % Need this for fitting later models. 
    
    % Now we fit the data, passing in the fModFactor_fit structure with all the
    % precomputed stuff.
    % N.B. Gains are effectively random effects in a model - we don't care at 
    % point of fitting. 
    fModFactor_fit.options = optimset('MaxFunEvals',50000,'MaxIter',50000);   % 50000 is overkill.
    fModFactor_fit.fitted_fn = fminsearch(@(f) mtfErrFn(f,fModFactor_fit),fModFactor_fit.initModFn,fModFactor_fit.options);
    
    % Post processing involves recomputing all the gains for the best fit,
    % errors etc.
    fModFactor_fit = postTheoreticalMaxFit(fModFactor_fit);
    
    % --------- Plot the results of the fitting ----------
    figure; 
    % This is the resulting function. 
    subplot(2,2,1);
    plot(fModFactor_fit.flist,fModFactor_fit.normalised_fitted_fn);
    xlabel('f_mod (Hz)'); ylabel('Normalised metric values');
    % The range of fitted "gains"
    subplot(2,2,2);
    hist(fModFactor_fit.fitted_normalised_gains,50);
    xlabel('Fitted gains'); ylabel('Number of records');
    % The shapes of the predicted classifier performance. 
    subplot(2,2,3); hold on;
    cellfun(@(f,y) plot(f,y),fModFactor_fit.amfreqs,fModFactor_fit.fitted_measures );
    xlabel('f_mod (Hz)'); ylabel('Predicted F1-score');
    % Measured vs. predicted performance.
    subplot(2,2,4);
    plot([fModFactor_fit.measureV(fModFactor_fit.nonnans)],[fModFactor_fit.predmeasureV(fModFactor_fit.nonnans)],'.');
    xlabel('Score (measured)'); ylabel('Score (pred)');
    % 80% of variance - slightly higher than before! 
    fModFactor_fit.R2 = corrcoef(fModFactor_fit.measureV(fModFactor_fit.nonnans),fModFactor_fit.predmeasureV(fModFactor_fit.nonnans)).^2;
    text(0.05,0.9,['R^2:' num2str(fModFactor_fit.R2(2))]);
