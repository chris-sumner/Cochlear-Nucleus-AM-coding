function fModFactor_fit = fitTheoreticalMaxModel_linear(stats_stacked_4models,meanfModFn)
% fitTheoreticalMaxModel_linear does what it says. 
%

    % REMOVE 1ST FREQUENCY - FIXING GAIN TO 1.

% List of all possible modulation frequencies (rationalised steps).
flist = unique(stats_stacked_4models.usedamfreqs_dprime_rationalised);

    % Some meta parameters for the fit.
    gains = [0:0.1:25];  % Total range of gains to consider. 
    
    % To speed up the fit we precompute as much as possible. 
    [fModFactor_fit] = prepTheoreticalMaxFit_linear(stats_stacked_4models,gains, ...
        meanfModFn);
    
    % Now we fit the data, passing in the fModFactor_fit structure with all the
    % precomputed stuff.
    % N.B. Gains are effectively random effects in a model - we don't care at 
    % point of fitting. 
    fModFactor_fit.options = optimset('MaxFunEvals',50000,'MaxIter',50000);   % 50000 is overkill.
    fModFactor_fit.fitted_fn = fminsearch(@(f) mtfErrFn_linear(f,fModFactor_fit),fModFactor_fit.initModFn,fModFactor_fit.options);
    
    % Post processing involves recomputing all the gains for the best fit,
    % errors etc.
    fModFactor_fit = postTheoreticalMaxFit_linear(fModFactor_fit);
    
    % % --------- Plot the results of the fitting ----------
    % figure; 
    % % This is the resulting function. 
    % subplot(2,2,1);
    % plot(fModFactor_fit.flist,fModFactor_fit.normalised_fitted_fn);
    % xlabel('f_mod (Hz)'); ylabel('Normalised metric values');
    % % The range of fitted "gains"
    % subplot(2,2,2);
    % hist(fModFactor_fit.fitted_normalised_gains,50);
    % xlabel('Fitted gains'); ylabel('Number of records');
    % % The shapes of the predicted classifier performance. 
    % subplot(2,2,3); hold on;
    % cellfun(@(f,y) plot(f,y),fModFactor_fit.amfreqs,fModFactor_fit.fitted_measures );
    % xlabel('f_mod (Hz)'); ylabel('Predicted F1-score');
    % % Measured vs. predicted performance.
    % subplot(2,2,4);
    % plot([fModFactor_fit.measureV(fModFactor_fit.nonnans)],[fModFactor_fit.predmeasureV(fModFactor_fit.nonnans)],'.');
    % xlabel('Score (measured)'); ylabel('Score (pred)');
    % % 80% of variance - slightly higher than before! 
    % fModFactor_fit.R2 = corrcoef(fModFactor_fit.measureV(fModFactor_fit.nonnans),fModFactor_fit.predmeasureV(fModFactor_fit.nonnans)).^2;
    % text(0.05,0.9,['R^2:' num2str(fModFactor_fit.R2(2))]);
