function fModFactor_fit = fitTheoreticalMaxModel_softmax(stats_stacked_4models,meanfModFn)
% fitTheoreticalMaxModel_softmax does what it says. 
%

% List of all possible modulation frequencies (rationalised steps).
flist = unique(stats_stacked_4models.usedamfreqs_dprime_rationalised);

% Code specific to logistic model
    
    % Some meta parameters for the fit.
    gains = [0:0.1:10 10.2:0.2:20 20.5:0.5:50 51:100 102:2:200];  % Total range of gains to consider. Needs to be more than you need.
    % Here we set the range based on a previous fit - and knowing there are
    % relatively fewer very high gains. Equally we don't need superdense
    % set of very low gains. 
    
    % To speed up the fit we precompute as much as possible. 
    [fModFactor_fit, ClassIndex] = prepTheoreticalMaxFit_softmax(stats_stacked_4models,gains, ...
        meanfModFn - min(meanfModFn));

    % Already in fModFactor_fit now.
    %fModFactor_fit.ClassIndex_stacked = ClassIndex; % Need this for fitting later models. 
    %fModFactor_fit.DataIndex_stacked = DataIndex; % Need this for fitting later models. 
    
    % Now we fit the data, passing in the fModFactor_fit structure with all the
    % precomputed stuff.
    % N.B. Gains are effectively random effects in a model - we don't care at 
    % point of fitting. 
    fprintf('Fitting underlying function');
    fModFactor_fit.options = optimset('MaxFunEvals',50000,'MaxIter',50000);   % 50000 is overkill.
    fModFactor_fit.fitted_fn = fminsearch(@(f) mtfErrFn_softmax(f,fModFactor_fit),fModFactor_fit.initModFn,fModFactor_fit.options);
    
    % Post processing involves recomputing all the gains for the best fit,
    % errors etc.
    fModFactor_fit = postTheoreticalMaxFit_softmax(fModFactor_fit);
    
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
    xlabel('f_mod (Hz)'); ylabel('Predicted score');
    % Measured vs. predicted performance.
    subplot(2,2,4);
    plot([fModFactor_fit.measureV(fModFactor_fit.nonnans)],[fModFactor_fit.predmeasureV(fModFactor_fit.nonnans)],'.');
    xlabel('Score (measured)'); ylabel('Score (pred)');
    % 83% of variance - slightly higher than before! 
    fModFactor_fit.R2 = corrcoef(fModFactor_fit.measureV(fModFactor_fit.nonnans),fModFactor_fit.predmeasureV(fModFactor_fit.nonnans)).^2;
    text(0.05,0.9,['R^2:' num2str(fModFactor_fit.R2(2))]);
