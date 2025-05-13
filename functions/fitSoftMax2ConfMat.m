function [Z, B, err, fit] = fitSoftMax2ConfMat(confmat,smoothw,noisew,niter,plims,Zmax)
% fitSoftMax2ConfMat   = fits a softmax m-class logistic model to a confusion matrix
%
% Fits a softmax (m-class logistical) model, with boas to a confusion matrix.
% Works by first analytically computing the Z and B values based on
% simplifying assumptions. These become initial values in fit. 
% Fits can be run multiple times. In practice, 1st fit is usually the best
% to the raw confustion matrix. 
%
% Becuase Z is poorly defined very high proabilities the fits can be very
% sensitive to small changes in high probailities (reflecting the
% oversampling of the confusion matrix). Adding a very small smoothness 
% constraint helps to stabilise the fit. Adding noise top the initial paramters 
% and averaging also helps sometimes. A Bayesian fit would probably be better. 
%
%  [Z, B, err, fit] = fitSoftMax2ConfMat(confmat,smoothw,noisew,niter,plims)
%       confmat: confusion matrix. Each column is #/probability choices for each class. 
%       smoohtw: weight which penalises fit according the sum of the
%             squared 1st differential of the parameters (so imposes a 
%             smoothness constraint between neighbouring classes).  
%       noisew: amplitude of uniform noise to add (units of Z/B) to initial parameters
%             for each fit. 
%       niter: number of times to fit the model.
%       plims: [low high] limits on the maximum and minimum confusion
%               matrix probabilty. This did not help stablise the fit.
%       Z: final metrics for peformance in each class.
%       B: final bias towards choosing each class.
%       err: the final RMS error on the confusion matrix (units: probability choice|class)   
%       fit: full matrix of all the fitting stuff
%
% IMPORTANT: rows are choices, columns are classes. 

% Intialise default parameters
if nargin<2
    smoothw = 0;
end;

if nargin<3
    noisew = 0;
end;

% If nboot is zero it returns the initial parameters (does no fit). 
if nargin<4
    niter = 1;
end;

% Set limits on maximum and minimum probabilities. 
if nargin<5
    pmax = 1;
    pmin = 0;
else 
    pmax = max(plims);
    pmin = min(plims);
end;  

if nargin<6
    Zmax = 20;
end;

% Set this to 1 as I don't think it helps. 
% Increasing the number of iterations is probably better. 
maxfiti = 1;

% Prepare some variables.
fit.m = length(confmat);    % Number of classes
npars = fit.m*2;            % Number of model parameters.
nmat = fit.m^2;

% This places a limit on the maximum and minimum probablity.
% Does does it on choices so we still end up with probability density 
% after the conversion to prob.
nsamplesperclass = sum( confmat );
meansamplesperclass = mean(nsamplesperclass);
maxchoices = pmax*meansamplesperclass;
minchoices = pmin*meansamplesperclass;
confmat = max( min(confmat,maxchoices), minchoices );

pconfmat = confmat./(ones(fit.m,1)*nsamplesperclass);  % Convert confusion matrix to probability.
rootmeanfn = @(x,l) sqrt( x/l );

fprintf('fitSoftMax2ConfMat: found %d classes\n',fit.m);

% --- Calcuate the initial bias function (Bs) ---
% This is an approximation of B derived by assuming class does not matter 
% and Z is zero. 
pchoice_confmat = mean(pconfmat,2);   % Mean probability of making each choice irrrespective of choice. 
pchoice_confmat(pchoice_confmat==0) = 1e-5;

% Difficult to decide what the best way to reference initial parameters is. 
[minp, B0_ind] = min(pchoice_confmat);
%fit.B_init = (log(pchoice_confmat) - log(pchoice_confmat(1)))';        % Fixing B as zero for 1st class gives all the other Bs.
%fit.B_init = (log(pchoice_confmat) - log(pchoice_confmat(end)))';       % Fixing B as zero for final class gives all the other Bs.
%fit.B_init = (log(pchoice_confmat) - log(mean(pchoice_confmat)) )';    % Set mean bias to zero. 
%fit.B_init = (log(pchoice_confmat) )';                                 % Leave untouched. This tends to lead to negative bias.
fit.B_init = (log(pchoice_confmat) - log(pchoice_confmat(B0_ind)) )';   % Reference to the minimum value.

% --- Calculate the inital Z function - the by classficiation metric ----
% We should be able to compute the zs given the Bs?
pi = sum(pconfmat.*eye(length(pconfmat))); % Probability of choosing the right answer - as a matrix
pi(pi==0) =  1e-5;                         % Very small number in place of 0, to prevent log(0).
pi(pi==1) =  1- 1e-5;                      % Subtract very small number off   
fit.Z_init = log(pi) - log(1-pi) + log( sum( exp(fit.B_init))) - fit.B_init;
fit.pmodel_init = softMaxEval(fit.Z_init, fit.B_init);
% The error from the initial parameters. 
fit.rmserr_init = rootmeanfn(softMaxErr([fit.Z_init, fit.B_init],pconfmat),nmat); 
fprintf('\tRMS error for initial parameters: %g\n',fit.rmserr_init);

% Preallocate fitted parameters. 
fit.Z = nan([niter,fit.m]);
fit.B = nan([niter,fit.m]);
fit.err = nan([niter,1]);
fit.fitstoconverge = nan([niter,1]);

if niter>0
   opts = optimset('MaxFunEvals',1000000,'MaxIter',1000000); 
   fprintf('\t%d iterations: ',niter);

   % loop around nboot times
   for fi = 1:niter 
       fprintf('.');
       % Set up initial parameters for this iteration.
       % noisew adds uniform noise to the initial parameters.
       % Fixing one value of B:
       initpars = [fit.Z_init fit.B_init([1:fit.m]~=B0_ind)] + noisew*(rand(1,npars-1)-0.5);
       
       % Run the fit. Attempt it up to maxfiti times. 
       fiti = 1;
       exitflag = 0;
       while exitflag == 0 & fiti<=maxfiti
            [fitpars,fval,exitflag] = fminsearch(@(pars) softMaxErrM_1(pars,pconfmat,smoothw,Zmax,B0_ind),initpars,opts);
            fiti = fiti+1;
       end;

      % Reconstruct the actual model pars from those fitted. 
      fitpars = [fitpars(1:fit.m+B0_ind-1) 0 fitpars(fit.m+B0_ind:2*fit.m-1)];  % If fixing any B value.

      % Hard limit to Zmax. This does very little as few data push the fit
      % above the soft limit imposed by the penaliser in the error
      % function.
      fitpars = min(fitpars,Zmax);
      
      % Store the fitted parameter values.
      fit.Z(fi,:) = fitpars(1:fit.m);
      fit.B(fi,:) = fitpars(fit.m+1:npars);

       % Recalculate the mode (without noise or smooth weighting).
       fit.err(fi) = rootmeanfn(softMaxErr(fitpars,pconfmat),nmat);

       % Return an exit code. -1 indicates did not converge Otherwise
       % returns the number of times it needed to run. 
       if exitflag ~=1
          fit.fitstoconverge(fi) = -1;
          fprintf('\tNo converged fit.\n')
       else
          fit.fitstoconverge(fi) = fiti;
       end;
   end;    
    % endloop
  fprintf('\n');

  % compute the mean parameters for the final fit.
  Z = mean(fit.Z,1); fit.Z_final = Z;
  B = mean(fit.B,1); fit.B_final = B;
else
  Z = fit.Z_init; fit.Z_final = Z;
  B = fit.B_init; fit.B_final = B;
  err = fit.rmserr_init ;
  fprintf('\tReturning intial parameters as fit.\n');
end;

fit.pmodel_final = softMaxEval(Z,B);
err = rootmeanfn(softMaxErr([Z B],pconfmat),nmat);
fit.rmserr_final = err;
fprintf('\tFinal RMS error: %g\n',err);

% Set final average return parameters. 

% Return the calling parameters.
fit.function_pars.smoothw = smoothw;
fit.function_pars.noisew = noisew;
fit.function_pars.niter = niter;




