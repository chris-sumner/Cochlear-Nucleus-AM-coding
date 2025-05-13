function [statsmodels dPrime_reBMFdp dPrime_reBMFdp_reBestDP stats_stacked_4models nltables] = statsModels_general(statsmat_selected,coreoptions,measure,linkfn)
% statsModels_general fits all the statistical models in the paper. 
%
%  [statsmodels nltables] = statsModels_general(statsmat_selected,measure,linkfn)
%
%
% Produces the same format of output for whatever model variant you choose.
% 


%% ------------------- Prepare the data  ---------------------

% Selection of data. Exclude those with modulation depth <.5.
selection_criteria_4models = { ...
    'fn','depthMod>=.5'...
    };

statsmat_selected_4models = select_Datasets(statsmat_selected, selection_criteria_4models{:});
statsmat_selected_4models.dprimes = cellfun(@(x) real(x),statsmat_selected_4models.dprimes,'uni',false);

stats_stacked_4models = stack_Datasets(statsmat_selected_4models,'dprime');

% Make a separate set without any ChS neurons. 
selection_criteria_4models_noChS = { ...
    'fn','depthMod>=.5', ...
    'fn','rationalisedType~="ChS"' ...
    };
statsmat_selected_4models_noChS = select_Datasets(statsmat_selected, selection_criteria_4models_noChS{:});
statsmat_selected_4models_noChS.dprimes = cellfun(@(x) real(x),statsmat_selected_4models_noChS.dprimes,'uni',false);
stats_stacked_4models_noChS = stack_Datasets(statsmat_selected_4models_noChS,'dprime');


%% ---- Derive a mean transfer function across all neuron types ------

% This is the frequency dependence by type.
dPrime_reBMFdp_reBestDP = makePopnTable(stats_stacked_4models,'rationalisedType', ...
    'usedamfreqs_dprime_rationalised','dprimes_pcBestDP', ...
    coreoptions{:},'fn','Rayleigh>13.8','roworder',[1 2 4 5 6 3]);

% Also produce this for plotting.
dPrime_reBMFdp = makePopnTable(stats_stacked_4models,'rationalisedType', ...
    'usedamfreqs_dprime_rationalised','dprimes', ...
    coreoptions{:},'fn','Rayleigh>13.8','roworder',[1 2 4 5 6 3]);

% The mean dependence of normalised d' on modulation frequency.
% This weights the mean according the number of each type.
% It is useful as an intial set of parameters of the frequency function 
meanfModFn = nansum(dPrime_reBMFdp_reBestDP.n.*dPrime_reBMFdp_reBestDP.mean)./nansum(dPrime_reBMFdp_reBestDP.n,1);

% Run the fitting process to generate the "fModFactor" which fits all the
% data best. 

% Specific to model 
switch linkfn
    case ('logistic')
        % This is the logistic version. 
        disp('Link function is logistic');
        fModFactor_fit = fitTheoreticalMaxModel_logistic(stats_stacked_4models,meanfModFn);

        % This gives the p0 value for the fitting. 
        stats_stacked_4models.p0 = fModFactor_fit.p0_stacked;


    case ('linear')
        % Assumes that the underlying function is of the same form, and
        % linearly related, to metric to be predicted. 
        disp('Link function is linear');

        % We fix the in-going function to be 1 for it's lowest frequency 
        % And zero for the lowest score. 
        meanfModFn = (meanfModFn-min(meanfModFn))/(meanfModFn(1)-min(meanfModFn));
        % Run the fit
        fModFactor_fit = fitTheoreticalMaxModel_linear(stats_stacked_4models,meanfModFn);
        
        stats_stacked_4models.fitted_normalised_gains =  nan(size(stats_stacked_4models.dprimes));  
        tmpNormGains =  cellfun(@(x,y) x*ones(1,length(y)), ...
            arrayfun(@(x) x,fModFactor_fit.fitted_normalised_gains,'uni',false), ...
             fModFactor_fit.amfreqs,'uni',false);
        stats_stacked_4models.fitted_normalised_gains([fModFactor_fit.inds{:}]) = [tmpNormGains{:}];
                     
    otherwise
        error('No link function specified');
end;
% End of specific to model.     
    
% Copy to nlm for later plots etc.
nlm.stim_theoreticalmax = fModFactor_fit;
            
% For further model fitting: vector of the correct value depending on fmod:
% This gives the function value for each datapoint (every MTF and frequency). 
fModFactor = fModFactor_fit.normalised_fitted_fn(arrayfun(@(x) find(fModFactor_fit.flist==x), ...
    stats_stacked_4models.usedamfreqs_dprime_rationalised));


%% ---- Derive a mean transfer function without any ChS in it ------

% Mean modulation function with no ChopS neurons
meanfModFn_noChS = nansum(dPrime_reBMFdp_reBestDP.n.*dPrime_reBMFdp_reBestDP.mean)./nansum(dPrime_reBMFdp_reBestDP.n,1);

% Specific to model:
switch linkfn
    case ('logistic')

        % Logistic fit of theoretical best fit to generate best fit underlying function 
        fModFactor_fit_noChS = fitTheoreticalMaxModel_logistic(stats_stacked_4models_noChS,meanfModFn_noChS);
        % This gives the p0 value for the fitting. 
        stats_stacked_4models_noChS.p0 = fModFactor_fit_noChS.p0_stacked;

    case ('linear')
        % Assumes that the underlying function is of the same form, and
        % linearly related, to metric to be predicted. 

        % We fix the in-going function to be 1 for it's lowest frequency 
        % And zero for the lowest score. 
        meanfModFn_noChS = meanfModFn_noChS-min(meanfModFn_noChS)/(meanfModFn_noChS(1)-min(meanfModFn_noChS));
        % Run the fit
        fModFactor_fit_noChS = fitTheoreticalMaxModel_linear(stats_stacked_4models_noChS,meanfModFn_noChS);
                   
    otherwise
        error('No link function specified');
end;
% end of specific to model.    
    
% Extract a list of the factors for every data point.
fModFactor_noChS =  fModFactor_fit_noChS.normalised_fitted_fn(arrayfun(@(x) find(fModFactor_fit_noChS.flist==x), ...
    stats_stacked_4models_noChS.usedamfreqs_dprime_rationalised));

% Copy to nlm. In case we want to look at it. 
nlm.stim_theoreticalmax_noChS = fModFactor_fit_noChS;



%% ----- Make the tables for the statistical models to fit to ----

% This is not very efficient. It grew...

% This is one to check the fitting - it should produce a coefficient close
% to one on the gain parameter.
tables.reproduceMax = table(stats_stacked_4models.fitted_normalised_gains', ...
                    fModFactor',stats_stacked_4models.dprimes', ...
           'VariableNames',{'normGain','fFactor','dprimes'});  
       
% Sitmulus only.
tables.stimonly = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', ...
           fModFactor', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','dprimes'});  
       
% For one additional continuous variable the table format is the same.
% Using this we can have the stimulus plus any one predictor we want. 
tables.var1 = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', ...
           fModFactor', stats_stacked_4models.VS', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','Var1','dprimes'});  

% For two additional continuous variables.
tables.var2 = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', ...
           fModFactor', stats_stacked_4models.VS',stats_stacked_4models.VS', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','dprimes'});  

% For three additional continuous variables.
tables.var3 = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', ...
           fModFactor', stats_stacked_4models.VS',stats_stacked_4models.VS', stats_stacked_4models.VS', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','dprimes'});  

% For four additional continuous variables.
tables.var4 = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', ...
           fModFactor', stats_stacked_4models.VS',stats_stacked_4models.VS', stats_stacked_4models.VS', ...
           stats_stacked_4models.VS', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','dprimes'});  

% For five additional continuous variables.
tables.var5 = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', ...
           fModFactor', stats_stacked_4models.VS',stats_stacked_4models.VS', stats_stacked_4models.VS', ...
           stats_stacked_4models.VS', stats_stacked_4models.VS', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5','dprimes'});  
       
% For TEN additional continuous variables (which is all the ones you could want).
tables.var10 = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', fModFactor', ...
           stats_stacked_4models.VS', stats_stacked_4models.VS', stats_stacked_4models.VS', stats_stacked_4models.VS', ...
           stats_stacked_4models.VS', stats_stacked_4models.VS',stats_stacked_4models.VS', stats_stacked_4models.VS', ...
           stats_stacked_4models.VS',stats_stacked_4models.VS', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5', ...
           'Var6','Var7','Var8','Var9','Var10','dprimes'});  
       
% Type is a bit involved to make it work. Each type is a level in the
% factor. 
tables.type = tables.var1;   
tables.type.Var1(1:length(statsmat_selected_4models.rationalisedType)) =   (   1*strcmp(statsmat_selected_4models.rationalisedType,'ChS') +  2*strcmp(statsmat_selected_4models.rationalisedType,'ChT') +  ...
                          3*strcmp(statsmat_selected_4models.rationalisedType,'On') +   4*strcmp(statsmat_selected_4models.rationalisedType,'PBU') +  ...
                          5*strcmp(statsmat_selected_4models.rationalisedType,'PL') +   6*strcmp(statsmat_selected_4models.rationalisedType,'PLN') )';

% We make individual tables for each predictor variable.  
% Basic statistics where predictors are drawn from the spike train for each
% modulation freqency. Not CV obviously - that is the CV of the pure tone
% response.
tables.Z = tables.var1;   tables.Z.Var1 =  stats_stacked_4models.Z';
tables.VS = tables.var1;   tables.VS.Var1 =  stats_stacked_4models.VS';
tables.CI = tables.var1;   tables.CI.Var1 =  stats_stacked_4models.CI';
tables.CV = tables.var1;   tables.CV.Var1 =  stats_stacked_4models.CV';
tables.SPP = tables.var1;   tables.SPP.Var1 =  stats_stacked_4models.spikesperperiod';
tables.SPP = tables.var1;   tables.SPP.Var1 =  stats_stacked_4models.spikesperperiod';
tables.spikerate = tables.var1;   tables.spikerate.Var1 =  stats_stacked_4models.spikesperperiod'.* ...
                                                            stats_stacked_4models.usedamfreqs_dprime';
tables.numSACpeaks = tables.var1; tables.numSACpeaks.Var1 =  stats_stacked_4models.numSACpeaks_p0001';
tables.reliability = tables.var1; tables.reliability.Var1 =  stats_stacked_4models.reliability_R_0p24ms';
tables.envFluct = tables.var1; tables.envFluct.Var1 =  stats_stacked_4models.envFluct_0p48ms';
% For the SAC salience (peak-to-trough) we run a selection of the possible models.
tables.SACpeak1Sal_p0001 = tables.var1; tables.SACpeak1Sal_p0001.Var1 =  stats_stacked_4models.SACpeak1sal_p0001';
tables.SACpeak234Sal_p0001 = tables.var3; tables.SACpeak234Sal_p0001.Var1 =  stats_stacked_4models.SACpeak2sal_p0001';
    tables.SACpeak234Sal_p0001.Var2 =  stats_stacked_4models.SACpeak3sal_p0001';
    tables.SACpeak234Sal_p0001.Var3 =  stats_stacked_4models.SACpeak4sal_p0001';
tables.SACpeak1234Sal_p0001 = tables.var4; tables.SACpeak1234Sal_p0001.Var1 =  stats_stacked_4models.SACpeak1sal_p0001';
    tables.SACpeak1234Sal_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_p0001';
    tables.SACpeak1234Sal_p0001.Var3 =  stats_stacked_4models.SACpeak3sal_p0001';
    tables.SACpeak1234Sal_p0001.Var4 =  stats_stacked_4models.SACpeak4sal_p0001';
tables.SACpeak12345Sal_p0001 = tables.var5; tables.SACpeak12345Sal_p0001.Var1 =  stats_stacked_4models.SACpeak1sal_p0001';
    tables.SACpeak12345Sal_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_p0001';
    tables.SACpeak12345Sal_p0001.Var3 =  stats_stacked_4models.SACpeak3sal_p0001';
    tables.SACpeak12345Sal_p0001.Var4 =  stats_stacked_4models.SACpeak4sal_p0001';
    tables.SACpeak12345Sal_p0001.Var5 =  stats_stacked_4models.SACpeak5sal_p0001';
tables.SACpeak2345Sal_p0001 = tables.var4; tables.SACpeak2345Sal_p0001.Var1 =  stats_stacked_4models.SACpeak2sal_p0001';
    tables.SACpeak12345Sal_p0001.Var2 =  stats_stacked_4models.SACpeak3sal_p0001';
    tables.SACpeak12345Sal_p0001.Var3 =  stats_stacked_4models.SACpeak4sal_p0001';
    tables.SACpeak12345Sal_p0001.Var4 =  stats_stacked_4models.SACpeak5sal_p0001';
    
% "150Hz or thereabouts statistics" - 
% The statistics are taken from the rationalised 150Hz responses.
tables.Z_150ish = tables.var1;   tables.Z_150ish.Var1 =  stats_stacked_4models.Z_150ish';
tables.VS_150ish = tables.var1;   tables.VS_150ish.Var1 =  stats_stacked_4models.VS_150ish';
tables.CI_150ish = tables.var1;   tables.CI_150ish.Var1 =  stats_stacked_4models.CI_150ish';
tables.spikerate_150ish = tables.var1;   tables.spikerate_150ish.Var1 =  stats_stacked_4models.spikesperperiod_150ish'.* ...
                                                                             stats_stacked_4models.usedamfreqs_dprime';
tables.numSACpeaks_150ish = tables.var1; tables.numSACpeaks_150ish.Var1 =  stats_stacked_4models.numSACpeaks_p0001_150ish';
tables.reliability_150ish = tables.var1; tables.reliability_150ish.Var1 =  stats_stacked_4models.reliability_R150ish_0p24ms';
tables.envFluct_150ish = tables.var1; tables.envFluct_150ish.Var1 =  stats_stacked_4models.envFluct_R150ish_0p48ms';

%  ----- models looking at the salience of peaks at 150-ish Hz -----

% Zero peak (~correlation index).
tables.SACpeak0Sal_150ish_p0001 = tables.var1; tables.numSACpeaks_150ish.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';

% Zero peak + peak 1. 
tables.SACpeak01Sal_150ish_p0001 = tables.var2; 
tables.SACpeak01Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak01Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';        

% Zero peak + peak 1. 
tables.SACpeak012Sal_150ish_p0001 = tables.var3; 
tables.SACpeak012Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak012Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';        
tables.SACpeak012Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        

% Zero peak + peak 1,2,3. 
tables.SACpeak0123Sal_150ish_p0001 = tables.var4; 
tables.SACpeak0123Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak0123Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';        
tables.SACpeak0123Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak0123Sal_150ish_p0001.Var4 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        

% Zero peak + peak 1,2,3,4. 
tables.SACpeak01234Sal_150ish_p0001 = tables.var5; 
tables.SACpeak01234Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak01234Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';        
tables.SACpeak01234Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak01234Sal_150ish_p0001.Var4 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        
tables.SACpeak01234Sal_150ish_p0001.Var5 =  stats_stacked_4models.SACpeak4sal_150ish_p0001';        

% For the salience of all 6 peaks we make a table from scratch as it has no other use. 
tables.SACpeakAllSals_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.SACpeak0sal_150ish_p0001', stats_stacked_4models.SACpeak1sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak2sal_150ish_p0001', stats_stacked_4models.SACpeak3sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak4sal_150ish_p0001', stats_stacked_4models.SACpeak5sal_150ish_p0001', ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','SACsal_0','SACsal_1','SACsal_2','SACsal_3','SACsal_4','SACsal_5','dprimes'});  

% Starting with peak 1, since 0 and 1 are highly correlated. 
tables.SACpeak1Sal_150ish_p0001 = tables.var1; 
tables.SACpeak1Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';        

%  peak 1,2. 
tables.SACpeak12Sal_150ish_p0001 = tables.var2; 
tables.SACpeak12Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';        
tables.SACpeak12Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        

%  peak 1,2,3. 
tables.SACpeak123Sal_150ish_p0001 = tables.var3; 
tables.SACpeak123Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';        
tables.SACpeak123Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak123Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        

%  peak 1,2,3,4. 
tables.SACpeak1234Sal_150ish_p0001 = tables.var4; 
tables.SACpeak1234Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';        
tables.SACpeak1234Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak1234Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        
tables.SACpeak1234Sal_150ish_p0001.Var4 =  stats_stacked_4models.SACpeak4sal_150ish_p0001';        

% peaks 1-5.
tables.SACpeak12345Sal_150ish_p0001 = tables.var5; 
tables.SACpeak12345Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak1sal_150ish_p0001';
tables.SACpeak12345Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak12345Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        
tables.SACpeak12345Sal_150ish_p0001.Var4 =  stats_stacked_4models.SACpeak4sal_150ish_p0001';        
tables.SACpeak12345Sal_150ish_p0001.Var5 =  stats_stacked_4models.SACpeak5sal_150ish_p0001';        

% peaks 2,3,4 (no peak at the period). 
tables.SACpeak23Sal_150ish_p0001 = tables.var2; 
tables.SACpeak23Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak23Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        

% peaks 2,3,4 (no peak at the period). 
tables.SACpeak234Sal_150ish_p0001 = tables.var3; 
tables.SACpeak234Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak234Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        
tables.SACpeak234Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak4sal_150ish_p0001';        

% peaks 2,3,4,5 (no peak at the period). 
tables.SACpeak2345Sal_150ish_p0001 = tables.var4; 
tables.SACpeak2345Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak2345Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        
tables.SACpeak2345Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak4sal_150ish_p0001';        
tables.SACpeak2345Sal_150ish_p0001.Var4 =  stats_stacked_4models.SACpeak5sal_150ish_p0001';        



% Zero peak + peak 2. Peak 1 is adding nothing so we omit it. 
tables.SACpeak0Plus2Sal_150ish_p0001 = tables.var2; 
tables.SACpeak0Plus2Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak0Plus2Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        

% Zero peak + peaks 2,3. 
tables.SACpeak0Plus23Sal_150ish_p0001 = tables.var3; 
tables.SACpeak0Plus23Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak0Plus23Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak0Plus23Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        

% Zero peak + peaks 2 +4. 
tables.SACpeak0Plus24Sal_150ish_p0001 = tables.var3; 
tables.SACpeak0Plus24Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak0Plus24Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak0Plus24Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak4sal_150ish_p0001';        

% Zero peak + peaks 2,3,4. 
tables.SACpeak0Plus234Sal_150ish_p0001 = tables.var4; 
tables.SACpeak0Plus234Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak0Plus234Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak0Plus234Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        
tables.SACpeak0Plus234Sal_150ish_p0001.Var4 =  stats_stacked_4models.SACpeak4sal_150ish_p0001';        

% Zero peak + peaks 2,3,4,5. 
tables.SACpeak0Plus2345Sal_150ish_p0001 = tables.var5; 
tables.SACpeak0Plus2345Sal_150ish_p0001.Var1 =  stats_stacked_4models.SACpeak0sal_150ish_p0001';
tables.SACpeak0Plus2345Sal_150ish_p0001.Var2 =  stats_stacked_4models.SACpeak2sal_150ish_p0001';        
tables.SACpeak0Plus2345Sal_150ish_p0001.Var3 =  stats_stacked_4models.SACpeak3sal_150ish_p0001';        
tables.SACpeak0Plus2345Sal_150ish_p0001.Var4 =  stats_stacked_4models.SACpeak4sal_150ish_p0001';        
tables.SACpeak0Plus2345Sal_150ish_p0001.Var5 =  stats_stacked_4models.SACpeak5sal_150ish_p0001';        


% This is the model with peaks 1-4 in it, but no sustained choppers. 
tables.SACpeak12345Sal_150ish_p0001_noChS = table(stats_stacked_4models_noChS.modLvl_rationalised', ...
            stats_stacked_4models_noChS.depthMod', fModFactor_noChS', ...
            stats_stacked_4models_noChS.SACpeak1sal_150ish_p0001', ...
            stats_stacked_4models_noChS.SACpeak2sal_150ish_p0001', stats_stacked_4models_noChS.SACpeak3sal_150ish_p0001', ...
            stats_stacked_4models_noChS.SACpeak4sal_150ish_p0001', stats_stacked_4models_noChS.SACpeak5sal_150ish_p0001',  ...
            stats_stacked_4models_noChS.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor', ...
             'Var1','Var2','Var3','Var4','Var5','dprimes'});  

tables.Z_150ish_noChS = table(stats_stacked_4models_noChS.modLvl_rationalised', ...
            stats_stacked_4models_noChS.depthMod', fModFactor_noChS', ...
            stats_stacked_4models_noChS.Z_150ish',  ...
            stats_stacked_4models_noChS.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor', ...
             'Var1','dprimes'});  
         
tables.reliability_150ish_noChS = table(stats_stacked_4models_noChS.modLvl_rationalised', ...
            stats_stacked_4models_noChS.depthMod', fModFactor_noChS', ...
            stats_stacked_4models_noChS.reliability_R150ish_0p24ms',  ...
            stats_stacked_4models_noChS.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor', ...
             'Var1','dprimes'});  
         

% ----- A model with all the high performing individual statistics -----

% For the salience of all 4 SAC peaks + Z, NuFl and releiability. 
tables.AllSensibleVars_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', ...
            stats_stacked_4models.envFluct_R150ish_0p48ms',stats_stacked_4models.reliability_R150ish_0p24ms', ...
            stats_stacked_4models.SACpeak1sal_150ish_p0001', stats_stacked_4models.SACpeak2sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak3sal_150ish_p0001', stats_stacked_4models.SACpeak4sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Z','NuFl','rel' ...
            'SACsal_1','SACsal_2','SACsal_3','SACsal_4','dprimes'});  

% ------ models with 1 parameter removed  ---------        

tables.AllSensibleVars_minusSACPeaks_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', ...
            stats_stacked_4models.envFluct_R150ish_0p48ms',stats_stacked_4models.reliability_R150ish_0p24ms', ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','dprimes'});  


% For the salience of the period peak, Z, NuFl and reliability. 
tables.AllSensibleVars_minusExtraSACPeaks_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', ...
            stats_stacked_4models.envFluct_R150ish_0p48ms',stats_stacked_4models.reliability_R150ish_0p24ms', ...
            stats_stacked_4models.SACpeak1sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','dprimes'});  
        
% For the salience of the extra peaks + Z, NuFl and releiability. 
tables.AllSensibleVars_minusSACpeak0_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', ...
            stats_stacked_4models.envFluct_R150ish_0p48ms',stats_stacked_4models.reliability_R150ish_0p24ms', ...
            stats_stacked_4models.SACpeak2sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak3sal_150ish_p0001', stats_stacked_4models.SACpeak4sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5','Var6','dprimes'});          

% For the salience of all the sensible parameters minus Z 
tables.AllSensibleVars_minusZ_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.envFluct_R150ish_0p48ms',stats_stacked_4models.reliability_R150ish_0p24ms', ...
             stats_stacked_4models.SACpeak1sal_150ish_p0001',  stats_stacked_4models.SACpeak2sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak3sal_150ish_p0001', stats_stacked_4models.SACpeak4sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5','Var6','dprimes'});          

% For the salience of all the sensible parameters minus Z 
tables.AllSensibleVars_minusNuFl_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', ...
            stats_stacked_4models.reliability_R150ish_0p24ms', ...
            stats_stacked_4models.SACpeak1sal_150ish_p0001',  stats_stacked_4models.SACpeak2sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak3sal_150ish_p0001', stats_stacked_4models.SACpeak4sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5','Var6','dprimes'});          
        
% For the salience of all the sensible parameters minus Z 
tables.AllSensibleVars_minusReliability_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', stats_stacked_4models.envFluct_R150ish_0p48ms',...
            stats_stacked_4models.SACpeak1sal_150ish_p0001',  stats_stacked_4models.SACpeak2sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak3sal_150ish_p0001', stats_stacked_4models.SACpeak4sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5','Var6','dprimes'});          
        
% ----------- pairs of parameters -----------

% SAC peaks and Z
tables.TwoVars_PeaksPlusZ_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', ...
            stats_stacked_4models.SACpeak1sal_150ish_p0001', stats_stacked_4models.SACpeak2sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak3sal_150ish_p0001', stats_stacked_4models.SACpeak4sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5','dprimes'});  

% SAC peaks plus neural fluctuation        
tables.TwoVars_PeaksPlusNuFl_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.envFluct_R150ish_0p48ms', ...
            stats_stacked_4models.SACpeak1sal_150ish_p0001', stats_stacked_4models.SACpeak2sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak3sal_150ish_p0001', stats_stacked_4models.SACpeak4sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5','dprimes'});  

% SAC peaks plus reliability        
tables.TwoVars_PeaksPlusReliability_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.reliability_R150ish_0p24ms', ...
            stats_stacked_4models.SACpeak1sal_150ish_p0001', stats_stacked_4models.SACpeak2sal_150ish_p0001', ...
            stats_stacked_4models.SACpeak3sal_150ish_p0001', stats_stacked_4models.SACpeak4sal_150ish_p0001',  ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','Var3','Var4','Var5','dprimes'});  

% Z plus neural fluctuation        
tables.TwoVars_ZPlusNuFl_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', ...
            stats_stacked_4models.envFluct_R150ish_0p48ms', ...            
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','dprimes'});  
        
% Z plus reliability        
tables.TwoVars_ZPlusReliability_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.Z_150ish', ...
            stats_stacked_4models.reliability_R150ish_0p24ms', ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','dprimes'});  

%neural fluctuation plus reliability        
tables.TwoVars_NuFlPlusReliability_150ish = table(stats_stacked_4models.modLvl_rationalised', ...
            stats_stacked_4models.depthMod', fModFactor', ...
            stats_stacked_4models.envFluct_R150ish_0p48ms', ...
            stats_stacked_4models.reliability_R150ish_0p24ms', ...
            stats_stacked_4models.dprimes', ...
            'VariableNames',{'modLvl','depthMod','fFactor','Var1','Var2','dprimes'});  

             
        
%% --------- Make the equation forms to use in the models -------             
% N.B. removed the highest mod-level and smallest depth from these, as they
% for three values there should only be 2 degrees of freedom (ANOVA style).
% Likewise there are 7 types but only 6 levels. 
% This has virtually no impact on r^2
% These are the "core" linear regression parts of the model.

eqns.reproduceMax = 'v1*normGain';
eqns.var1 = 'v1*Var1 + l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)'; 
eqns.var2 = 'v1*Var1 +  v2*Var2 + l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)'; 
eqns.var3 = 'v1*Var1 +  v2*Var2 +  v3*Var3 + l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)'; 
eqns.var4 = 'v1*Var1 +  v2*Var2 +  v3*Var3 +  v4*Var4 + l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)'; 
eqns.var5 = 'v1*Var1 +  v2*Var2 +  v3*Var3 +  v4*Var4 +  v5*Var5 + l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)'; 
eqns.var6 = 'v1*Var1 +  v2*Var2 +  v3*Var3 +  v4*Var4 +  v5*Var5 +  v6*Var6 + l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)'; 

eqns.var10 = ['v1*Var1 +  v2*Var2 +  v3*Var3 +  v4*Var4 +  v5*Var5 + ' ...
                'v6*Var6 +  v7*Var7 +  v8*Var8 +  v9*Var9 +  v5*Var10 + ' ...
                'l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)']; 
eqns.stimonly = 'l1*(modLvl==30) + l2*(modLvl==50)   + d2*(depthMod==1) + d3*(depthMod==2)'; 
eqns.type = ['l1*(modLvl==30) + l2*(modLvl==50)  + d2*(depthMod==1) + d3*(depthMod==2) +' ...
               't1*(Var1==1) +  t2*(Var1==2) + t3*(Var1==3) + t4*(Var1==4) +  t5*(Var1==5) + ' ...
               't6*(Var1==6)'];

% Model for SAC peak salience - long so split
eqns.AllSalsStr = 's0*SACsal_0 + s1*SACsal_1 + s2*SACsal_2 + s3*SACsal_3 + s4*SACsal_4 + s5*SACsal_5 + ';
eqns.SACpeakAllSals =[ eqns.AllSalsStr ' l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)']; 

% Model for all the sensible (vaguely promising) variables - long so split.
eqns.AllSensibleVarsStr = [ 'z*Z + nf*NuFl + r*rel + s1*SACsal_1 + s2*SACsal_2 + s3*SACsal_3 + s4*SACsal_4 + '];
eqns.AllSensibleVars = [ eqns.AllSensibleVarsStr  ' l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2)']; 

% Starting parameter fits to use. 
beta0.reproduceMax = [1];
beta0.stimonly = [3 2 2 1.8];
beta0.type = [3 2 2 1.8 2 2 1 1 1 1];
beta0.var1 = [1 3 2 2 1.8];
beta0.var2 = [1 2 3 2 2 1.8];
beta0.var3 = [1 1 2 3 2 2 1.8];
beta0.var4 = [1 1 1 1 3 2 2 1.8];
beta0.var5 = [1 1 1 1 1 3 2 2 1.8];
beta0.var6 = [1 1 1 1 1 1 3 2 2 1.8];
beta0.var7 = [1 1 1 1 1 1 1 3 2 2 1.8];
beta0.var9 = [1 1 1 1 1 1 1 1 1 3 2 2 1.8];
beta0.var10 =[1 1 1 1 1 1 1 1 1 1 3 2 2 1.8];

%% -------------- Fitting the models ---------------

% Code which is specific to the logistic model
switch linkfn
    case ('logistic')
        % Logistic is slightly ad-hoc, as it makes a correction for chance,
        % dependent on the number of classes (via p0 parameter) but it is
        % built on a logistic function which assumes a 2-class task. It was
        % a warm-up really.

        % Function which does any prepraration required for the table. 
        % This makes it simply to alter the tables without changing their
        % definition.
        % For logistic function it adds p0 to the table.
        tablePrepFn = @(intable) addvars(intable,stats_stacked_4models.p0','NewVariableNames','p0');
        tablePrepFn_noChS = @(intable) addvars(intable,stats_stacked_4models_noChS.p0','NewVariableNames','p0');

        % Preparing the tables. Adding p0. 
        tablenames = fieldnames(tables);
        for ti = 1:length(tablenames)
            if ~isempty(strfind(tablenames{ti},'noChS'))
                nltables.(tablenames{ti}) = tablePrepFn_noChS(tables.(tablenames{ti}));
            else
                nltables.(tablenames{ti}) = tablePrepFn(tables.(tablenames{ti}));
            end;
        end;

        % Function prepares equation, taking the linear regression part of the
        % equation as input.tbales
        % This makes it easy to change the non-linear parts, and the link function.
        eqnPrepFn = @(leqn) ['dprimes ~ logisticLinkFn(fFactor*(' leqn '),p0,1)'];
        nleqns = structfun(@(s) eqnPrepFn(s),eqns,'UniformOutput',false);
        
    case ('linear')
        % Linear version 
        % This is logical to use directly on transforms of the clasifier
        % performance which should in theory correct for the stimulus set. 
        % So directly on the softmax or d' calculations for each dataset. 
        
        % There are no new variables to add but we need to check for inf and nan values.
        tablenames = fieldnames(tables);
        for ti = 1:length(tablenames)
            nltables.(tablenames{ti}) = removeNanInfRows(tables.(tablenames{ti}));
        end;
        
        % Function prepares equations. Only need add the frequency factor
        % for the linear model. 
        eqnPrepFn = @(leqn) ['dprimes ~ fFactor*(' leqn ')'];
        nleqns = structfun(@(s) eqnPrepFn(s),eqns,'UniformOutput',false);

    otherwise
        error('No link function specified');
end;
% end of code specific to different link functions.     

% Try to reproduce the theoretical maximum
% If the fitting is working it produces a similar fit to theoretical max, and coefficient of 1.
nlm.reproduceMax = fitnlm(nltables.reproduceMax,nleqns.reproduceMax,beta0.reproduceMax);

% Stimulus conditions only
nlm.stimonly = fitnlm(nltables.stimonly,nleqns.stimonly,beta0.stimonly);

% Stimulus conditions + type. 
nlm.type = fitnlm(nltables.type,nleqns.type,beta0.type)

% Stimulus + predictor statistic drawn from the same spike train (except
% CV).
nlm.CI = fitnlm(nltables.CI,nleqns.var1,beta0.var1);
nlm.VS = fitnlm(nltables.VS,nleqns.var1,beta0.var1);
nlm.Z = fitnlm(nltables.Z,nleqns.var1,beta0.var1);
nlm.SPP = fitnlm(nltables.SPP,nleqns.var1,beta0.var1);
%nlm.spikerate = fitnlm(nltables.spikerate,nleqns.var1,beta0.var1);
nlm.CV = fitnlm(nltables.CV,nleqns.var1,beta0.var1);
nlm.numSACpeaks = fitnlm(nltables.numSACpeaks,nleqns.var1,beta0.var1);
nlm.reliability = fitnlm(nltables.reliability,nleqns.var1,beta0.var1);
nlm.envFluct = fitnlm(nltables.envFluct,nleqns.var1,beta0.var1);
% We don't run all the peak salience models for the statistics drawn from the same spike train - only key examples.
nlm.SACpeaks1Sal_p0001 = fitnlm(nltables.SACpeak1Sal_p0001,nleqns.var1,beta0.var1);
nlm.SACpeaks234Sal_p0001 = fitnlm(nltables.SACpeak234Sal_p0001,nleqns.var3,beta0.var3);
nlm.SACpeaks1234Sal_p0001 = fitnlm(nltables.SACpeak12345Sal_p0001,nleqns.var4,beta0.var4);
nlm.SACpeaks12345Sal_p0001 = fitnlm(nltables.SACpeak12345Sal_p0001,nleqns.var5,beta0.var5);
nlm.SACpeaks2345Sal_p0001 = fitnlm(nltables.SACpeak2345Sal_p0001,nleqns.var4,beta0.var4);

% Single values near 150Hz.
nlm.Z_150ish = fitnlm(nltables.Z_150ish,nleqns.var1,beta0.var1);
nlm.VS_150ish = fitnlm(nltables.VS_150ish,nleqns.var1,beta0.var1);
nlm.CI_150ish = fitnlm(nltables.CI_150ish,nleqns.var1,beta0.var1);
nlm.spikerate_150ish = fitnlm(nltables.spikerate_150ish,nleqns.var1,beta0.var1);
nlm.numSACpeaks_150ish = fitnlm(nltables.numSACpeaks_150ish,nleqns.var1,beta0.var1);
nlm.reliability_150ish = fitnlm(nltables.reliability_150ish,nleqns.var1,beta0.var1);
nlm.envFluct_150ish = fitnlm(nltables.envFluct_150ish,nleqns.var1,beta0.var1);

% Progressively adding the salience of additional peaks. 
nlm.SACpeak0Sal_150ish_p0001 = fitnlm(nltables.SACpeak0Sal_150ish_p0001,nleqns.var1,beta0.var1);
nlm.SACpeaks01Sal_150ish_p0001 = fitnlm(nltables.SACpeak01Sal_150ish_p0001,nleqns.var2,beta0.var2);
nlm.SACpeaks012Sal_150ish_p0001 = fitnlm(nltables.SACpeak012Sal_150ish_p0001,nleqns.var3,beta0.var3);
nlm.SACpeaks0123Sal_150ish_p0001 = fitnlm(nltables.SACpeak0123Sal_150ish_p0001,nleqns.var4,beta0.var4);
nlm.SACpeaks01234Sal_150ish_p0001 = fitnlm(nltables.SACpeak01234Sal_150ish_p0001,nleqns.var5,beta0.var5);
nlm.SACpeakAllSals_150ish = fitnlm(nltables.SACpeakAllSals_150ish,nleqns.SACpeakAllSals,beta0.var6);

% Omitting peak 0, which is highly correlated with peak 1.
nlm.SACpeak1Sal_150ish_p0001 = fitnlm(nltables.SACpeak1Sal_150ish_p0001,nleqns.var1,beta0.var1);
nlm.SACpeak12Sal_150ish_p0001 = fitnlm(nltables.SACpeak12Sal_150ish_p0001,nleqns.var2,beta0.var2);
nlm.SACpeak123Sal_150ish_p0001 = fitnlm(nltables.SACpeak123Sal_150ish_p0001,nleqns.var3,beta0.var3);
nlm.SACpeak1234Sal_150ish_p0001 = fitnlm(nltables.SACpeak1234Sal_150ish_p0001,nleqns.var4,beta0.var4);
nlm.SACpeak12345Sal_150ish_p0001 = fitnlm(nltables.SACpeak12345Sal_150ish_p0001,nleqns.var5,beta0.var5);

% Alternative sequential adding of peaks with 1 (highly correlated with 0) removed. 
nlm.SACpeak0Plus2Sal_150ish_p0001 = fitnlm(nltables.SACpeak0Plus2Sal_150ish_p0001,nleqns.var2,beta0.var2);
nlm.SACpeak0Plus23Sal_150ish_p0001 = fitnlm(nltables.SACpeak0Plus23Sal_150ish_p0001,nleqns.var3,beta0.var3);
nlm.SACpeak0Plus24Sal_150ish_p0001 = fitnlm(nltables.SACpeak0Plus24Sal_150ish_p0001,nleqns.var3,beta0.var3);
nlm.SACpeak0Plus234Sal_150ish_p0001 = fitnlm(nltables.SACpeak0Plus234Sal_150ish_p0001,nleqns.var4,beta0.var4);
% This last model does not fit well. 
nlm.SACpeak0Plus2345Sal_150ish_p0001 = fitnlm(nltables.SACpeak0Plus2345Sal_150ish_p0001,nleqns.var5,beta0.var5);

% And removing the first two (zero lag and at the period).
nlm.SACpeak23Sal_150ish_p0001 = fitnlm(nltables.SACpeak23Sal_150ish_p0001 ,nleqns.var2,beta0.var2);
nlm.SACpeak234Sal_150ish_p0001 = fitnlm(nltables.SACpeak234Sal_150ish_p0001,nleqns.var3,beta0.var3);
nlm.SACpeak2345Sal_150ish_p0001 = fitnlm(nltables.SACpeak2345Sal_150ish_p0001,nleqns.var4,beta0.var4);

% However adding all the good indivudal variables together makes a still better model.  
nlm.AllSensibleVars_150ish = fitnlm(nltables.AllSensibleVars_150ish,nleqns.AllSensibleVars,beta0.var7);

% Since mode-locking perse only exists in ChS neurons, we ask how these
% key predictors fare when there are no sustained choppers in the dataset.
nlm.SACpeak12345Sal_150ish_p0001_noChS = fitnlm(nltables.SACpeak12345Sal_150ish_p0001_noChS,nleqns.var5,beta0.var5);
nlm.Z_150ish_noChS = fitnlm(nltables.Z_150ish_noChS,nleqns.var1,beta0.var1);
nlm.reliability_150ish_noChS = fitnlm(nltables.reliability_150ish_noChS,nleqns.var1,beta0.var1);



% ----- Removing parameters from the "All sensible" one at at time ------ 
% None of these are critical. 
% Removing the extra peaks (no difference)
nlm.AllSensibleVars_minusExtraSACPeaks_150ish = ...
    fitnlm(nltables.AllSensibleVars_minusExtraSACPeaks_150ish,nleqns.var4,beta0.var4);
% removing the period peak (no difference)  
nlm.AllSensibleVars_minusSACpeak0_150ish = ...
    fitnlm(nltables.AllSensibleVars_minusSACpeak0_150ish,nleqns.var6,beta0.var6);
% Removing all the peaks. 
nlm.AllSensibleVars_minusSACPeaks_150ish = ...
    fitnlm(nltables.AllSensibleVars_minusSACPeaks_150ish,nleqns.var3,beta0.var3);
% removing Z (very small difference)  
nlm.AllSensibleVars_minusZ_150ish = ...
    fitnlm(nltables.AllSensibleVars_minusZ_150ish,nleqns.var6,beta0.var6);
% removing neural fluctuation (very small difference)  
nlm.AllSensibleVars_minusNuFl_150ish = ...
    fitnlm(nltables.AllSensibleVars_minusNuFl_150ish,nleqns.var6,beta0.var6);
% removing reliability (very, very  small difference)  
nlm.AllSensibleVars_minusReliability_150ish = ...
    fitnlm(nltables.AllSensibleVars_minusReliability_150ish,nleqns.var6,beta0.var6);

% ---- Models with pairs of parameters ----
% Any two parameters pretty much falls between any 3. 

% SAC peaks + Z
nlm.TwoVars_PeaksPlusZ_150ish = ...
    fitnlm(nltables.TwoVars_PeaksPlusZ_150ish,nleqns.var5,beta0.var5);
% SAC peaks + NuFl
nlm.TwoVars_PeaksPlusNuFl_150ish = ...
    fitnlm(nltables.TwoVars_PeaksPlusNuFl_150ish,nleqns.var5,beta0.var5);
% SAC peaks + rel. 
nlm.TwoVars_PeaksPlusReliability_150ish = ...
    fitnlm(nltables.TwoVars_PeaksPlusReliability_150ish,nleqns.var5,beta0.var5);
% Z + Nufl
nlm.TwoVars_ZPlusNuFl_150ish = ...
    fitnlm(nltables.TwoVars_ZPlusNuFl_150ish,nleqns.var2,beta0.var2);
% Z + rel
nlm.TwoVars_ZPlusReliability_150ish = ...
    fitnlm(nltables.TwoVars_ZPlusReliability_150ish,nleqns.var2,beta0.var2);
% NuFl + rel
nlm.TwoVars_NuFlPlusReliability_150ish = ...
    fitnlm(nltables.TwoVars_NuFlPlusReliability_150ish,nleqns.var2,beta0.var2);

% Make the output structure. 
statsmodels = struct('linkfn',linkfn,'measure',measure, ...
    'meanfModFn',meanfModFn,'nleqns',nleqns,'nlm',nlm,'nltables',nltables, ...
    'dPrime_reBMFdp',dPrime_reBMFdp, 'dPrime_reBMFdp_reBestDP',dPrime_reBMFdp_reBestDP, ...
    'stats_stacked_4models',stats_stacked_4models);

% Saving for github we remove the tables of the data to be fit. They are
% not needed. 
statsmodels = rmfield(statsmodels,'nltables');


function tableout = removeNanInfRows(tablein)

inds = find( prod(~table2array(varfun(@isnan, tablein))' & ...
            ~table2array(varfun(@isinf, tablein) )'  ) );
tableout  = tablein(inds,:);






