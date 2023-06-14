%% ------------------- stats models ---------------------

% New selection of data. Exclude those with modulation depth <.5.
selection_criteria_4models = { ...
    'fn','depthMod>=.5'...
    };
statsmat_selected_4models = select_Datasets(statsmat_selected, selection_criteria_4models{:});
stats_stacked_4models = stack_Datasets(statsmat_selected_4models,'dprime');

% This is a non-linear model:
% d' = d'max * f(fMod)
% Whereas we have been using a linear model, which adds a different static 
% amount dependent on frequency.
% d' = x + a(fMod)

%% ---- Derive a mean transfer function across all neuron types ------

% This is the frequency dependence by type.
dPrime_reBMFdp_reBestDP = makePopnTable(stats_stacked_4models,'rationalisedType', ...
    'usedamfreqs_dprime_rationalised','dprimes_pcBestDP', ...
    coreoptions{:},'fn','Rayleigh>13.8','roworder',[1 2 4 5 6 3]);

% Also produce this for plotting.
dPrime_reBMFdp = makePopnTable(stats_stacked_4models,'rationalisedType', ...
    'usedamfreqs_dprime_rationalised','dprimes', ...
    coreoptions{:},'fn','Rayleigh>13.8','roworder',[1 2 4 5 6 3]);

% This is the mean dependence of normalised d' on modulation frequency.
% This weights the mean according the number of each type.
meanfModFn = nansum(dPrime_reBMFdp_reBestDP.n.*dPrime_reBMFdp_reBestDP.mean)./nansum(dPrime_reBMFdp_reBestDP.n,1);

% Now you need a vector of the correct value depending on fmod:
flist = unique(stats_stacked_4models.usedamfreqs_dprime_rationalised);
fModFactor = meanfModFn(arrayfun(@(x) find(flist==x),stats_stacked_4models.usedamfreqs_dprime_rationalised));


%% ----- Make the tables for the statistical models to fit to ----

% Sitmulus only.
nltables.stimonly = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', ...
           fModFactor', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','dprimes'});  
       
% For one additional continuous variable the table format is the same.
nltables.var1 = table(stats_stacked_4models.modLvl_rationalised',stats_stacked_4models.depthMod', ...
           fModFactor', stats_stacked_4models.VS', stats_stacked_4models.dprimes', ...
           'VariableNames',{'modLvl','depthMod','fFactor','Var1','dprimes'});  
        
% Type is a bit involved to make it work. Each type is a level in the
% factor. 
nltables.type = nltables.var1;   
nltables.type.Var1(1:length(statsmat_selected_4models.rationalisedType)) =   (   1*strcmp(statsmat_selected_4models.rationalisedType,'ChS') +  2*strcmp(statsmat_selected_4models.rationalisedType,'ChT') +  ...
                          3*strcmp(statsmat_selected_4models.rationalisedType,'On') +   4*strcmp(statsmat_selected_4models.rationalisedType,'PBU') +  ...
                          5*strcmp(statsmat_selected_4models.rationalisedType,'PL') +   6*strcmp(statsmat_selected_4models.rationalisedType,'PLN') )';

% Equation forms to use.             
% N.B. removed the highest mod-level and smallest depth from these, as they
% for three values there should only be 2 degrees of freedom (ANOVA style).
% Likewise there are 7 types but only 6 levels. 
% This has virtually no impact on r^2
nleqns.var1 = 'dprimes ~ fFactor*( v1*Var1 + l1*(modLvl==30) + l2*(modLvl==50) +  d2*(depthMod==1) + d3*(depthMod==2) )'; 
nleqns.stimonly = 'dprimes ~ fFactor*( l1*(modLvl==30) + l2*(modLvl==50)   + d2*(depthMod==1) + d3*(depthMod==2) )'; 
nleqns.type = ['dprimes ~ fFactor*( l1*(modLvl==30) + l2*(modLvl==50)  +' ...
                                'd2*(depthMod==1) + d3*(depthMod==2) +' ...
                                't1*(Var1==1) +  t2*(Var1==2) + ' ...
                                't3*(Var1==3) + t4*(Var1==4) +  t5*(Var1==5) + ' ...
                                't6*(Var1==6) )'];
                                            
                            
% Starting parameter fits to use. 
beta0.stimonly = [3 2 2 1.8];
beta0.type = [3 2 2 1.8 2 2 1 1 1 1];

%% -------------- Fitting the models ---------------

nlm.stimonly = fitnlm(nltables.stimonly,nleqns.stimonly,beta0.stimonly)
nlm.type = fitnlm(nltables.type,nleqns.type,beta0.type)

% An upper bound model where each dataset is separately fitted to the transfer function.  

% Take the stimonly fits and fit a single free parameter for each dataset.
nleqns.stimonly = 'dprimes ~ fFactor*( l1*(modLvl==30) + l2*(modLvl==50) + l3*(modLvl==70) + d1*(depthMod==0.5) + d2*(depthMod==1) + d3*(depthMod==2) )'; 

tmp.vars = nltables.stimonly(:,1:2)
tmp.est = nlm.stimonly.Coefficients.Estimate
tmp.fudge = tmp.est(3).*(tmp.vars.modLvl==30) + ...
            tmp.est(4).*(tmp.vars.modLvl==50)  + ...
            tmp.est(1).*(tmp.vars.depthMod==2) + ...
            tmp.est(2).*(tmp.vars.depthMod==3);
tmp.dprimecheck =  nltables.stimonly.fFactor .* tmp.fudge;  % Sums to zero if you added it correctly.
  
nleqns.stim_upperbound = 'dprimes ~ fFactor*(stimScalar + u1)';   % Equation incorportes the fitted stimonly model.
nltables.stim_upperbound = table( nltables.stimonly.fFactor, tmp.fudge, nltables.stimonly.dprimes, ...
            'VariableNames',{'fFactor','stimScalar','dprimes'});

% This fits the models to each dataset - so level and modulation depth
% separately. 
tmp.unitlist = unique([stats_stacked_4models.unitid' stats_stacked_4models.modLvl_rationalised'  stats_stacked_4models.depthMod'], 'rows' ); % Loop through each dataset.
tmp.coeff = nan(size(nltables.stim_upperbound,1),1);
tmp.Rsqrd = nan(size(nltables.stim_upperbound,1),1);
for ui = 1:length(tmp.unitlist)
    tmp.inds = find(stats_stacked_4models.unitid == tmp.unitlist(ui,1) & ...
                    stats_stacked_4models.modLvl_rationalised == tmp.unitlist(ui,2) & ...
                    stats_stacked_4models.depthMod == tmp.unitlist(ui,3) );
    tmp.subtable = nltables.stim_upperbound(tmp.inds,:);
    tmp.nlfit = fitnlm(tmp.subtable,nleqns.stim_upperbound,1);
    tmp.coeff(tmp.inds) = tmp.nlfit.Coefficients.Estimate;
    tmp.Rsqrd(tmp.inds) = tmp.nlfit.Rsquared.Ordinary;    
end;

% Now evaluate the function for each unit.
tmp.fitted = nltables.stim_upperbound.fFactor.*(nltables.stim_upperbound.stimScalar + tmp.coeff)
plot(nltables.stim_upperbound.dprimes,tmp.fitted,'.');
tmp.totalCorrCoef = corrcoef(nltables.stim_upperbound.dprimes,tmp.fitted);
tmp.overallR = tmp.totalCorrCoef(2);
nlm.stim_theoreticalmax = tmp;

clear tmp;
