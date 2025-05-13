function variables = makeFigure9AS6S8DataSheet(fname,stats_stacked_4models,nlm,option)
% makeFigure9AS6S8DataSheet  writes all the raw values for correlations to
%                            a spreadsheet. 

% ----------- A sheet with all the variables used in the models -------

% Reconstruct the normalised gain for each data point. 
fFactor = nan(size(stats_stacked_4models.unitid'));
% Reconstruct the frequency function value for each data point. Fixed for
% all models.
[lgc f_ind] = ismember(stats_stacked_4models.allamfreqs_rationalised,nlm.stim_theoreticalmax.flist);
fFactor = nlm.stim_theoreticalmax.normalised_fitted_fn(f_ind);

figure9sheet = table( ...
  (stats_stacked_4models.unitid'), (stats_stacked_4models.modLevel'), ...
  (stats_stacked_4models.depthMod'), (stats_stacked_4models.usedamfreqs_rationalised'), ...
   stats_stacked_4models.rationalisedType', (stats_stacked_4models.modLvl_rationalised'), ...
  (stats_stacked_4models.usedamfreqs_rationalised'),  fFactor', ...
  stats_stacked_4models.VS_150ish',   stats_stacked_4models.CV' , ...
  stats_stacked_4models.Z_150ish',   stats_stacked_4models.reliability_R150ish_0p24ms' , ...
  stats_stacked_4models.envFluct_R150ish_0p48ms',  stats_stacked_4models.CI_150ish', ...
  stats_stacked_4models.SACpeak0sal_150ish_p0001',   ...
  stats_stacked_4models.SACpeak1sal_150ish_p0001',   ...    
 stats_stacked_4models.SACpeak2sal_150ish_p0001',    ...     
 stats_stacked_4models.SACpeak3sal_150ish_p0001',    ...      
 stats_stacked_4models.SACpeak4sal_150ish_p0001',    ...     
 stats_stacked_4models.SACpeak5sal_150ish_p0001', ...
  stats_stacked_4models.VS',   ...
  stats_stacked_4models.Z',   stats_stacked_4models.reliability_R_0p24ms' , ...
  stats_stacked_4models.envFluct_0p48ms',  stats_stacked_4models.CI', ...
  stats_stacked_4models.SACpeak0sal_p0001',   ...
  stats_stacked_4models.SACpeak1sal_p0001',   ...    
 stats_stacked_4models.SACpeak2sal_p0001',    ...     
 stats_stacked_4models.SACpeak3sal_p0001',    ...      
 stats_stacked_4models.SACpeak4sal_p0001',    ...     
 stats_stacked_4models.SACpeak5sal_p0001', ...
 stats_stacked_4models.dprimes');     

colnames = {'NeuronID','SoundLevel','ModulationDepth', 'ModulationFreq', ...
            'Type','SoundLevelRationalised','ModulationFreqRationalised', ...
            'fModFactor','VS_150Hz','CV','Z_150Hz','NeuralReliability_150Hz', ...
            'EnvelopeFluctuation_150Hz','CI_150Hz', ...
            'SAC_peak0_150Hz','SAC_peak1_150Hz','SAC_peak2_150Hz', ...
            'SAC_peak3_150Hz','SAC_peak4_150Hz','SAC_peak5_150Hz', ...
            'VS','Z','NeuralReliability', ...
            'EnvelopeFluctuation','CI', ...
            'SAC_peak0','SAC_peak1','SAC_peak2', ...
            'SAC_peak3','SAC_peak4','SAC_peak5', ...
            'C_prime'};


if strcmp(option,'9A')
    % ---------- All the prediction parameters -----------
    % Only store these in the Figure 9A datasheet.
    
    figure9sheet = renamevars(figure9sheet, figure9sheet.Properties.VariableNames, colnames);              
    xlswrite(fname,[colnames; table2cell(figure9sheet)],'Variables');
end;
        
if strcmp(option,'9A') || strcmp(option,'S6') 

    % ---------- Sheets for each of the models in 9A -----------

    % Predictors in clear. 
    writeSheet(nlm.stimonly,fname,'stimuliOnly','C_prime','Predicted_C_prime');

    % Predictors in black.
    writeSheet(nlm.type,fname,'addType','C_prime','Predicted_C_prime');
    writeSheet(nlm.CV,fname,'addCV','C_prime','Predicted_C_prime');

    % Predictors in yellow
    writeSheet(nlm.VS_150ish,fname,'addVS_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak1Sal_150ish_p0001,fname,'addSACpeak1_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak2345Sal_150ish_p0001,fname,'addSACpeak2to5_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.Z_150ish,fname,'addZ_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak12345Sal_150ish_p0001,fname,'addSACpeak1to5_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.reliability_150ish,fname,'addReliability_150Hz','C_prime','Predicted_C_prime');

    % Green in Figure 9
    writeSheet(nlm.TwoVars_PeaksPlusZ_150ish,fname,'addZandSACpeaks_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.TwoVars_ZPlusReliability_150ish,fname,'addZandReliability_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.TwoVars_PeaksPlusReliability_150ish,fname,'addSACPeaksandReliability_150Hz','C_prime','Predicted_C_prime');

    % Blue in Figure 9
    writeSheet(nlm.AllSensibleVars_minusNuFl_150ish,fname,'add3best_150Hz','C_prime','Predicted_C_prime');

    % Purple in Figure 9
    writeSheet(nlm.Z_150ish_noChS,fname,'addZ_150Hz_noChS','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak12345Sal_150ish_p0001_noChS,fname,'addSACpeak1to5_150Hz_noChS','C_prime','Predicted_C_prime');
    writeSheet(nlm.reliability_150ish_noChS,fname,'addReliab_150Hz_noChS','C_prime','Predicted_C_prime');

end;

if strcmp(option,'S6') 

    % Additional yellow in S6
    writeSheet(nlm.envFluct_150ish,fname,'addNeuralFluct_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.CI_150ish,fname,'addCI_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.numSACpeaks_150ish,fname,'addNumSACpeaks_150Hz','C_prime','Predicted_C_prime');


    % Purple in Figure S6
    writeSheet(nlm.VS,fname,'addVS','C_prime','Predicted_C_prime');
    writeSheet(nlm.CI,fname,'addCI','C_prime','Predicted_C_prime');
    writeSheet(nlm.numSACpeaks,fname,'addNumSACpeaks','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeaks1Sal_p0001,fname,'addSACpeak1','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeaks2345Sal_p0001,fname,'addSACpeak2to5','C_prime','Predicted_C_prime');
    writeSheet(nlm.Z,fname,'addZ','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeaks12345Sal_p0001,fname,'addSACpeak1to5','C_prime','Predicted_C_prime');
    writeSheet(nlm.reliability,fname,'addNeuralReliability','C_prime','Predicted_C_prime');
    writeSheet(nlm.envFluct,fname,'addNeuralFluctuation','C_prime','Predicted_C_prime');

    % Blue Figure S6
    writeSheet(nlm.AllSensibleVars_minusNuFl_150ish,fname,'add3_notNeuralFluct_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.AllSensibleVars_minusZ_150ish,fname,'add3_notZ_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.AllSensibleVars_minusSACPeaks_150ish,fname,'add3_notSACpeaks_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.AllSensibleVars_minusReliability_150ish,fname,'add3_notReliability_150Hz','C_prime','Predicted_C_prime');

    % Additional green in S6
    writeSheet(nlm.TwoVars_PeaksPlusNuFl_150ish,fname,'addNeuralFlandSACpeaks_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.TwoVars_NuFlPlusReliability_150ish,fname,'addNeuralFlandReliability_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.TwoVars_ZPlusNuFl_150ish,fname,'addNeuralFlandZ_150Hz','C_prime','Predicted_C_prime');

    % Red in S6
    writeSheet(nlm.AllSensibleVars_150ish,fname,'add4best_150Hz','C_prime','Predicted_C_prime');

end;


if strcmp(option,'S8') 

    % Additional SAC columns for S8.
    % Starting at zero
    writeSheet(nlm.SACpeak0Sal_150ish_p0001,fname,'addSACpeak0_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeaks01Sal_150ish_p0001,fname,'addSACpeak01_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeaks012Sal_150ish_p0001,fname,'addSACpeak012_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeaks0123Sal_150ish_p0001,fname,'addSACpeak0123_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeaks01234Sal_150ish_p0001,fname,'addSACpeak01234_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeakAllSals_150ish,fname,'addSACpeak012345_150Hz','C_prime','Predicted_C_prime');
    % Starting at 1
    writeSheet(nlm.SACpeak1Sal_150ish_p0001,fname,'addSACpeak1_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak12Sal_150ish_p0001,fname,'addSACpeak12_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak123Sal_150ish_p0001,fname,'addSACpeak123_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak1234Sal_150ish_p0001,fname,'addSACpeak1234_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak12345Sal_150ish_p0001,fname,'addSACpeak12345_150Hz','C_prime','Predicted_C_prime');
    % Omitting 1
    writeSheet(nlm.SACpeak0Plus2Sal_150ish_p0001,fname,'addSACpeak02_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak0Plus23Sal_150ish_p0001,fname,'addSACpeak023_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak0Plus234Sal_150ish_p0001,fname,'addSACpeak0234_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak0Plus2345Sal_150ish_p0001,fname,'addSACpeak02345_150Hz','C_prime','Predicted_C_prime');
    % Starting at 2
    writeSheet(nlm.SACpeak23Sal_150ish_p0001,fname,'addSACpeak23_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak234Sal_150ish_p0001,fname,'addSACpeak234_150Hz','C_prime','Predicted_C_prime');
    writeSheet(nlm.SACpeak2345Sal_150ish_p0001,fname,'addSACpeak2345_150Hz','C_prime','Predicted_C_prime');

end;

variables = figure9sheet;

function writeSheet(m,fname,sheetname,dataname,predname)

data4sheet = num2cell([m.Variables.dprimes m.Fitted]);
data4sheet = [ {dataname},{predname}; data4sheet ];
data4sheet{1,3} = 'R^2';
data4sheet{2,3} = corr(m.Variables.dprimes,m.Fitted)^2;

xlswrite(fname,data4sheet,sheetname);



