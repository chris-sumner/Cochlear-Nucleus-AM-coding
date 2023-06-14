function chi2tests = plotMTFTypeProportions(statsmat,figname)
% Look at the proportion of bandpass units across type and level.


mtftypes = unique(statsmat.VS_MTFType);
unittypes = unique(statsmat.rationalisedType);
unittypes = {'ChS','ChT','PBU','PL','PLN','On'};

% VS table.
for ui = 1:length(unittypes)
    vsmtftype_unittype = statsmat.VS_MTFType(strcmp(statsmat.rationalisedType,unittypes{ui}));
 
   vs_nos(ui,:) =  ...
   [ sum(strcmp(vsmtftype_unittype,'lowpass')) ...
    sum(strcmp(vsmtftype_unittype,'bandpass')) ...
    sum(~strcmp(vsmtftype_unittype,'lowpass') & ~strcmp(vsmtftype_unittype,'bandpass'))];

end;


% VS table.
for ui = 1:length(unittypes)
    dpmtftype_unittype = statsmat.dprime_MTFType(strcmp(statsmat.rationalisedType,unittypes{ui}));
 
   dp_nos(ui,:) =  ...
   [ sum(strcmp(dpmtftype_unittype,'lowpass')) ...
    sum(strcmp(dpmtftype_unittype,'bandpass')) ...
    sum(~strcmp(dpmtftype_unittype,'lowpass') & ~strcmp(dpmtftype_unittype,'bandpass'))];

end;
maxv = max([sum(vs_nos,2); sum(dp_nos,2)]);

% % proportions:
fprintf('Proportions of lowpass MTF-VS:\n');
arrayfun(@(x,y) fprintf('%s  \t%g\n',x{1},y),unittypes,(100*vs_nos(:,1)./sum(vs_nos(:,1:2),2))')

fprintf('Proportions of lowpass MTF-d'':\n');
arrayfun(@(x,y) fprintf('%s  \t%g\n',x{1},y),unittypes,(100*dp_nos(:,1)./sum(dp_nos(:,1:2),2))')

% -------- Chi2 calculations ---------------

% Is there a difference of MTF-VS type across the units?
fprintf('\nChi-squared  test of proportions for MTF-VS type across neurons.\n');
chi2tests.VS.p_neuron = chi2test(vs_nos);

% Is there a difference of MTF-VS type across the units?
fprintf('\nChi-squared  test of proportions for MTF-d'' type across neurons.\n');
chi2tests.dp.p_neuron = chi2test(dp_nos);

% Type by type comparison of MTF shape proporitions across VS vs. d'
fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs VS, for ChS:\n');
chi2tests.VS_vs_dp.p_ChS = chi2test([vs_nos(1,1:2); dp_nos(1,1:2)]);
fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs VS, for ChT:\n');
chi2tests.VS_vs_dp.p_ChT = chi2test([vs_nos(2,1:2); dp_nos(2,1:2)]);
fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs VS, for PBU:\n');
chi2tests.VS_vs_dp.p_PBU= chi2test([vs_nos(3,1:2); dp_nos(3,1:2)]);
fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs VS, for PL:\n');
chi2tests.VS_vs_dp.p_PL = chi2test([vs_nos(4,1:2); dp_nos(4,1:2)]);
fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs VS, for PLN:\n');
chi2tests.VS_vs_dp.p_PLN = chi2test([vs_nos(5,1:2); dp_nos(5,1:2)]);
fprintf('\nChi-squared  test of proportions MTF-shape: d'' vs VS, for On:\n');
chi2tests.VS_vs_dp.p_On = chi2test([vs_nos(6,1:2); dp_nos(6,1:2)]);


% --------------------------------------------------
figure; 
subplot(2,1,1);
bar(vs_nos,'stacked'); box off;
ylabel('# datasets');
set(gca,'xticklabel',[]);
text(0,maxv*1.2,'A. MTF-VS classification (all levels etc.)','fontweight','bold');
ylim([ 0 maxv*1.1]);

subplot(2,1,2);
bar(dp_nos,'stacked');
lh = legend('Lowpass','Bandpass','Unc');
set(lh,'box','off');
box off;
ylabel('# datasets');
set(gca,'xticklabel',unittypes);
xlabel('Unit type');
text(0,maxv*1.2,'B. MTF-d'' classification (all levels etc.)','fontweight','bold');
ylim([ 0 maxv*1.1]);

print('-dtiff','-r300',[figname '_alllvls']);

% ----------------- SPLIT BY SOUND LEVEL ---------------

levels = [30 50 70];

% VS table.
coli = 1; vs_nos = []; dp_nos = [];
for ui = 1:length(unittypes)
    for li = 1:length(levels)
       vsmtftype_unittype = statsmat.VS_MTFType( ...
            strcmp(statsmat.rationalisedType,unittypes{ui}) & ...
            statsmat.modLevel>=(levels(li)-10) & statsmat.modLevel<(levels(li)+10));
 
       vs_nos(coli,:) =  ...
           [ sum(strcmp(vsmtftype_unittype,'lowpass')) ...
            sum(strcmp(vsmtftype_unittype,'bandpass')) ...
            sum(~strcmp(vsmtftype_unittype,'lowpass') & ~strcmp(vsmtftype_unittype,'bandpass'))];
        
      dpmtftype_unittype = statsmat.dprime_MTFType( ...
          strcmp(statsmat.rationalisedType,unittypes{ui}) & ...
          statsmat.modLevel>=(levels(li)-10) & statsmat.modLevel<(levels(li)+10));
 
       dp_nos(coli,:) =  ...
           [ sum(strcmp(dpmtftype_unittype,'lowpass')) ...
            sum(strcmp(dpmtftype_unittype,'bandpass')) ...
            sum(~strcmp(dpmtftype_unittype,'lowpass') & ~strcmp(dpmtftype_unittype,'bandpass'))];

        coli = coli+1;
    end;
    % create a gap between types.
    coli=coli+1;
end;


% -------- Chi2 calculations ---------------

% Is there a difference of MTF-VS type across the level?
fprintf('\nChi-squared  test of proportions for MTF-VS type across levels.\n');
chi2tests.VS.p_level = chi2test([  sum(vs_nos([1:4:end],:),1); ...
                                   sum(vs_nos([2:4:end],:),1); ...
                                   sum(vs_nos([3:4:end],:),1)  ]);

% Is there a difference of MTF-VS type across the levels?
fprintf('\nChi-squared  test of proportions for MTF-d'' type across levels.\n');
chi2tests.dp.p_level = chi2test([  sum(dp_nos([1:4:end],:),1); ...
                                   sum(dp_nos([2:4:end],:),1); ...
                                   sum(dp_nos([3:4:end],:),1)  ]);



maxv = max([sum(vs_nos,2); sum(dp_nos,2)]);

figure; 
subplot(2,1,1);
bar(vs_nos,'stacked'); box off;
ylabel('# datasets');
set(gca,'xticklabel',[],'xtick',[2 6 10 14 18 22]);
text(0,maxv*1.2,'A. MTF-VS classification','fontweight','bold');
ylim([0 maxv*1.1])

subplot(2,1,2);
bar(dp_nos,'stacked'); box off;
lh = legend('Lowpass','Bandpass','Unc');
set(lh,'box','off');
box off;
ylabel('# datasets');
set(gca,'xticklabel',unittypes,'xtick',[2 6 10 14 18 22]);
xlabel('Unit type');
text(0,maxv*1.2,'B. MTF-d'' classification ','fontweight','bold');
text(4.5,90,'30'); text(5.5,90,'50'); text(6.5,90,'70 dB SPL');
ylim([0 maxv*1.1])

print('-dtiff','-r300',[figname '_bylvl']);




function [p chistat] = chi2test(vs_nos)

% Get the numbers of observed of each type.
nLowpass = vs_nos(:,1);
nBandpass = vs_nos(:,2);
nEachType = nLowpass + nBandpass;

% pooled proportion of lowpass
p0 = sum(nLowpass)./sum( nEachType );
% expected numbers for each unit type.
nLowpass_expected = p0 * nEachType;
nBandpass_expected = (1-p0) * nEachType;

% Compile all expected and observed numbers.
observed = [nLowpass; nBandpass];
expected = [nLowpass_expected; nBandpass_expected];

% calculate the chi2 statistic.
chistat = sum( (observed - expected).^2 ./expected );

% Calculate the p value for 1 degree of freedom (2 MTF types).
p = 1- chi2cdf(chistat,1);

% Report this. 
fprintf('Chi-squared statistic: %g, p-value = %g\n', chistat, p);








