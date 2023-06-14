function output = VSControlModel(pars,classifier)

% -------------------------------------------------------------------------

pars.srate = 10020;                       % Odd sample rate prevents artifacts in period histogram
phbins = 0:1/32:(1-1/32);


% This make a simple set of Poisson spike trains. Just like the expt. 
% This uses von Mises distribution as the basis for the time varying
% proability. Based on Kessler et al. 2021. 
for fi = 1:size(pars.amfreqs,1)    
    spiketimes(1:pars.nsweeps,fi) = simVS(pars.amfreqs(fi), ...
        pars.phaselocking(fi,1),pars.phaselocking(fi,2), ...
        pars.dcrate(fi),pars.srate,pars.dur_s,pars.nsweeps);
    
    % Implement a 0.8ms deadtime.
    %dt = diff(spiketimes{ri,fi});
    %spiketimes{ri,fi} = spiketimes{ri,fi}(find(dt>0.8e-3)+1);

    % Compute relative phases. N.B. spike times in milliseconds. 
    spikephase{fi} = rem([spiketimes{1:pars.nsweeps,fi}],1e3/pars.amfreqs(fi))/ ...
        (1e3/pars.amfreqs(fi));
    % Compute the period histogram.
    ph{fi} = hist(spikephase{fi},phbins);
    % Mean polar coords of phase. 
    polar(fi) = mean(cos( 2*pi*spikephase{fi} ) + i*sin(2*pi*spikephase{fi}));
    
end;

%spiketrainset{1} = spiketimes;
stimset = ones(size(spiketimes,1),1)*pars.amfreqs';

% WR classifier.
if nargin<2
    classifier.tau =1;
    classifier.iterations = 1000;
    classifier.store_distances = 0;
end;

% Run the model.
wr_test1 = WR_Classifier(spiketimes,stimset,classifier);

% Assign outputs.
output.wr = wr_test1;
output.ph = ph; 
output.vs(:,1) = abs(polar);
output.vs(:,2) = angle(polar);


% % Plot the results
% figure; 
% subplot(2,2,1);
% hbins = [0:0.25:80];
% hist([spiketimes{:,3}],hbins)
% subplot(2,2,2);
% plot(wr_test1.stimulusconditions,wr_test1.dprime);
% hold on; 
% plot(wr_test1.stimulusconditions,output.vs(:,1));
% subplot(2,2,3);
% imagesc(wr_test1.confusionmatrix);
% subplot(2,2,4);
% plot(phbins,ph{3});



