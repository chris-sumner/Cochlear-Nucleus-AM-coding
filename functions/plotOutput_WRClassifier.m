function plotOutput_WRClassifier(unitoutput,ahlist);

if nargin<2
    figure; set(gcf,'position',[100 100 500 400]);
    ahlist(1) = subplot(2,2,1);
    ahlist(2) = subplot(2,2,2);
    ahlist(3) = subplot(2,2,3);
end;

ploti = 1;

if isfield(unitoutput,'wr_spktrdist')
    subplot(ahlist(ploti));
    set(gca,'fontsize',8);
    imagesc(unitoutput.stimulusset(:),unitoutput.stimulusset(:),unitoutput.wr_spktrdist);
    xlabel('Stimulus value'); ylabel('Stimulus value');
    ploti=ploti+1;
end;    

subplot(ahlist(ploti));
set(gca,'fontsize',8);
imagesc(unitoutput.stimulusconditions,unitoutput.stimulusconditions,unitoutput.confusionmatrix);
xlabel('Stimulus value'); ylabel('Stimulus choice');
ploti=ploti+1;

subplot(ahlist(ploti));
set(gca,'fontsize',8);
plot(unitoutput.stimulusconditions,unitoutput.hitrate,'b'); hold on;
plot(unitoutput.stimulusconditions,unitoutput.falsealarmrate,'r');
%plot(unitoutput.stimulusconditions,unitoutput.dprime/max(unitoutput.dprime),'k');
chancerate = 1/length(unitoutput.stimulusconditions);
xlm = [min(unitoutput.stimulusconditions) max(unitoutput.stimulusconditions)];
line(xlm,[chancerate chancerate], ...
    'color','k','linestyle',':');
xlabel('Stimulus value'); %legend('HR','FAR','d''_N');
ylabel('Rate');
xlim(xlm); ylim([-.1 1.1]);


