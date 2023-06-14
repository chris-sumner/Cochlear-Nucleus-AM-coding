function WREgFig(thesespikes,amf1,amf2,amfreqlist,alpha,ah)

greyshade = [.7 .7 .7];

% Find the indexes for the frequencies.
am1ind = find( amfreqlist == amf1 );
am2ind = find( amfreqlist == amf2 );

% Take 1 spike train from each frequency.
% If the same frequency - takes different spikes. 
st_am1 = thesespikes{am1ind}{1};
st_am2 = thesespikes{am2ind}{2};

% Compute the traces.
wrtraces.am1 = WRTrace(st_am1,st_am2,alpha,[20 110])

% Time series access to 
t = wrtraces.am1.Tlim(2)*[1:length(wrtraces.am1.diffSqTrace)]/length(wrtraces.am1.diffSqTrace);

if nargin<6
    figure; 
else 
    axes(ah);
end;

hold on;
plot(t,wrtraces.am1.yS,'color',greyshade,'linewidth',2)
plot(t,-wrtraces.am1.yT,'color',greyshade,'linewidth',2)
plot(t,wrtraces.am1.diffSqTrace,'color',[1 .3 .3],'linewidth',2);
plot(st_am1,0.8*ones(size(st_am1)),'o', ...
    'color','k','markerface',greyshade);
plot(st_am2,-0.8*ones(size(st_am2)),'o', ...
    'color','k','markerface',greyshade);
xlim([20 110])
ylim([-1 1]);

xlabel('Time (ms)','fontweight','bold');
