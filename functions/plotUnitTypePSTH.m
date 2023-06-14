
% Clear all dimension variables.
clear   panelGap leftMargin bottomMargin  rastGap  bottomSideMargin rastX rastY rastAc ...
        phX  phY phAc  phGap isiX isiY  phY isiAc isiGap rastAx toneAx phAx isihAx ...
        raphLeftMargin thisUp graphSizeX graphSizeY graphYGap1  graphYGap2  imageSizeX  ...
        graphSizeX graphYGap3 dprimeX  dprimeY dprimeAc psthX psthY psthAc tmtfX tmtfY ...
        tmtfAc discUp discAc discX discY tmtfUp dprimeUp psthUp;

unit2runstr = num2str(unit2run);
%load(bngFn('rhode_modelocking\fromrhode\cn88data\allUnitEgs'),['unit' unit2runstr 'list'],['unit' unit2runstr 'spikes']);
load('rawdata\allUnitEgs',['unit' unit2runstr 'list'],['unit' unit2runstr 'spikes']);
thislist = eval(['unit' unit2runstr 'list']);
thesespikes = eval(['unit' unit2runstr 'spikes']);

% The unit number associated with these data.
this_eg_ind = find(unitids == unit2run); 

% A final specifier that Will make sure we get the right one - hand coded.
this_eg_ind = this_eg_ind(unitDataSetIndex);

% first collect everything that will need to be displayed
% if unit2run > 90000000
%     thisDir = bngFn('rhode_modelocking\fromrhode\cn91data');
% else
%     thisDir = bngFn('rhode_modelocking\fromrhode\cn88data');
% end
% thisDir = currentDir;
thisDir = 'rawdata';

expNo = unit2runstr(1:5);unitNo = unit2runstr(end-1:end);
eval(['load( ''' thisDir '\' num2str(expNo) 'unitStructs\Exp' num2str(expNo) 'U' num2str(unitNo) '.mat'' )']);
% eval(['load ' thisDir '\Exp' num2str(expNo) 'U' num2str(unitNo) '.mat']);
eval(['thisStruct = Exp' num2str(expNo) 'U' num2str(unitNo) 'struct;']);
%eval(['thisModStruct = Exp' num2str(expNo) 'U' num2str(unitNo) 'ModStruct;']);
%modStruct = thisModStruct(modStrChoice);

% PSTHs - normalised per second
psthLevs = thisStruct.psth.levels(psthChoice);
psthVals = thisStruct.psth.psthTenthms.values(psthChoice,:);
psthTimes = thisStruct.psth.psthTenthms.times;
psthTimes = psthTimes*1000; % in ms

colorgrey = [.4 .4 .4];
lightgrey = [.7 .7 .7];

% PSTH AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
%axes('position',[psthAc/xSize psthUp/ySize psthX/xSize psthY/ySize]);
bar(psthTimes,psthVals,'k');
maxY = round2(max(psthVals),100,'ceil');
set(gca,'xlim',[0 60],'ytick',[0 maxY],'yticklabel',[0 maxY],'ylim',[0 maxY],'FontSize',10);
%set(gca,'units','centimeters','box','off');
innPos = get(gca,'position');outPos = get(gca,'outerposition');furthLeft = outPos(1) - innPos(1);
xlabel('Time (ms)','fontweight','bold');yh=ylabel('Spikes/s','fontweight','bold');
box off;
%set(yh,'units','centimeters');
% if length(freqConds)==4
%     oldPos = get(yh,'position');
%     set(yh,'position',[furthLeft-.02 oldPos(2:3)])
% end;



