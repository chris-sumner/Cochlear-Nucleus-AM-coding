function DPvsVSegs_paperfig_v2(statsmat_oneperunit)

% Hand selected examples for each type. 
examples.ChS_unitid = [88328009   88335020    88355025    88399028    89001030    ]
examples.ChT_unitid = [88328049    88374013    88399003    88401005    89001011 ];
examples.PL_unitid = [88349020    88349055    88393079    88393105    88393108];
examples.PLN_unitid = [88299042    88310023    88340046    90275121    90340011 ];

examples.ChS = find(ismember([statsmat_oneperunit.unitid{:}],examples.ChS_unitid));
examples.ChT = find(ismember([statsmat_oneperunit.unitid{:}],examples.ChT_unitid));
examples.PL = find(ismember([statsmat_oneperunit.unitid{:}],examples.PL_unitid));
examples.PLN = find(ismember([statsmat_oneperunit.unitid{:}],examples.PLN_unitid))

figh = figure('position',[100 100 1200 800],'paperposition',[.5 .5 24 16]);

psthsize = [ 0.19 0.14];
mtfsize = [0.19,0.2];

%% ----------------- Plot PSTHS -----------------------
global bngFn;
load('datafiles\WR_results','unitoutputs');
% Extract unit outputs and a unitid (indentifier with which to find the
% unit).
unitinfo = [unitoutputs(:).unitinfo];
unitids = [unitinfo(:).unitid];


% ChS ***************************************************************
% 88299021, 90275099, 91057069, 91060018
% % CARR F = 24000, MOD LEVEL = 30, THR = 16
unitType = 'ChS';
unit2run = 88299021;
unitDataSetIndex = 3;
psthChoice = 2;

axes('position',[0.05 0.77 psthsize]);
plotUnitTypePSTH;
text(-5,6000,'a. Example pure tone responses','fontweight','bold','fontsize',12);
text(50,4900*.9,'ChS','fontweight','bold','fontsize',14);
ylp = get(yh,'Position');
ylp(1) = ylp(1) + 5;
set(yh,'Position',ylp)

%% ChT ***************************************************************
% 91016074
% % CARR F = 12000, MOD LEVEL = 50, THR = 18
unitType = 'ChT';
unit2run = 91016074;
% levCondInds = 27:52; % 50 dB
unitDataSetIndex =1;
psthChoice = 1;%modStrChoice=1;% 30 dB PSTH

axes('position',[0.29 0.77 psthsize]);
plotUnitTypePSTH;
ylabel('');
text(50,2000*.9,'ChT','fontweight','bold','fontsize',14);

% PL ****************************************************************
unitType = 'PL';
unit2run = 88340053;
unitDataSetIndex = 1;
% 60 dB PSTH
psthChoice = 2;
axes('position',[0.53 0.77 psthsize]);
plotUnitTypePSTH;
ylabel('');
text(50,900*.9,'PL','fontweight','bold','fontsize',14);

% PLN ***************************************************************
% 91019022, 91057055
% % CARR F = 25500, MOD LEVEL = 50, THR = 29
% THIS DID NOT WORK WHEN USIING THE 30dB RESPONSES - no responses.
unitType = 'PLN';
unit2run = 91019022;
unitDataSetIndex =1;
psthChoice = 2; %modStrChoice=4;% 60 dB PSTH
axes('position',[0.77 0.77 psthsize]);
plotUnitTypePSTH;
ylabel('');
text(50,1400*.9,'PLN','fontweight','bold','fontsize',14);


%% MTFs

dph_ChS = subplot('position',[0.05,0.45,mtfsize]); %2,4,1);
vsh_ChS = subplot('position',[0.05,0.13,mtfsize]);  %2,4,5);
plotMTFsOneType(statsmat_oneperunit,'ChS',examples.ChS,[dph_ChS vsh_ChS]);
subplot(dph_ChS); axis([0 800 -0.1 7]); box off; ylabel('d prime','fontweight','bold');
text(-100,8,'b. Examples of classifier performance (MTF-d'')','fontweight','bold','fontsize',12);
set(dph_ChS,'xtick',[0:200:800],'ytick',[1:7],'fontsize',10); title('');

subplot(vsh_ChS); axis([0 800 -0.02 1]); box off; ylabel('VS','fontweight','bold');
text(-100,1.2,'c. Phase-locking in the same examples (MTF-VS)','fontweight','bold','fontsize',12);
set(vsh_ChS,'xtick',[0:200:800],'ytick',[0:.2:1],'fontsize',10); title(''); xlabel('f_m_o_d(Hz)','fontweight','bold');

dph_ChT = subplot('position',[0.29,0.45,mtfsize]); %2,4,1);
vsh_ChT = subplot('position',[0.29,0.13,mtfsize]);  %2,4,5);
plotMTFsOneType(statsmat_oneperunit,'ChT',examples.ChT,[dph_ChT vsh_ChT]);
subplot(dph_ChT); 
axis([0 800 -0.1 7]); box off;
set(dph_ChT,'xtick',[0:200:800],'ytick',[1:7],'fontsize',10); title('');
subplot(vsh_ChT); axis([0 800 -0.02 1]); box off;
set(vsh_ChT,'xtick',[0:200:800],'ytick',[0:.2:1],'fontsize',10); title(''); xlabel('f_m_o_d(Hz)','fontweight','bold');

dph_PL = subplot('position',[0.53,0.45,mtfsize]); %2,4,1);
vsh_PL = subplot('position',[0.53,0.13,mtfsize]);  %2,4,5);
plotMTFsOneType(statsmat_oneperunit,'PL',examples.PL,[dph_PL vsh_PL]);
subplot(dph_PL); 
axis([0 1800 -0.1 7]); box off; 
set(dph_PL,'xtick',[0:400:2000],'ytick',[1:7],'fontsize',10); title('');
subplot(vsh_PL); 
axis([0 1800 -0.02 1]); box off;
set(vsh_PL,'xtick',[0:400:2000],'ytick',[0:.2:1],'fontsize',10); title(''); xlabel('f_m_o_d(Hz)','fontweight','bold');

dph_PLN = subplot('position',[0.77,0.45,mtfsize]); %2,4,1);
vsh_PLN = subplot('position',[0.77,0.13,mtfsize]);  %2,4,5);
plotMTFsOneType(statsmat_oneperunit,'PLN',examples.PLN,[dph_PLN vsh_PLN]);
subplot(dph_PLN); 
axis([0 1800 -0.1 7]); box off; 
set(dph_PLN,'xtick',[0:400:2000],'ytick',[1:7],'fontsize',10); title('');
subplot(vsh_PLN); 
axis([0 1800 -0.02 1]); box off;
set(vsh_PLN,'xtick',[0:400:2000],'ytick',[0:.2:1],'fontsize',10); title(''); xlabel('f_m_o_d(Hz)','fontweight','bold');
