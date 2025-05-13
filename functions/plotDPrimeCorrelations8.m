function dPrimeCorrelations(stats_stacked,stackinds,scoremax,scorename,xlsfilename)

% ------------ Discrmination performance vs. other measures -----------

xfactor = 1.7;
xoffset = 0.2
x0 = xoffset +xfactor*.15; x1 = xoffset + xfactor*.21; x2 = xoffset + xfactor*.34; x3 = xoffset +.57; x4 = .79; 
w1 = 0.2; xfactor*.2; h1 = .27;
h0 = 0.1; w0 = 0.08; xfactor*0.07;
y0 = 0.17+.15; y1 = 0.4; 
y0 = 0.068; y1 = 0.18; y2 = 0.56;  y3 = 0.68;
yoffset = 0;

a1 = subplot('position',[x1 y3+yoffset w1 h1]);
[dh lh R optable] = scatterPlotStack(stats_stacked,'Z','dprimes','rationalisedType',stackinds);
axis([0 15 -.1 scoremax]); xlabel(''); ylabel('');
set(gca,'xticklabel',[],'yticklabel',[]);
text(0.5,scoremax*0.95,['R^2=' sprintf('%.2g',R(2).^2)]);
text(0,scoremax*1.1,'b. Z_I_S_I at f_m_o_d=150Hz','fontweight','bold');

% Make an Excel sheet with the data values in it. 
if nargin>4
    colnames = optable.Properties.VariableNames;
    xlswrite(xlsfilename,[colnames; table2cell(optable)],'Z_ISI');
end;

a2 = subplot('position',[x2 y3+yoffset w1 h1]);
[dh lh R optable] = scatterPlotStack(stats_stacked,'VS','dprimes','rationalisedType',stackinds);
axis([0 1 -.1 scoremax]); xlabel('');ylabel(''); 
set(gca,'xticklabel',[],'yticklabel',[]);
legend off;
text(0.07,scoremax*0.95,['R^2=' sprintf('%.2g',R(2).^2)]);
text(0,scoremax*1.1,'c. Vector Strength at f_m_o_d=150Hz','fontweight','bold');

% Make an Excel sheet with the data values in it. 
if nargin>4
    colnames = optable.Properties.VariableNames;
    xlswrite(xlsfilename,[colnames; table2cell(optable)],'VS');
end;

a3 = subplot('position',[x1 y1+yoffset w1 h1]);
[dh lh R optable] = scatterPlotStack(stats_stacked,'SACpeaks1234sum_150ish_p0001','dprimes','rationalisedType',stackinds,0);
axis([0 7 -.1 scoremax]); xlabel('');  ylabel('');
set(gca,'xticklabel',[],'yticklabel',[]);
legend off;
text(.25,scoremax*0.95,['R^2=' sprintf('%.2g',R(2).^2)]);
text(0,scoremax*1.1,'e. Sum of SAC peaks at f_m_o_d=150Hz','fontweight','bold');

% Make an Excel sheet with the data values in it. 
if nargin>4
    colnames = optable.Properties.VariableNames;
    xlswrite(xlsfilename,[colnames; table2cell(optable)],'Sum of SAC peaks');
end;


a4 = subplot('position',[x2 y1+yoffset w1 h1]);
[dh lh R optable] = scatterPlotStack(stats_stacked,'reliability_R_0p24ms','dprimes','rationalisedType',stackinds,0,true);
axis([-0.1 1 -.1 scoremax]); xlabel('');  ylabel('');
set(gca,'xticklabel',[],'yticklabel',[]);
legend off;
text(0.0,scoremax*0.95,['R^2=' sprintf('%.2g',R(2).^2) ' (log d'')']);
text(-0.1,scoremax*1.1,'f. PSTH reliabilty at f_m_o_d=150Hz','fontweight','bold');

% Make an Excel sheet with the data values in it. 
if nargin>4
    colnames = optable.Properties.VariableNames;
    xlswrite(xlsfilename,[colnames; table2cell(optable)],'Neural reliability');
end;

% Lower bar plots
subplot('position',[x1 y2+yoffset w1 h0]);
plotHorizontalSideBars(stats_stacked,'Z','rationalisedType',stackinds)
axis([0 15 -2.5 0]); box off;
xlabel('Z_I_S_I','fontweight','bold'); 

subplot('position',[x2 y2+yoffset w1 h0]);
plotHorizontalSideBars(stats_stacked,'VS','rationalisedType',stackinds)
axis([0 1 -2.5 0]); box off;
xlabel('VS','fontweight','bold'); 

subplot('position',[x1 y0+yoffset w1 h0]);
plotHorizontalSideBars(stats_stacked,'SACpeaks1234sum_150ish_p0001','rationalisedType',stackinds)
axis([0 7 -2.5 0]); box off;
xlabel('Sum of SAC peaks','fontweight','bold'); 

subplot('position',[x2 y0+yoffset w1 h0]);
plotHorizontalSideBars(stats_stacked,'reliability_R_0p24ms','rationalisedType',stackinds)
axis([-0.1 1 -2.5 0]); box off;
xlabel('Reliability','fontweight','bold'); 

% Side bar plot of d'.
subplot('position',[x0 y3+yoffset w0 h1]);
plotVerticalSideBars(stats_stacked,'dprimes','rationalisedType',stackinds)
axis([-3.5 0 -.1 scoremax]); 
ylabel(scorename,'fontweight','bold');

% Repeated
subplot('position',[x0 y1+yoffset w0 h1]);
plotVerticalSideBars(stats_stacked,'dprimes','rationalisedType',stackinds)
axis([-3.5 0 -.1 scoremax]); 
ylabel(scorename,'fontweight','bold');

