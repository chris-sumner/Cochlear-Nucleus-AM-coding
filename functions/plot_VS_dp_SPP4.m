function po = plot_VS_dp_SPP4(stats_stacked,row,mod1opt,coreoptions,plotopts)
%printtitle,roworder)
% VS,dp & SPP vs. AM re: BMF as figure 4.

% ------------- process plotting options ---------

po.printtitle = [];
po.roworder = [1 2 4 5 6 3];
po.rowtitle = [];
po.figh = [];
po.ahs = [];
po.axes_y = [0.15 0.75];    % Normal y coord method. 
po.axes_x = [0.2 0.14];     % 1st number is the total gap between panels.
po.colorlist = {'m','g','r','k','b','c'};
po.axeslist = [];
po.offset_x = .05;
po.sem = 'sem';

i=1;
while i<length(plotopts)
    po.(plotopts{i}) = plotopts{i+1};
    i=i+2;
end;



% ---------- set up the figure ------------

if isempty(po.figh)
    po.figh  = figure('position',[100 100 1800 300],'paperposition',[.5 .5 36 6]);
end;
%figure(po.figh);
set(po.figh,'DefaultAxesFontSize',8); 



% --------- Plot the tables ------------

axesfn = @(ai)  [po.offset_x+po.axes_x(1)*(ai-1) po.axes_y(1) po.axes_x(2) po.axes_y(2) ];


for ai = 1:length(po.axeslist)
    po.ahs(ai) = subplot('position',axesfn(ai));
   
    switch (po.axeslist{ai})
        case 'VS'
            % VS - make the table
            rowandcol =  {row,'usedamfreqs_reBMF_rationalised'};
            stats.VS_reBMF = makePopnTable(stats_stacked,rowandcol{:},'VS',coreoptions{:}, ...
                mod1opt{:},'fn','Rayleigh>13.8','roworder',po.roworder);
            % N.B. Could linearly interpolate VS measurements to resmaple AM freq.
            % better. 
            
            % Plot
            plotTable2(stats.VS_reBMF,po.sem,'colorlist',po.colorlist); 
            ylim([0.1 .9]); xlim([-500 2000]);
            xlabel('f_m_o_d re:BMF-VS (Hz)','fontweight','bold');
            ylabel('VS','fontweight','bold');
        case 'DP'
            % D-PRIME - normalise to the best dprime (not BMF)
            rowandcol =  {row,'usedamfreqs_reBMFdp_rationalised'};
            stats.dPrime_reBMFdp = makePopnTable(stats_stacked,rowandcol{:},'dprimes',coreoptions{:} ...
                ,mod1opt{:},'fn','Rayleigh>13.8','roworder',po.roworder);
            %stats.dPrime_reBMF = makePopnTable(stats_stacked,rowandcol{:},'dprimes',coreoptions{:},mod1opt{:},'fn','Rayleigh>13.8');
            
            % Plot
            plotTable2(stats.dPrime_reBMFdp,po.sem,'colorlist',po.colorlist);
            ylim([0 6]); xlim([-500 2000]);
            xlabel('f_m_o_d re:BMF (Hz)','fontweight','bold');
            ylabel('d''','fontweight','bold');
        case 'DP_refMod'
            % D-PRIME - normalise to the best dprime (not BMF)
            rowandcol =  {row,'usedamfreqs_rationalised'};
            stats.dPrime_reBMFdp = makePopnTable(stats_stacked,rowandcol{:},'dprimes',coreoptions{:} ...
                ,mod1opt{:},'fn','Rayleigh>13.8','roworder',po.roworder);
            %stats.dPrime_reBMF = makePopnTable(stats_stacked,rowandcol{:},'dprimes',coreoptions{:},mod1opt{:},'fn','Rayleigh>13.8');
            
            % Plot
            plotTable2(stats.dPrime_reBMFdp,po.sem,'colorlist',po.colorlist);
            ylim([0 5]); xlim([0 2000]);
            xlabel('f_m_o_d (Hz)','fontweight','bold');
            ylabel('d''','fontweight','bold');
        case 'SPP'
            % Spikes per period. make the table 
            rowandcol =  {row,'usedamfreqs_reBMF_rationalised'};
            stats.SPP_reBMF = makePopnTable(stats_stacked,rowandcol{:},'spikesperperiod', ...
                coreoptions{:},mod1opt{:},'fn','Rayleigh>13.8','roworder',po.roworder);
            
            % Plot.
            plotTable2(stats.SPP_reBMF,'colorlist',po.colorlist);
            ylim([0 6]);
            xlabel('f_m_o_d re:BMF-VS (Hz)','fontweight','bold');
        case 'SPP_reBMFdp'
            % SPP but relative to BMFdp - make the table. 
            rowandcol =  {row,'usedamfreqs_reBMFdp_rationalised'};
            stats.SPP_reBMFdp = makePopnTable(stats_stacked,rowandcol{:},'spikesperperiod', ...
                coreoptions{:},mod1opt{:},'fn','Rayleigh>13.8','roworder',po.roworder);

            plotTable2(stats.SPP_reBMFdp,'colorlist',po.colorlist);
            ylim([0 6]);
            xlabel('f_m_o_d re:BMF-d'' (Hz)','fontweight','bold');
        case 'VSn'
            % Normalise the VS to the peak VS 
            rowandcol =  {row,'usedamfreqs_reBMF_rationalised'};
            stats.VS_reBMF_reBestVS = makePopnTable(stats_stacked,rowandcol{:}, ...
                'VS_pcBestVS',coreoptions{:},mod1opt{:},'fn','Rayleigh>13.8', ...
                'roworder',po.roworder);

            % Plot.
            plotTable2(stats.VS_reBMF_reBestVS,'colorlist',po.colorlist);
            ylim([0.1 1.1]); xlim([-1000 1500]);
            xlabel('f_m_o_d re:BMF-VS (Hz)','fontweight','bold');
        case 'DPn'
            % D prime this time normalised to the best dp - so difference from peak. 
            rowandcol =  {row,'usedamfreqs_reBMFdp_rationalised'};
            stats.dPrime_reBMFdp_reBestDP = makePopnTable(stats_stacked,rowandcol{:},'dprimes_pcBestDP', ...
                coreoptions{:},mod1opt{:},'fn','Rayleigh>13.8','roworder',po.roworder);

            % Plot.    
            plotTable2(stats.dPrime_reBMFdp_reBestDP,'colorlist',po.colorlist);
            ylim([0 1.1]); xlim([-500 2000]);
            xlabel('f_m_o_d re:BMF-d'' (Hz)','fontweight','bold');
        case 'DPn_refMod'
            % D prime this time normalised to the best dp - so difference from peak. 
            rowandcol =  {row,'usedamfreqs_rationalised'};
            stats.dPrime_reBMFdp_reBestDP = makePopnTable(stats_stacked,rowandcol{:},'dprimes_pcBestDP', ...
                coreoptions{:},mod1opt{:},'fn','Rayleigh>13.8','roworder',po.roworder);

            % Plot.    
            plotTable2(stats.dPrime_reBMFdp_reBestDP,'colorlist',po.colorlist);
            ylim([0 1.1]); xlim([0 2000]);
            xlabel('f_m_o_d (Hz)','fontweight','bold');
        case 'SACpeaks_reBMFdp'
             % D prime this time normalised to the best dp - so difference from peak. 
            rowandcol =  {row,'usedamfreqs_reBMFdp_rationalised'};
            stats.SACpeaks_reBMFdp_reBestDP = makePopnTable(stats_stacked,rowandcol{:},'numberofpeaks', ...
                coreoptions{:},mod1opt{:},'fn','CI_p>.99','roworder',po.roworder);

            % Plot.    
            plotTable2(stats.SACpeaks_reBMFdp_reBestDP,'colorlist',po.colorlist);
            ylim([0 15]); xlim([-500 2000]);
            xlabel('f_m_o_d re:BMF-d'' (Hz)','fontweight','bold');
            
            
        case 'Z'
            rowandcol =  {row,'usedamfreqs_rationalised'};         
            stats.Z = makePopnTable(stats_stacked,rowandcol{:},'Z', ...
                coreoptions{:},mod1opt{:},'fn','Rayleigh>13.8','fn','VS>.1','roworder',po.roworder);
            
           %      Plot.    
            plotTable2(stats.Z,'colorlist',po.colorlist);
            %ylim([0 15]); xlim([-500 2000]);
            %xlabel('f_m_o_d re:BMF-d'' (Hz)','fontweight','bold');
        
         case 'CI'
           rowandcol =  {row,'usedamfreqs_reBMFdp_rationalised'};         
            stats.CI = makePopnTable(stats_stacked,rowandcol{:},'CI', ...
                coreoptions{:},mod1opt{:},'fn','CI_p>.99','roworder',po.roworder);
            
           %      Plot.    
            plotTable2(stats.CI,'colorlist',po.colorlist);
            %ylim([0 15]); xlim([-500 2000]);
            %xlabel('f_m_o_d re:BMF-d'' (Hz)','fontweight','bold');
            

    end;            
    
    if ai==1 & ~isempty(po.rowtitle)
        text(min(xlim),diff(ylim)*1.1 + min(ylim),po.rowtitle,'fontweight','bold');
    end;
    
    if ai<length(po.axeslist)
        legend off;
    end;

end;


% Option to print the figure.
if ~isempty(po.printtitle)
    print('-dtiff','-r150',po.printtitle);
end;



