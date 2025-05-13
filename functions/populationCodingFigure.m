function populationCodingFigure(growingPopnByTypeLvl,options,populationCodingMeasure,maxscore,typeList,levels,scorestr)




%% -------------- Figure  ---------------

popnSize2Disp = 10; % The number of neurons for the population mtfs.
extraPopnSize2Disp = 30; % An additional one for the highest sound level.
                         % SHows how this differs between types. 
growXLim = [35 20 35 42];                         

typesToDisplay = {'ChS','ChT','PL','PLN'}
figh = figure('position',[100 100 1200 800],'paperposition',[.5 .5 24 16]);
typeColorMap;

for ci = 1:length(typesToDisplay)

    ti = find(strcmp(typeList, typesToDisplay{ci}));  % 1; % Type index.
    colorlist = min(1,[1.2 0.8 0.5]'*findRatTypeCol(typeList{ti}));  % Get the list of colours.

    mtf_h(ci) = subplot('position',[0.05+0.22*(ci-1) 0.63 0.18 0.27 ]);
    ngrow_h(ci) = subplot('position',[0.05+0.22*(ci-1) 0.2 0.18 0.27 ]);

    for li =1:3
        totalN = length(growingPopnByTypeLvl(li,ti).best1st.addedUnit);    
        totalNworst = length(growingPopnByTypeLvl(li,ti).worst1st.addedUnit);    

        % Some analysis of these results.     
        fullDataIndex_addedUnit = ...
            cellfun(@(x)  x.WRi{1}.fullDataMUIndexes, growingPopnByTypeLvl(li,ti).best1st.addedUnit(1:totalN));
        meanDPrime_addedUnit = ...
            cellfun(@(x) mean(x.WRi{1}.(populationCodingMeasure)(1:options.numFreqs4Mean)), growingPopnByTypeLvl(li,ti).best1st.addedUnit(1:totalN));
        meanDPrime_eachUnit = ...
            cellfun(@(x) mean(x.WRi{1}.(populationCodingMeasure)(1:options.numFreqs4Mean)),  ...
            growingPopnByTypeLvl(li,ti).best1st.eachUnit(  growingPopnByTypeLvl(li,ti).best1st.options.sampleOrder ));
        meanDPrime_addedUnit_rev = ...
            cellfun(@(x) mean(x.WRi{1}.(populationCodingMeasure)(1:options.numFreqs4Mean)), growingPopnByTypeLvl(li,ti).worst1st.addedUnit(1:totalNworst));
        meanDPrime_eachUnit_rev = ...
            cellfun(@(x) mean(x.WRi{1}.(populationCodingMeasure)(1:options.numFreqs4Mean)),  ...
            growingPopnByTypeLvl(li,ti).worst1st.eachUnit(  growingPopnByTypeLvl(li,ti).worst1st.options.sampleOrder ));
        fullDataMUIndex_lastUnit = ...
            cellfun(@(x) x.WRi{1}.fullDataMUIndexes(end), growingPopnByTypeLvl(li,ti).best1st.popns(1:totalN))
        nofFreqs_popn = ...
            cellfun(@(x)  length(x.WRi{1}.stimulusconditions), growingPopnByTypeLvl(li,ti).best1st.popns(1:totalN));
        maxFreq_popn = ...
            cellfun(@(x)  max(x.WRi{1}.stimulusconditions), growingPopnByTypeLvl(li,ti).best1st.popns(1:totalN));
        meanDPrime_popn = ...
            cellfun(@(x)  mean(x.WRi{1}.(populationCodingMeasure)(1:options.numFreqs4Mean)), growingPopnByTypeLvl(li,ti).best1st.popns(1:totalN));
        meanDPrime_popn_rev = ...
            cellfun(@(x)  mean(x.WRi{1}.(populationCodingMeasure)(1:options.numFreqs4Mean)), growingPopnByTypeLvl(li,ti).worst1st.popns(1:totalNworst));

        subplot(mtf_h(ci)); hold on;
        tmp_i = find(growingPopnByTypeLvl(li,ti).best1st.sampleN_bins(1:totalN)==popnSize2Disp);

        cellfun(@(x) plot( x.WRi{1}.stimulusconditions, x.WRi{1}.(populationCodingMeasure), 'color',colorlist(li,:),'linewidth',1), ...
                    growingPopnByTypeLvl(li,ti).best1st.popns(tmp_i))
        
        if li == 3
          tmp_ii = find(growingPopnByTypeLvl(li,ti).best1st.sampleN_bins(1:totalN)==extraPopnSize2Disp);
          if ~isempty(tmp_ii)
                cellfun(@(x) plot( x.WRi{1}.stimulusconditions, x.WRi{1}.(populationCodingMeasure),'--','color',colorlist(li,:),'linewidth',1), ...
                    growingPopnByTypeLvl(li,ti).best1st.popns(tmp_ii));
                        % PLot the sound level close to the line.
            stimf = growingPopnByTypeLvl(li,ti).best1st.popns{tmp_ii}.WRi{1}.stimulusconditions;                
            xperf = growingPopnByTypeLvl(li,ti).best1st.popns{tmp_ii}.WRi{1}.(populationCodingMeasure);
            text(stimf(li)-30,xperf(1+li),['n=' num2str(extraPopnSize2Disp)],'fontweight','bold','color',colorlist(li,:));

           end;
        end;
                
        ylim([-maxscore*0.03 maxscore]); xlim([50 950]); box off;
        xlabel('Modulation frequency (Hz)','fontweight','bold');
        %ylabel('d''','fontweight','bold');
        ylabel(['Ensemble classifier score ' scorestr],'fontweight','bold');
        if li == 1
            text(700,maxscore*0.95,typeList{ti},'fontweight','bold','fontsize',14);
        end;
        
        if ci == length(typesToDisplay) && li == 3 && ~isempty(strfind(scorestr,'(C)')) 
            % Add a "proportion correct" scale.
            pToSoftmaxZnoBiasFn = @(p,m) log(p) + log(m-1) - log(1-p)
            pvals = [0.2 0.6 0.9 0.98 0.997]
            pctickszvals = pToSoftmaxZnoBiasFn(pvals,9)
            line( [940 950]'*ones(1, length(pctickszvals)),[1 1]'*pctickszvals,'color','k'); 
            text(970*ones(1,length(pvals)),pctickszvals,num2str(pvals'),'fontsize',9);
             text(1200,6,'Equvialent hit rate','fontsize',9,'fontweight','bold','rotation',-90);

        end;
            
        % PLot the sound level close to the line.
        stimf = growingPopnByTypeLvl(li,ti).best1st.popns{tmp_i}.WRi{1}.stimulusconditions;                
        xperf = growingPopnByTypeLvl(li,ti).best1st.popns{tmp_i}.WRi{1}.(populationCodingMeasure);
        text(stimf(2+li)-30,xperf(2+li),[num2str(levels(li))],'fontweight','bold','color',colorlist(li,:));

        % Growth of performance with the the number of neurons. 
        subplot(ngrow_h(ci));      
        ngrow = growingPopnByTypeLvl(li,ti).best1st.sampleN_bins(1:totalN);
        plot(ngrow,meanDPrime_popn(1:totalN), 'color',colorlist(li,:) ,'linewidth',1); hold on; 
        plot(ngrow,meanDPrime_eachUnit(1:totalN),':', 'color',colorlist(li,:),'linewidth',1);

        % PLot the sound level close to the line.
        text(ngrow(5+2*li)-.5,meanDPrime_popn(5+2*li),[num2str(levels(li))],'fontweight','bold','color',colorlist(li,:));

        xlabel('Number of neurons','fontweight','bold');
        ylabel(['Mean ensemble score ' scorestr],'fontweight','bold');
        ylim([-maxscore*0.03 maxscore]); xlim([0 growXLim(ci)]); box off;
        
        if ci == length(typesToDisplay) && li == 3 && ~isempty(strfind(scorestr,'(C)')) 
            % Add a "proportion correct" scale.
            line( [growXLim(ci)-0.5 growXLim(ci)]'*ones(1, length(pctickszvals)),[1 1]'*pctickszvals,'color','k'); 
            text(1.02*growXLim(ci)*ones(1,length(pvals)),pctickszvals,num2str(pvals'),'fontsize',9);
            text(53,6,'Equvialent hit rate','fontsize',9,'fontweight','bold','rotation',-90);
        end;
        
        if li == 1
            text(7/9.5*max(xlim),maxscore*0.95,typeList{ti},'fontweight','bold','fontsize',14);
        end;
               
    end;        
end;

subplot(  mtf_h(1) );
text(0,maxscore*1.15,'a. Modulation transfer functions for small neural populations (n=10 or 30)','fontweight','bold','fontsize',10);
subplot(  ngrow_h(1) );
text(0,maxscore*1.15,'b. Average performance for f_m_o_ds up to 600Hz as a function of population size','fontweight','bold','fontsize',10);

%print -dtiff -r600 figures\Figure8.tif
%hgsave(figh,'figures\Figure8');