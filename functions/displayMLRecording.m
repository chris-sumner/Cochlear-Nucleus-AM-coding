function displayMLRecording(ii,ML,VS,CI,duration,fh)
% Displays the mode-locking analysis from a single recording.

figure(fh);

subplot(3,1,1); hold off;
subplot(3,1,2); hold off;
subplot(3,1,3); hold off;

phbins = [0:.02:1]; 
isibins = [0:0.5:25];
psthbins = [20:0.5:duration];

subplot(3,1,1);
if isempty(ML{ii})
    title(['Record #' num2str(ii)]);
    return;
end;
hist(ML{ii}.st,psthbins);
xlim([20 duration]);
hold on;
[t0,x0] = hist(ML{ii}.PHsur.st,psthbins);
plot(x0,t0,'r');
if ~isempty(ML{ii}.NHPP.fixedDt)
    [t0,x0] = hist(ML{ii}.NHPP.fixedDt.st,psthbins);
    plot(x0,t0,'m');
end;
xlabel('Time (ms)');
legend('Data','Phase shuffled','NHPP','ISI shuffled');
title(['Record #' num2str(ii) ' f_m_o_d = ' num2str(ML{ii}.modFreq) 'Hz' ...
        ' VS:' num2str(VS{ii}.VSvalue) ' CI:' num2str(CI{ii}.CIvalue_0lag)  ]);

subplot(3,1,2);
hist(ML{ii}.ph,phbins); title('Period histogram');
[t1,x1]  = hist(ML{ii}.PHsur.ph,phbins);
hold on; plot(x1,t1,'r');
[t1b,x1b]  = hist(ML{ii}.Isur.ph,phbins);
hold on; plot(x1b,t1b,'g');
if ~isempty(ML{ii}.NHPP.fixedDt)
    [t1c x1c] = hist(ML{ii}.NHPP.fixedDt.ph,phbins); title('PH');
    plot(x1c, t1c,'m'); 
end;
legend(['Data VS:' num2str(ML{ii}.stats.vals(1)) '/' num2str(ML{ii}.stats.vals(3))], ...
 'Phase shuffled', ...
['ISI shuffled VS:' num2str(ML{ii}.stats.vals(2)) '/' num2str(ML{ii}.stats.vals(4))], ...
 'NHPP ');
xlim([0 1]);

subplot(3,1,3);
hist(ML{ii}.ISI,isibins)
[t3 x3] = hist(ML{ii}.PHsur.ISI,isibins);
hold on; plot(x3, t3,'r'); 
[t3b x3b] = hist(ML{ii}.Isur.ISI,isibins);
hold on; plot(x3b, t3b,'g');
if ~isempty(ML{ii}.NHPP.fixedDt)
    [t3c x3c] = hist(ML{ii}.NHPP.fixedDt.ISI,isibins); 
    plot(x3c, t3c,'m'); 
    legstr =  ['NHPP Z:' num2str(ML{ii}.NHPP.fixedDt.Z)];
else
    legstr = '';
end;
title('ISI');
legend('Data', ...
 ['Phase shuffled Z:' num2str(ML{ii}.stats.vals(5))], ...
  'ISI shuffled ', ...
    legstr);
    xlim([0 25]);

    