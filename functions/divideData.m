function [newdatainds dividingdata] = divideData(data,EXPLOGLIST,AMtype)
% divideData  - mark out where individual MTFs start. 
%             N.B. assumption is that this also marks the end of the
%             previous.

% Which columns to use to divide the data.
cols2dividedata = [  find(strcmp(EXPLOGLIST,'uniqueID')) ...
                     find(strcmp(EXPLOGLIST,'modLevel')) ...
                     find(strcmp(EXPLOGLIST,'amToneDur')) ...
                     find(strcmp(EXPLOGLIST,'carrFreq')) ...
                     find(strcmp(EXPLOGLIST,'depthMod')) ...
                     find(strcmp(EXPLOGLIST,'modFreq'))];
% Find how each changes. 
dividingdata = data(:,cols2dividedata);
% Any nans must be replaced by -1 for this to work with them. 
dividingdata(isnan(dividingdata)) = -1;

diffdivdata = [ones(1,length(cols2dividedata)-1) -1; diff(dividingdata)];

% Also look for changes in AM type.
amdiff(1) = 0;
for i=2:length(AMtype)
    amdiff(i) = strcmp(AMtype{i},AMtype{i}-1);
end;
                 
% Now use all these to demarkate the beginnings of the next MTF. 
newdatainds = find( diffdivdata(:,1)~=0 |  diffdivdata(:,2)~=0 | ...
                    diffdivdata(:,3)~=0 |  diffdivdata(:,4)~=0 | ...
                    diffdivdata(:,5)~=0 |  diffdivdata(:,6)<0 | ...
                    amdiff'~=0 ); 
                % Indexes to the first row of individual mtfs.

                