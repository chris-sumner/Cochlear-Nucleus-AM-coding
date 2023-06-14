function popntable = makePopnTable(stackdata,rowname,colname,valuename, varargin)
% makePopnTables  make tables of all statistics.  

% Find the table dimensions and values 
rowvalues = unique(stackdata.(rowname));  % All the column values.
if isnumeric(rowvalues)
    rowvalues = rowvalues( ~isnan(rowvalues) );% Remove nans.
end;
if iscell(rowvalues) && isstr(rowvalues{1})
    rowvalues = rowvalues(cellfun(@(x) ~isempty(x),rowvalues)); %remove empty strings.
end;

colvalues = unique(stackdata.(colname));  % All the row values. 
if isnumeric(colvalues)
    colvalues = colvalues( ~isnan(colvalues) ); % remove nans
end;
if iscell(colvalues) && isstr(colvalues{1})
    colvalues = colvalues(cellfun(@(x) ~isempty(x),colvalues)); %remove empty strings.
end;


% Store these. 
popntable.rowName = rowname;
popntable.rowValues = rowvalues;
popntable.colName = colname;
popntable.colValues = colvalues;
popntable.valueName = valuename;

% ------------ Defaults ----------------

% Don't give any statistics.
norm_flag = false;
mean_flag = false;
sd_flag = false;
hist_flag = false;
field_fns = {};
n_field_fns = 0;

% options.
min_n = 1;


% Initialise 'n' (always do this). 
popntable.n = nan(length(rowvalues),length(colvalues));


ind = find(strcmp('normalise',varargin));      % Normalise histogram?
if ~isempty(ind)
    norm_flag = true;
end;

ind = find(strcmp('mean',varargin));           % Compute mean
if ~isempty(ind)
    mean_flag = true;
    popntable.mean = nan(length(rowvalues),length(colvalues)); % initialise output
end;

ind = find(strcmp('median',varargin));           % Compute mean
if ~isempty(ind)
    median_flag = true;
    popntable.median = nan(length(rowvalues),length(colvalues)); % initialise output
end;

ind = find(strcmp('sd',varargin));             % s.d.
if ~isempty(ind)
    sd_flag = true;
    popntable.sd = nan(length(rowvalues),length(colvalues));% initialise output
end;

ind = find(strcmp('minn',varargin));           % minimum number of units to include condition.
if ~isempty(ind)
    min_n = varargin{ind+1};
end;

ind = find(strcmp('fn',varargin));
if ~isempty(ind)
    for i=1:length(ind)
        n_field_fns = n_field_fns +1;
        field_fns{n_field_fns} = varargin{ind(i)+1};
    end;
end;

ind = find(strcmp('hist',varargin));           % Make a histogram.
if ~isempty(ind)
    hist_flag = true;
    hist_bins = varargin{ind+1};
end;

ind = find(strcmp('roworder',varargin));           % Make a histogram.
if ~isempty(ind)
    roworder = varargin{ind+1};
else 
    roworder = 1:length(rowvalues);
end;
popntable.rowValues = rowvalues(roworder);


% Go through each row.
for ri = 1:length(rowvalues)

    % Go through each column
    for ci = 1:length(colvalues)
        
        % Find the values that satisfy the row values. 
        if iscell(rowvalues) && ischar(rowvalues{1})
            r = strcmp(rowvalues{ri},stackdata.(rowname));  % cells of strings
        else
            r = (stackdata.(rowname) == rowvalues(ri));     % numbers
        end;
            
        % Find the values that satisfy the column values. 
        if iscell(colvalues) && ischar(colvalues{1})
            c = strcmp(colvalues{ci},stackdata.(colname)); % cells of strings
        else
            c = (stackdata.(colname) == colvalues(ci));    % numbers
        end;
        
%         % Check any other inclusion criteria.
%         if n_field_fns>0
%             for fi=1:length(field_fns)
%                 f(fi,:) = eval(['stackdata.' field_fns{fi} ]);
%             end;
%         else 
%             f = ones(size(r));
%         end;

        % Check any other inclusion criteria.
        if n_field_fns>0
            for fi=1:length(field_fns)
                fname = fieldnames(stackdata);
                fieldchk = find(cellfun(@(x) ~isempty(findstr(field_fns{fi},x)), fname));
                if  length(fieldchk)>0
                    % Always select out the longest field name which matches.
                    [maxfield maxi] = max(cellfun(@(x) length(x),fname(fieldchk)));
                    fieldchk = fieldchk(maxi);

                    %fname =  fieldList_common{find(fieldchk)};
                    ind = findstr(field_fns{fi},fname{fieldchk});
                    if ind==1
                        cmdstr = ['stackdata.' field_fns{fi}];
                    else
                        cmdstr = [field_fns{fi}(1:ind-1) 'stackdata.' field_fns{fi}(ind:end)];
                    end;
                    f(fi,:) = eval(cmdstr);
                end;
            end;
        else 
            f = ones(size(r));
        end;

        
        
        % Compute where both are satisfied.
        inds = ( r & c & prod(f,1)==1);
        n = sum(inds);              % How many are found.
                
                                    % Check if min_n is met.
        if n>=min_n
             popntable.n(find(roworder==ri),ci) = n;        % Store n.
            
            % Calculate mean.
            if mean_flag 
                popntable.mean(find(roworder==ri),ci) = nanmean( stackdata.(valuename)(inds) );
            end;

            % Calculate mean.
            if median_flag 
                popntable.median(find(roworder==ri),ci) = nanmedian( stackdata.(valuename)(inds) );
            end;

            % Calculate sd.
            if sd_flag 
                popntable.sd(find(roworder==ri),ci) = nanstd( stackdata.(valuename)(inds) );
            end;

            % Calculate mean.
            if hist_flag 
                popntable.hist(find(roworder==ri),ci,[1:length(hist_bins)]) = histc( stackdata.(valuename)(inds), hist_bins );
            end;
        end;
    
    end;
end;    




