function report = compareUnitInfo(info1,info2,dontcarefields)
% compareUnitInfo comapres two unitinfo structures. 

% Get the list of fields in 1
i1fields = fieldnames(info1);

% Get the list of fields in 2.
i2fields = fieldnames(info2);

% Find the intersect and unique fields.
commonfields = intersect(i1fields,i2fields);
exclfields = setxor(i1fields,i2fields);

% Find the intersect and unique fields.
commonfields = intersect(i1fields,i2fields);
exclfields = setxor(i1fields,i2fields);

% Exclude dontcarefields from further analysis.
if nargin>2
    fieldswecareabout = ~ismember(commonfields,dontcarefields);
    commonfields = commonfields(fieldswecareabout);
end;

% Report any fields which are not the same.
if ~isempty(exclfields) && sum(strcmp(dontcarefields,'ignoreuniquefields'))==0
    fprintf('Found fields unique to either structure:\n');
    fprintf(' %s\n',exclfields{:});
end;

classmatch = nan(size(commonfields));
sizematch =  nan(size(commonfields));
valuematch = nan(size(commonfields));

% Loop through the intersec
for i=1:length(commonfields)

    if ~isempty(info1.(commonfields{i})) && ~isempty(info2.(commonfields{i}))
            % Can't compare if one isempty. But if neither are empty...
        %fprintf('%s',commonfields{i});
        
        % check the type
        if strcmp(class(info1.(commonfields{i})),class( info2.(commonfields{i})))
            classmatch(i) = 1;
        else
            classmatch(i) = 0;
            fprintf('%s... are different classes.\n',commonfields{i});
        end;

        % Check the size.
        if size(info1.(commonfields{i})) == size( info2.(commonfields{i}))
            sizematch(i) = 1;
            fieldsize = size(info1.(commonfields{i}));
        else
            sizematch(i) = 0;    
            fieldsize = nan;
            fprintf('%s... are different sizes.\n',commonfields{i});
        end;

    else 
        if isempty(info1.(commonfields{i})) && isempty(info2.(commonfields{i}))
            sizematch(i) = 1;
        else 
            sizematch(i) = 0;
        end;
         classmatch(i) = 0;
    end;
    
    if sizematch(i) == 0 || classmatch(i) ==0
        % If not the same size can't compare.
        continue;
    end;
        
    % Check the values are the same
    err = nan(fieldsize);
    for j =1:prod(fieldsize)
        if iscell(  info1.(commonfields{i}) )
            err(j)  = comVal( info1.(commonfields{i}){j}, info2.(commonfields{i}){j} );
        elseif isnumeric( info1.(commonfields{i}) ) || ischar(info1.(commonfields{i}) )
            err(j)  = comVal( info1.(commonfields{i})(j), info2.(commonfields{i})(j) );
        end;
    end;
    
    if sum(err)==0 
        valuematch(i) = 1;
        %fprintf('... are the same.\n',commonfields{i});
    else 
        valuematch(i) = 0;
        fprintf('%s have different values.\n',commonfields{i});
        fprintf('%d,',info1.(commonfields{i})); fprintf('\n');
        fprintf('%d,',info2.(commonfields{i})); fprintf('\n');
    end;
end;

notnan = ~isnan(valuematch);

report.fieldwhichdontmatch = commonfields(classmatch(notnan) & sizematch(notnan) & ~valuematch(notnan));
report.fieldsofdifferentsize = commonfields(classmatch & ~sizematch);
report.commonfields = commonfields;
report.classmatch = classmatch;
report.sizematch = sizematch;
report.valuematch = valuematch;
report.exclusivefields = exclfields;




% Compare one element. 
function err = comVal(el1,el2)

    % Check the values are the same
    if isnumeric( el1 ) 
        err = abs(el1-el2);
        if isnan(el1) && isnan(el2)
            err = 0;
        end;
    elseif isstr( el1 )
        if strcmp(el1, el2)
            err = 0;
        else
            err = 1;
        end;
    end;



        
        
        
