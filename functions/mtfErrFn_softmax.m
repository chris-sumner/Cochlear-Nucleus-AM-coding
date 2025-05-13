function err = mtfErrFn_softmax(f_fn,tmp)
% mtfErrFn_softmax

% Set up functions and tables.
gainfn = @(g) g*f_fn;   % Matrix: set of all gain values for each frequencies.
ngain = length(tmp.gains);
pcfn = @(lt,zs,cis) softMaxLinkFn(lt,zs,cis);        % Matrix: put any post gain non-linearity here. Does notthing for linear. 
gaintable = gainfn(tmp.gains');  % Generic table of functions for all neurons 

lU = length(tmp.unitlist);

% Preallocate outputs in tmp.
tmp.best_ind = nan(lU,1);
tmp.fitted_err = nan(lU,1);

% Go through all neurons and compute the error. 
for ui = 1:lU   
    % A table of possible functions - for the frequencies used in this dataset.
    mclass = length(tmp.amfreq_inds{ui});
    tmp.measuretable = pcfn( ...
        gaintable(:,tmp.amfreq_inds{ui}), ...
         [1:mclass], tmp.DataIndex{ui}(1));      

    % Mean sqr difference between the actual MTF values and each row in the
    % table.
    tmp.difference = nanmean((ones(ngain,1)*tmp.measure{ui} - tmp.measuretable).^2,2);
    % Find the best gain for this dataset.
    best_ind = find(tmp.difference == min(tmp.difference) );
    tmp.best_ind(ui) = best_ind(1);
    % Store the error for this dataset.
    tmp.fitted_err(ui) = tmp.difference(tmp.best_ind(ui));
end;
 % Sum of all the mean squared errors for each dataset.
err = sum( tmp.fitted_err );

