function err = mtfErrFn_linear(f_fn,tmp)
% mtfErrFn_linear

% Set up functions and tables.
gainfn = @(g) g*f_fn;   % Matrix: set of all gain values for each frequencies.
ngain = length(tmp.gains);
pcfn = @(lt) lt;        % Matrix: put any post gain non-linearity here. Does notthing for linear. 
gaintable = gainfn(tmp.gains');  % Generic table of functions for all neurons 

% Go through all neurons and comute the error. 
for ui = 1:length(tmp.unitlist)   
    % A table of possible functions - for the frequencies used in this dataset.
    tmp.measuretable = pcfn(gaintable(:,tmp.amfreq_inds{ui}));        
    % Mean sqr difference between the actual MTF values and each row in the
    % table.
    tmp.difference = nanmean((ones(ngain,1)*tmp.measure{ui} - tmp.measuretable).^2,2);
    % Find the best gain for this dataset.
    tmp.best_ind(ui) =find(tmp.difference == min(tmp.difference) );
    % Store the error for this dataset.
    tmp.fitted_err(ui) = tmp.difference(tmp.best_ind(ui));
end;
 % Sum of all the mean squared errors for each dataset.
err = sum( tmp.fitted_err );

