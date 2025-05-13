function score = logLinkFn(metric,p0,x0)
% logLinkFn classifier score from metric.

score = p0 +(1-p0)./(1 + exp(-(metric-x0)));

% Could add a global look up table for exp to speed it up. 