function p = softMaxEval(pars,B)
% softMaxEval computes a confusion matrix from softmax function pars.
%
% Two forms:
%   p = softMaxEval(w,B)
%       w is a vector of weights - a channel for detecting each
%       class.
%       B is a vector of bias values.
%
%   p = softMaxEval(pars)
%       where pars = [ w B ]

if nargin<2
    l = length(pars);
    l2 = l/2;
    x = pars(1:l2);
    B = pars(l2+1:l);
elseif strcmp(B,"nobias")
    l = length(pars);   
    l2 = l;
    x = pars;
    B = zeros(1,l);
else
    x = pars;
    l2 = length(pars);
end;

sizeB = size(B);
if sizeB(1) ~= 1
    error('B must be a row vector');
end;

% Compute the confusion matrix from the softmax function.  
Z = x.*eye(length(x)) + B';
p = exp(Z)./(ones(l2,1)*sum(exp(Z)));
