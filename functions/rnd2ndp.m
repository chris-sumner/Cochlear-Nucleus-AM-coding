function y = rnd2ndp(x,n)
% rnd2ndp    round a number or matrix to n decimal places

% Calculate conversion factor to make all decimal places to
% be kept whole numbers.
dpfactor = 10^n;

% Multiply by the factor
bob = x*dpfactor;

% round to whole numbers
bob = round(bob);

% divide by the factor to regain the correct numbers
bob = bob/dpfactor;


y=bob;
