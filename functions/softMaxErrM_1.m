function err = softMaxErrM_1(pars,pdata,smoothw,Z_max,B0ind)


% This version sets one of the B values to zero. 
% which ever is specified by B0ind.

l = length(pars);
l2 = ceil(l/2);
x = [pars(1:l2)];
l3 = length(x);
% If not specified - B0 is the final value.
if nargin<5 
    B0ind = l2;
end;
% Insert a zero for B0
if B0ind == 1
    B = [0 pars(l2+1:l)];
else
    B = [pars(l2+1:l2+B0ind-1) 0 pars(l2+B0ind:l)];
end;

Z = x.*eye(l3) + B';
pmodel = exp(Z)./(ones(l3,1)*sum(exp(Z)));

errmat = (pdata - pmodel).^2;
if nargin<3
    err = sum(errmat(:));
else
    nB = sum(B<0);
    nP = sum(abs(pars)>Z_max);    
    penaliser = smoothw*( sum(diff(pars).^2) ) + nP  + nB;
    err = sum(errmat(:)) + penaliser ;
end;




