function err = softMaxErr(pars,pdata,smoothw,Z_max)

l = length(pars);
l2 = l/2;
x = pars(1:l2);
B = pars(l2+1:l);
Z = x.*eye(length(x)) + B';
pmodel = exp(Z)./(ones(length(x),1)*sum(exp(Z)));

errmat = (pdata - pmodel).^2;
if nargin<3
    err = sum(errmat(:));
else
    penaliser = smoothw*( sum(diff(pars).^2) ) + sum(abs(pars)>Z_max);
    err = sum(errmat(:)) + penaliser ;
end;

