function y = logisticPDF(x,mu)
% logisticPDF compute a logistic PDF
%
% y = logisticPDF(x,mu)

x2 = x - mu;
y = exp(-x2)./((1 + exp(-x2)).^2);
