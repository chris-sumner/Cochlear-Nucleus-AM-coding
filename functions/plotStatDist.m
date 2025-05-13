function plotStatDist(stat,labelstr)

[v,b] = hist(stat,100);
bar(b,v); hold on;
csnorm = cumsum(v)/sum(v);
plot(b,max(v)*csnorm);
line([1 1]*b(max(find(csnorm<0.5))),ylim,'linestyle','--','color','r');
xlabel(labelstr);
ylabel('Count');
title(labelstr);


end

