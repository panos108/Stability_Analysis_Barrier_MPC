function obj=objectived_min(x)
%  [c,~]=con3(x);
% obj=.5*log((c));%sum(x.^2);%-0.001*log(c);
[c,~]=con2(x);


obj=min(eig(c));