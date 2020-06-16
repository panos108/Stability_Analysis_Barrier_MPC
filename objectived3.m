function obj=objectived3(x)
%  [c,~]=con3(x);
% obj=.5*log((c));%sum(x.^2);%-0.001*log(c);
[c,~]=con4(x);
obj=-max(-eig(c));
