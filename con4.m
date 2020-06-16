function [c,ceq]=con4(x)
global mu H

h2=hessiand3(x);
 hU=hessiandU3(x);

% c=det(((h2-eye(6)/0.1)));
% h2=hessiand(x);
h1=pinv(H+mu*hU)*(-h2*mu+eye(size(H,1)));%+h2'-2*eye(6)/mu;
h=h1+h1';
c=(-h);
ceq=[];

