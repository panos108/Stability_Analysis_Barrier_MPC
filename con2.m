function [c,ceq]=con2(x)
global mu
h2=hessiand2(x);
% c=det(((h2-eye(6)/0.1)));
    % h2=hessiand(x);
h=h2'+h2;%-2*eye(2)/0.0001;
c=(h);
ceq=[];
% ceq=[];