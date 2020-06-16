function J=Barrier2(U)
global y W N H  G d_relax mu
n_in = 1;%size(U,1)
LU=eye(N*n_in);

G=[LU;-LU];
S=zeros(1);
x=y;

b=0;

for i=1:size(W,1)
% if (W(i)-(G(i,:)*U))<d_relax(i)
%     b=b+0.5*(((-G(i,:)*U+W(1)-2*d_relax(i))/d_relax(i))^2-1)-log(d_relax(i))+log(W(i))-1/W(i)*(G(i,:)*U);
% else
    b=b+(-log(-G(i,:)*U+W(i))+log(W(i))-1/W(i)*(G(i,:)*U));
% end
end

% 
% if (W(4)-(G(4,:)*U+S(4,:)*x))<d
%     b=b+0.5*(((-G(4,:)*U-S(4,:)*x+W(4)-2*d)/d)^2-1)-log(d)+log(W(4))-1/W(1)*(G(4,:)*U+S(4,:)*x);
% else
%     b=b+(-log(-G(4,:)*U-S(4,:)*x+W(4))+log(W(4))-1/W(4)*(G(4,:)*U+S(4,:)*x));
% end
% 
%     b=b+(-log(-G(2,:)*U-S(2,:)*x+W(2))+log(W(2))-G(2,:)/W(2)*(U));
% 
% % 
%     b=b+(-log(-G(3,:)*U-S(3,:)*x+W(3))+log(W(3))-G(3,:)/W(3)*(U));
%     
% %     b=b+(-log(-G(4,:)*U-S(4,:)*x+W(4))+log(W(4))-G(4,:)/W(4)*(U));
% 
%     b=b+(-log(-G(5,:)*U-S(5,:)*x+W(5))+log(W(5))-G(5,:)/W(5)*(U));
%    
%  b=b-log(-[ones(1,6)]*U+.2*W(6)) +log(.2*W(6))-[ones(1,6)]*U/0.2;
%  b=b-log([ones(1,6)]*U+.2*W(6))+log(.2*W(6))+[ones(1,6)]*U/0.2;
% 
 J=0.5*U'*H*U-x'*U+mu*b;



    