function h2=hessiand2(x)
global b G
L=G;
h2=zeros(size(L,2),size(L,2));
for i=1:size(L,1)
%     if b(i)-G(i,:)*x'>=0.0001
        h2=h2+G(i,:)'*G(i,:)*1/(b(i)-G(i,:)*x')^2;
%     else
%         h2=h2+G(i,:)'*G(i,:)*1/(0.0001)^2;
%     end
end
