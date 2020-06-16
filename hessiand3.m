function h2=hessiand3(x)
global L b G S d
L=G;

x=x';
h2=zeros(size(L,2),size(L,2));
for i=1:size(L,1)
     if (b(i)-[G(i,:),S(i,:)]*x)>=d(i)
        h2=h2+G(i,:)'*S(i,:)*1/(b(i)-[G(i,:),S(i,:)]*x)^2;
    else
        h2=h2+G(i,:)'*S(i,:)*1/(d(i))^2;
    end
end
