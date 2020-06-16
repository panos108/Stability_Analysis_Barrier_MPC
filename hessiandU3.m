function h2=hessiandU3(x)
global L G S d b LU
L=G;
x=x';
h2=zeros(size(L,2),size(L,2));
L1=[LU];ones(1,6);%eye(6);
for i=1:size(L1,1);
    
    if 1-L1(i,:)*x(1:size(G,2))>=.00001
        
        h2=h2+L1(i,:)'*L1(i,:)*1/(1-L1(i,:)*x(1:size(G,2)))^2;
    else
        h2=h2+L1(i,:)'*L1(i,:)*1/(0.00001)^2;
    end
    if 1+L1(i,:)*x(1:size(G,2))>=0.00001
        h2=h2+L1(i,:)'*L1(i,:)*1/(1+L1(i,:)*x(1:size(G,2)))^2;
    else
        h2=h2+L1(i,:)'*L1(i,:)*1/(0.00001)^2;
    end
    
end

for i=1:size(L,1)
    if b(i)-[G(i,:),S(i,:)]*x>=d(i)
        h2=h2+G(i,:)'*G(i,:)*1/(b(i)-[G(i,:),S(i,:)]*[x])^2;
    else
        h2=h2+G(i,:)'*G(i,:)*1/(d(i))^2;
    end
end
