function [Lamda Phi]=largematrices(N,M,n_in,n_out,nb,A,B,C,D)
Phi = zeros((N)*n_out,(M)*n_in);
AB = zeros((N)*n_out,n_in);
% AB(1:n_out,:) = D;
% Lamda(1:n_out,:) = C;
AA=eye(nb);
for j = 1:N,
AB((j-1)*n_out+1:(j)*n_out,:) = C*AA*B; %C*A^(j-2)*B;
AA=AA*A;
Lamda(1+(j-1)*n_out:j*n_out,:) = C*AA; %C*A^(j-1);
end;
for j = 1:M,
Phi(1+(j-1)*n_out:end,1+(j-1)*n_in:(j)*n_in) = AB(1:(N-j+1)*n_out,:);
end;
Phi = sparse(Phi);

% Phi = zeros((N)*nb,(M)*n_in);
% AB = zeros((N)*nb,n_in);
% % AB(1:n_out,:) = D;
% % Lamda(1:n_out,:) = C;
% AA=eye(nb);
% for j = 1:N,
% AB((j-1)*nb+1:(j)*nb,:) = AA*B; %C*A^(j-2)*B;
% AA=AA*A;
% Lamda(1+(j-1)*nb:j*nb,:) = AA; %C*A^(j-1);
% end;
% for j = 1:M,
% Phi(1+(j-1)*nb:end,1+(j-1)*n_in:(j)*n_in) = AB(1:(N-j+1)*nb,:);
% end;
% Phi = sparse(Phi);