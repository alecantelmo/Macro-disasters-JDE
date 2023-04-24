function [ M ] = gaussian_moments( n_e )
%Calculate cross moments of a stochastic vector of size n_e-by-1
%distributed normally with mean zero(n_e,1) and variance matrix eye(n_e).
%The cross moments are defined:
%M2=E(kron(x,x))
%M3=E(kron(x,x,x))
%M4=E(kron(x,x,x,x))
%M5=E(kron(x,x,x,x,x))
%They are stored as fields of struct M.
%
% © Copyright, Oren Levintal, June 13, 2016.

M.M2=speye(n_e); M.M2=M.M2(:);
M.M3=sparse(n_e^3,1);
M.M4=sparse(n_e^4,1);
for i=1:n_e
    M.M4(sub2ind([n_e,n_e,n_e,n_e],i,i,i,i))=3;
end
for i=1:n_e-1
    for j=i+1:n_e
        tempP=perms([i,i,j,j]);
        M.M4(sub2ind([n_e,n_e,n_e,n_e],tempP(:,1),tempP(:,2),tempP(:,3),tempP(:,4)))=1;
    end
end

M.M5=sparse(n_e^5,1);

end

