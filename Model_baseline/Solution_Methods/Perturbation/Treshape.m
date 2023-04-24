function [ rmat ] = Treshape( mat,m,n )
%Treshape(A,m,n) transposes matrix A and then reshapes it into m,n
%dimensions.
%
% © Copyright, Oren Levintal, June 13, 2016.

[j0,i0,s]=find(mat);
[n0,m0]=size(mat);

ind=(j0-1)*m0+i0;
[i1,j1]=ind2sub([m,n],ind);
rmat=sparse(i1,j1,s,m,n);

end

