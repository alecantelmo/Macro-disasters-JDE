function R=AkronkC(A,C,k)
%Perform R=A*kron(C,...,C), where C is multiplied k times by the Kronecker operator.
%
% © Copyright, Oren Levintal, June 13, 2016.

n=size(A,1);
m=size(C,1);
for j=1:k
    A=reshape(A,numel(A)/m,m);
    A=(A*C)';
end

R=reshape(A,numel(A)/n,n)';


end