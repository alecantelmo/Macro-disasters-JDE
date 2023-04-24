function [C,ind]=AkronB1U2B2(A,B1,B2,ind,n_ind,maxload,sum)
%The function calculates C=A*kron(kron(B1,B1)*U2,B2), where A is a
%symmetric tensor with A.tsize=[l,m,m,m], B1 is a tensor of size m,n1 and
%B2 is a tensor of size m,n2.
%The result is returned as a 2D sparse tensor
%
% © Copyright, Oren Levintal, June 13, 2016.

if strcmp(sum,'sum')
    sum1='vec';
    sum2='sum';
elseif strcmp(sum,'vec')
    sum1='vec';
    sum2='vec';
end

if isempty(ind)
    ind=cell(5,1);
end
l=A.tsize(1);
m=A.tsize(2);
n1=B1.tsize(2);


A=ptr1d(col2ptr(A,1));

% do A*kron(B1,B1)*U2
[C,ind{2}]=AkronBU2(A,B1,ind{2},n_ind,maxload,sum1);


C=ptr2d(C,l,m);
C=ptr2col(C,2);
C=col2ptr(C,1);
C=ptr1d(C);

% second product

[C,ind{4}]=contraction1(C,B2,ind{4},n_ind,maxload,sum2);


C=ptr2d(C,l,nchoosek(n1+1,2));
C=unfold(ptr2col(C,2));

end