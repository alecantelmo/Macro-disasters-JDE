function [C,ind]=AkronB1U3B2(A,B1,B2,ind,n_ind,maxload,sum)
%The function calculates C=A*kron(kron(B1,B1,B1)*U3,B2), where A is a
%symmetric tensor with A.tsize=[l,m,m,m,m], B1 is a tensor of size m,n1 and
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
    ind=cell(6,1);
end

l=A.tsize(1);
m=A.tsize(2);
n1=B1.tsize(2);
n2=B2.tsize(2);

unique3=nchoosek(n1+2,3);

A=ptr1d(col2ptr(A,1));

% do A*kron(B1,B1,B1)*U3
[C,ind{2}]=AkronBU3(A,B1,ind{2},n_ind,maxload,sum1);


C=ptr2d(C,l,m);
C=ptr2col(C,2);
C=col2ptr(C,1);
C=ptr1d(C);

% second product

[C,ind{5}]=contraction1(C,B2,ind{5},n_ind,maxload,sum2);


C=ptr2d(C,l,unique3);
C=unfold(ptr2col(C,2));

end