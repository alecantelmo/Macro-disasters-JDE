function [R]=tvec(A)
% vectorize a sparse tensor of dimensions m,n1,...,n4 to dimensions
% m*n1*...*n4,1.
%
% © Copyright, Oren Levintal, June 13, 2016.

A=unfold(A);

A=ptr2d(A,1,A.tsize(1));
A=ptr2col(A,1);
A=fold(A,prod(A.tsize),1);
A=col2ptr(A,1);
A=ptr1d(A);
R=rmfield(A,'ptr2d');

end