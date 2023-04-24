function [ A ] = ttranspose( A )
%transpose a 2D sparse tensor
%
% © Copyright, Oren Levintal, June 13, 2016.

if length(A.tsize)>2
    error('tensors or rank higher than 2 are not supported')
end

m=A.tsize(1);
n=A.tsize(2);

A.ptr2d=intarray([1,m]);
A=ptr2d(A,1,m); 
A=ptr2col(A,1);
A=col2ptr(A,2);
A=ptr1d(A);

A=rmfield(A,'ptr2d');

end

