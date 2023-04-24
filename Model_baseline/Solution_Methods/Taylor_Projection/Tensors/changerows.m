function [A]=changerows(A,newrows,rowdim)
%Change the rows of tensor A
%
% © Copyright, Oren Levintal, June 13, 2016.

if isfield(A,'ptr2d')
    error('2D pointer not supported')
end
rowcount=zeros(rowdim,1);
rowcount(newrows)=diff(A.ptr);
A.ptr=intarray(cumsum([1;rowcount]));
A.tsize(1)=intarray(rowdim);


end