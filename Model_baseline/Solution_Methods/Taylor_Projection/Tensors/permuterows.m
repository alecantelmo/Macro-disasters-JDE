function [A]=permuterows(A,newrows)
%Permute the rows of tensor A by newrows.
%
% © Copyright, Oren Levintal, June 13, 2016.

if isfield(A,'ptr2d')
    error('2D pointer not supported')
end

newrows=intarray(newrows);

[A.ptr,A.cols,A.vals] = perm_rows(A.ptr,A.cols,A.vals,A.tsize(1),newrows);

end