function [ A,ind ] = extract(B,extractrows,extractcols,compress,ind)
%The function extract rows and columns from a tensor with symmetric column
%and returns a new sparse tensor. If compress=1 the output is in compressed
%form.
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(ind)
    [A,tempind1]=takerows(B,extractrows);
    [A,tempind2]=takecols(A,extractcols);
    if compress==1
        [A,tempind3]=takeunique(A);
        tempind2=tempind2(tempind3);
    end
    tempind=tempind1(tempind2);
    ind.ptr=A.ptr;
    ind.cols=A.cols;
    ind.tsize=A.tsize;
    ind.ind=tempind;
else
    A.vals=B.vals(:,ind.ind);
    A.ptr=ind.ptr;
    A.cols=ind.cols;
    A.tsize=ind.tsize;
end  

end

