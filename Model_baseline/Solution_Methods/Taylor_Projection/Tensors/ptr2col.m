function [ R ] = ptr2col( ten,coli )
%A=ptr2col(A,j) converts a 2D pointer to 1D pointer.
%
% © Copyright, Oren Levintal, June 13, 2016.

if coli>length(ten.tsize)
    error('column dimension too large')
end

l1=ten.ptr2d(1);
l2=ten.ptr2d(2);
[ newcol ] = ptr2col_mex( ten.ptr,ten.ptr(end)-1 );
R=sptensor;
R.vals=ten.vals;
R.ptr=[ten.ptr(:,1);ten.ptr(end)];
R.cols=[ten.cols(:,1:coli-1),newcol,ten.cols(:,coli:end)];
R.tsize=[l1,ten.tsize(3:coli+1),l2,ten.tsize(coli+2:end)];


end

