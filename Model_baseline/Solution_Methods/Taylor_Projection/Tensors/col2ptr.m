function [ R ] = col2ptr( ten,colptr )
%A=col2ptr(A,j) converts the pointer of A into a 2D pointer.
%A sparse tensor of size l,m1,m2,..,mk has a pointer for the row dimension
%l. This function makes a pointer for l,mj, which is similar to permuting
%and reshaping the tensor to size l*mj,m1,m2,..,mj-1,mj+1,..,mk.
%
% © Copyright, Oren Levintal, June 13, 2016.

if size(ten.ptr,2)>1
    error('pointer must be 1D')
end
l=ten.tsize(1);
coldim=ten.tsize(2:end);
colptr=intarray(colptr);
ptrcoldim=coldim(colptr)+1;

[ newptr,newcols,newvals ] = col2ptr_mex( ten.ptr,ten.cols,ten.vals,l,coldim,colptr,ptrcoldim );

R=sptensor;
R.ptr=newptr;
R.cols=newcols;
R.vals=newvals;
R.tsize=[l,coldim(colptr),coldim(1:colptr-1),coldim(colptr+1:end)];
R.ptr2d=R.tsize(1:2);

end

