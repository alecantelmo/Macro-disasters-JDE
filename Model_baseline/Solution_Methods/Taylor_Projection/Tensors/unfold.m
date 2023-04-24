function [R]=unfold(A)
% unfold a sparse tensor of dimensions m,n1,...,n4 to a matrix of
% dimensions m,n1*...*n4
% the pointer can be 1D or 2D
%
% © Copyright, Oren Levintal, June 13, 2016.

if size(A.ptr,2)==1 %1D pointer
    rowdim=A.tsize(1);
    coldim=A.tsize(2:end);
    rank=length(A.tsize);

else % 2D pointer
    rowdim=A.tsize(1:2);
    coldim=A.tsize(3:end);
    rank=length(A.tsize)-1; 
end
% check overflow on 32-bit system
n=prod(int64(coldim));
if intarray(n)==intarray(n-1)
    error('integer overflow')
end

R=A;
R.tsize=intarray([rowdim,n]);

if rank==3
    R.cols=sub2ind(coldim,A.cols(:,1),A.cols(:,2));
elseif rank==4
    R.cols=sub2ind(coldim,A.cols(:,1),A.cols(:,2),A.cols(:,3));
elseif rank==5
    R.cols=sub2ind(coldim,A.cols(:,1),A.cols(:,2),A.cols(:,3),A.cols(:,4));
elseif rank>5
    error('tensor of rank higher than 5 not supported')
end
R.cols=intarray(R.cols);

end