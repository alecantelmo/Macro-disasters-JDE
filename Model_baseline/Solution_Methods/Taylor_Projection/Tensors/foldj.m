function [R]=foldj(A,j,dim)
% [R]=foldj(A,j,dim) reshapes the j'th col dimension of a sparse tensor A
% do new dimension dim.
%
% © Copyright, Oren Levintal, June 13, 2016.

R=A;
ndim=length(dim);
dim=intarray(dim);
if ndim==2
    n1=dim(1);
    n2=dim(2);
    [cols1,cols2]=ind2sub([n1,n2],A.cols(:,j));
    R.cols=[A.cols(:,1:j-1),cols1,cols2,A.cols(:,j+1:end)];
    R.tsize=[A.tsize(1:j),n1,n2,A.tsize(j+2:end)];
elseif ndim==3
    n1=dim(1);
    n2=dim(2);
    n3=dim(3);
    [cols1,cols2,cols3]=ind2sub([n1,n2,n3],A.cols(:,j));
    R.cols=[A.cols(:,1:j-1),cols1,cols2,cols3,A.cols(:,j+1:end)];
    R.tsize=[A.tsize(1:j),n1,n2,n3,A.tsize(j+2:end)];
elseif ndim==4
    n1=dim(1);
    n2=dim(2);
    n3=dim(3);
    n4=dim(4);    
    [cols1,cols2,cols3,cols4]=ind2sub([n1,n2,n3,n4],A.cols(:,j));
    R.cols=[A.cols(:,1:j-1),cols1,cols2,cols3,cols4,A.cols(:,j+1:end)];
    R.tsize=[A.tsize(1:j),n1,n2,n3,n4,A.tsize(j+2:end)];
elseif ndim>4
    error('tensor of rank higher than 5 not supported')
end

R.cols=intarray(R.cols);
R.tsize=intarray(R.tsize);

end