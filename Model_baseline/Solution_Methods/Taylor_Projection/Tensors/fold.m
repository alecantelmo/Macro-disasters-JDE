function [R]=fold(A,varargin)
% reshape a sparse matrix of dimensions m,n1*...*n4 to a tensor of
% dimensions m,n1,...,n4
%
% © Copyright, Oren Levintal, June 13, 2016.

if length(A.tsize)>2
    A=unfold(A);
end
R=A;
if nargin==3
    n1=varargin{1};
    n2=varargin{2};
    [cols1,cols2]=ind2sub([n1,n2],A.cols);
    R.cols=[cols1,cols2];
    R.tsize=[A.tsize(1),n1,n2];
elseif nargin==4
    n1=varargin{1};
    n2=varargin{2};
    n3=varargin{3};
    [cols1,cols2,cols3]=ind2sub([n1,n2,n3],A.cols);
    R.cols=[cols1,cols2,cols3];
    R.tsize=[A.tsize(1),n1,n2,n3];
elseif nargin==5
    n1=varargin{1};
    n2=varargin{2};
    n3=varargin{3};
    n4=varargin{4};
    [cols1,cols2,cols3,cols4]=ind2sub([n1,n2,n3,n4],A.cols);
    R.cols=[cols1,cols2,cols3,cols4];
    R.tsize=[A.tsize(1),n1,n2,n3,n4];
elseif nargin>5
    error('tensor of rank higher than 5 not supported')
end

R.cols=intarray(R.cols);
R.tsize=intarray(R.tsize);

end