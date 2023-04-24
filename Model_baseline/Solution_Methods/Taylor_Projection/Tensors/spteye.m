function T=spteye(n,varargin)
%
% © Copyright, Oren Levintal, June 13, 2016.

if nargin==1
    n_s=1;
elseif nargin==2
    n_s=varargin{1};
else
    error('wrong number of arguments')
end
% creates eye(T) as sptensor
n=double(n);
T.vals=ones(n_s,n);
n=intarray(n);
T.cols=(1:n)';
T.ptr=(1:n+1)';
T.tsize=[n,n];

end

