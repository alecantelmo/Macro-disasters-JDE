function [ A,index ] = colsort( ten,varargin )
%sorts tensor columns
%
% © Copyright, Oren Levintal, June 13, 2016.

if size(ten.ptr,2)>1
    error('2D ptr not supported')
end

if nargin==2
    index=varargin{1};
elseif nargin==1
    index=[];
end


if isempty(index) % sorting index is not available
    if nargout==2 % regular sort with an index
        [i,j]=tfind(ten);
        [sortij,index]=sortrows([i,j]);
        A=sptensor(sortij(:,1),sortij(:,2:end),ten.vals(:,index),ten.tsize(1),ten.tsize(2:end)); 
    elseif nargout==1 % a mex sort without index
        for coli=length(ten.tsize)-1:-1:1 % sort each col dimension
            [ ten.cols,ten.vals ] = sortcol_mex( ten.ptr,ten.cols,ten.vals,ten.tsize(1),ten.tsize(2:end),intarray(coli),ten.tsize(coli+1) );
        end
        A=ten;
    end
else % sorting index is available
    A=ten;
    A.vals=ten.vals(:,index);
    A.cols=ten.cols(index,:);
end

A.colsorted=1;

end

