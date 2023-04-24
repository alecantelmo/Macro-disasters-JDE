function [ A ] = tplus2d( B1,B2,varargin )
%adds two tensors with two dimensional pointers
%
% © Copyright, Oren Levintal, June 13, 2016.

if size(B1.ptr,2)==1
    error('pointer should be 2-D')
end
rank1=length(B1.tsize);
rank2=length(B2.tsize);
if ~isequal(rank1,rank2)
    error('incompatible dimensions')
end


A=sptensor;
A.ptr2d=B1.ptr2d;
A.tsize=B1.tsize;

if nargin==2
    maxload=intarray(1);
elseif nargin==3
    maxload=intarray(varargin{1});
end

if rank1==3
    %     A=tplus1(B1,B2);
    l1=B1.tsize(1);
    l2=B1.tsize(2);
    m1=B1.tsize(3);
    IB1=B1.ptr;
    IB2=B2.ptr;
    JB1=B1.cols;
    JB2=B2.cols;
    NB1=B1.vals;
    NB2=B2.vals;
    [A.vals,A.ptr,A.cols] = tplus1_2D_mex(l1,l2,m1,IB1,JB1,NB1,IB2,JB2,NB2,maxload);
    A.vals=A.vals(:,1:A.ptr(end)-1);
    A.cols=A.cols(1:A.ptr(end)-1,:);
%     A.ptr2d=1;
% elseif rank1==3
%     l=B1.tsize(1);
%     m=B1.tsize(2:end);
%     IB1=B1.ptr;
%     IB2=B2.ptr;
%     JB1=B1.cols;
%     JB2=B2.cols;
%     NB1=B1.vals;
%     NB2=B2.vals;
%     [A.vals,A.ptr,A.cols] = tplus2_mex(l,m,IB1,JB1,NB1,IB2,JB2,NB2,maxload);
%     A.vals=A.vals(:,1:A.ptr(end)-1);
%     A.cols=A.cols(1:A.ptr(end)-1,:);
% elseif rank1==4
%     l=B1.tsize(1);
%     m=B1.tsize(2:end);
%     IB1=B1.ptr;
%     IB2=B2.ptr;
%     JB1=B1.cols;
%     JB2=B2.cols;
%     NB1=B1.vals;
%     NB2=B2.vals;
%     [A.vals,A.ptr,A.cols] = tplus3_mex(l,m,IB1,JB1,NB1,IB2,JB2,NB2,maxload);
%     A.vals=A.vals(:,1:A.ptr(end)-1);
%     A.cols=A.cols(1:A.ptr(end)-1,:);
% elseif rank1==5
%     l=B1.tsize(1);
%     m=B1.tsize(2:end);
%     IB1=B1.ptr;
%     IB2=B2.ptr;
%     JB1=B1.cols;
%     JB2=B2.cols;
%     NB1=B1.vals;
%     NB2=B2.vals;
%     [A.vals,A.ptr,A.cols] = tplus4_mex(l,m,IB1,JB1,NB1,IB2,JB2,NB2,maxload);
%     A.vals=A.vals(:,1:A.ptr(end)-1);
%     A.cols=A.cols(1:A.ptr(end)-1,:);
else
    error('tensors of rank larger than 3 not supported')
end
    

end

