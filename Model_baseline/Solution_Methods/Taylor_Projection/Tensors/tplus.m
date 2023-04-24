function [ A ] = tplus( B1,B2,varargin )
%adds two tensors
%
% © Copyright, Oren Levintal, June 13, 2016.

sB1=size(B1.vals,1);
sB2=size(B2.vals,1);

if sB1==1
    if sB2>1
        B1.vals=repmat(B1.vals,sB2,1);
    end
elseif sB1>1
    if sB2==1
        B2.vals=repmat(B2.vals,sB1,1);
    elseif sB2>1 && sB1~=sB2
        error('incompatible states')
    end
end
    

if size(B1.ptr,2)>1 % 2D pointers
    if nargin==2
        maxload=intarray(1);
    elseif nargin==3
        maxload=intarray(varargin{1});
    end
    [ A ] = tplus2d( B1,B2,maxload );
else

    rank1=length(B1.tsize);
    rank2=length(B2.tsize);
    if ~isequal(rank1,rank2)
        error('incompatible dimensions')
    end


    A=sptensor;
    A.tsize=B1.tsize;

    l=B1.tsize(1);
    m=B1.tsize(2:end);
    IB1=B1.ptr;
    IB2=B2.ptr;
    JB1=B1.cols;
    JB2=B2.cols;
    NB1=B1.vals;
    NB2=B2.vals;
    if ~isreal(NB1) || ~isreal(NB2)
        error('complex numbers not supported')
    end

    if nargin==2
        maxload=intarray(1);
    elseif nargin==3
        maxload=intarray(varargin{1});
    end

    if rank1==2
        %     A=tplus1(B1,B2);
        [A.vals,A.ptr,A.cols] = tplus1_mex(l,m,IB1,JB1,NB1,IB2,JB2,NB2,maxload);
        A.vals=A.vals(:,1:A.ptr(end)-1);
        A.cols=A.cols(1:A.ptr(end)-1,:);
    elseif rank1==3
        [A.vals,A.ptr,A.cols] = tplus2_mex(l,m,IB1,JB1,NB1,IB2,JB2,NB2,maxload);
        A.vals=A.vals(:,1:A.ptr(end)-1);
        A.cols=A.cols(1:A.ptr(end)-1,:);
    elseif rank1==4
        [A.vals,A.ptr,A.cols] = tplus3_mex(l,m,IB1,JB1,NB1,IB2,JB2,NB2,maxload);
        A.vals=A.vals(:,1:A.ptr(end)-1);
        A.cols=A.cols(1:A.ptr(end)-1,:);
    elseif rank1==5
        [A.vals,A.ptr,A.cols] = tplus4_mex(l,m,IB1,JB1,NB1,IB2,JB2,NB2,maxload);
        A.vals=A.vals(:,1:A.ptr(end)-1);
        A.cols=A.cols(1:A.ptr(end)-1,:);
    else
        error('tensors of rank larger than 5 not supported')
    end
    
end

end

