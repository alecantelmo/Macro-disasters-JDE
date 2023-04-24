function [ OMEGA ] = create_OMEGA( n_x,varargin )
%The function creates sparse matrices OMEGA1,..,OMEGA9, that are used for
%the high-order multivariate chain rules described in Levintal, Oren,
%"Fifth Order Perturbation Solution to DSGE Models", 2014.
%
% © Copyright, Oren Levintal, June 13, 2016.

if length(varargin)==1
    order=varargin{1}; % do only specified order
else
    order=0; % do orders 3,4,5
end

if order==3 || order==0
    ind=[1:n_x^3];
    M=reshape(ind,1,n_x,n_x,n_x);
    Ix=speye(n_x^3);
    OMEGA.OMEGA1=Ix(:,ipermute(M,[1,3,4,2]))+Ix(:,ipermute(M,[1,2,4,3]))+Ix(:,ipermute(M,[1,2,3,4]));
end

if order==4 || order==0
    ind=[1:n_x^4];
    M=reshape(ind,1,n_x,n_x,n_x,n_x);
    Ix=speye(n_x^4);
    OMEGA.OMEGA2=Ix(:,ipermute(M,[1,4,5,3,2]))+Ix(:,ipermute(M,[1,3,5,4,2]))+Ix(:,ipermute(M,[1,2,5,4,3]))...
        +Ix(:,ipermute(M,[1,3,4,5,2]))+Ix(:,ipermute(M,[1,2,4,5,3]))+Ix(:,ipermute(M,[1,2,3,5,4]));
    OMEGA.OMEGA3=Ix(:,ipermute(M,[1,3,4,5,2]))+Ix(:,ipermute(M,[1,2,4,5,3]))+Ix(:,ipermute(M,[1,2,3,5,4]))+Ix(:,ipermute(M,[1,2,3,4,5]));
    OMEGA.OMEGA4=Ix(:,ipermute(M,[1,3,4,2,5]))+Ix(:,ipermute(M,[1,2,4,3,5]))+Ix(:,ipermute(M,[1,2,3,4,5]));
end

if order==5 || order==0
    ind=[1:n_x^5];
    M=reshape(ind,1,n_x,n_x,n_x,n_x,n_x);
    Ix=speye(n_x^5);
    OMEGA.OMEGA5=Ix(:,ipermute(M,[1,5,6,4,3,2]))+Ix(:,ipermute(M,[1,4,6,5,3,2]))+Ix(:,ipermute(M,[1,3,6,5,4,2]))+Ix(:,ipermute(M,[1,2,6,5,4,3]))...
        +Ix(:,ipermute(M,[1,4,5,6,3,2]))+Ix(:,ipermute(M,[1,3,5,6,4,2]))...
        +Ix(:,ipermute(M,[1,2,5,6,4,3]))+Ix(:,ipermute(M,[1,3,4,6,5,2]))+Ix(:,ipermute(M,[1,2,4,6,5,3]))+Ix(:,ipermute(M,[1,2,3,6,5,4]));
    OMEGA.OMEGA6=Ix(:,ipermute(M,[1,4,5,6,3,2]))+Ix(:,ipermute(M,[1,3,5,6,4,2]))+Ix(:,ipermute(M,[1,2,5,6,4,3]))+Ix(:,ipermute(M,[1,3,4,6,5,2]))...
        +Ix(:,ipermute(M,[1,2,4,6,5,3]))+Ix(:,ipermute(M,[1,2,3,6,5,4]))+Ix(:,ipermute(M,[1,3,4,5,6,2]))...
        +Ix(:,ipermute(M,[1,2,4,5,6,3]))+Ix(:,ipermute(M,[1,2,3,5,6,4]))+Ix(:,ipermute(M,[1,2,3,4,6,5]));
    OMEGA.OMEGA7=Ix(:,ipermute(M,[1,4,5,3,6,2]))+Ix(:,ipermute(M,[1,3,5,4,6,2]))+Ix(:,ipermute(M,[1,2,5,4,6,3]))...
        +Ix(:,ipermute(M,[1,3,4,5,6,2]))+Ix(:,ipermute(M,[1,2,4,5,6,3]))+Ix(:,ipermute(M,[1,2,3,5,6,4]))...
        +Ix(:,ipermute(M,[1,4,5,2,6,3]))+Ix(:,ipermute(M,[1,3,5,2,6,4]))+Ix(:,ipermute(M,[1,2,5,3,6,4]))...
        +Ix(:,ipermute(M,[1,3,4,2,6,5]))+Ix(:,ipermute(M,[1,2,4,3,6,5]))+Ix(:,ipermute(M,[1,2,3,4,6,5]))...
        +Ix(:,ipermute(M,[1,3,4,2,5,6]))+Ix(:,ipermute(M,[1,2,4,3,5,6]))+Ix(:,ipermute(M,[1,2,3,4,5,6]));
    OMEGA.OMEGA8=Ix(:,ipermute(M,[1,3,4,5,6,2]))+Ix(:,ipermute(M,[1,2,4,5,6,3]))+Ix(:,ipermute(M,[1,2,3,5,6,4]))+Ix(:,ipermute(M,[1,2,3,4,6,5]))+Ix(:,ipermute(M,[1,2,3,4,5,6]));
    OMEGA.OMEGA9=Ix(:,ipermute(M,[1,3,4,5,2,6]))+Ix(:,ipermute(M,[1,2,4,5,3,6]))+Ix(:,ipermute(M,[1,2,3,5,4,6]))+Ix(:,ipermute(M,[1,2,3,4,5,6]))...
        +Ix(:,ipermute(M,[1,3,4,6,2,5]))+Ix(:,ipermute(M,[1,2,4,6,3,5]))+Ix(:,ipermute(M,[1,2,3,6,4,5]))...
        +Ix(:,ipermute(M,[1,2,5,6,3,4]))+Ix(:,ipermute(M,[1,3,5,6,2,4]))+Ix(:,ipermute(M,[1,4,5,6,2,3]));
end

end

