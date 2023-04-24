function [ GAMMA ] = create_GAMMA( n_x,n_theta,varargin )
%The function creates sparse matrices GAMMA1, GAMMA2, ... 
%
% © Copyright, Oren Levintal, June 13, 2016.

if length(varargin)==1
    order_vec=varargin{1}; % do only specified order
else
    order_vec=2:3; % do all orders
end

for order=order_vec

    if order==2
        ind=1:n_x*n_x;
        M=reshape(ind,n_x,n_x);
        Ix=speye(n_x*n_x);
        GAMMA.GAMMA1=Ix(:,ipermute(M,[2,1]))+Ix(:,ipermute(M,[1,2]));
    end

    if order==3
        ind=1:n_x^3;
        M=reshape(ind,n_x,n_x,n_x);
        Ix=speye(n_x^3);
        GAMMA.GAMMA2=Ix(:,ipermute(M,[3,2,1]))+Ix(:,ipermute(M,[2,3,1]))+Ix(:,ipermute(M,[1,3,2]));

        GAMMA.GAMMA3=Ix(:,ipermute(M,[2,3,1]))+Ix(:,ipermute(M,[1,3,2]))+Ix(:,ipermute(M,[1,2,3]));
        
        GAMMA.GAMMA4=Ix(:,ipermute(M,[2,3,1]))+Ix(:,ipermute(M,[1,3,2]))+Ix(:,ipermute(M,[1,2,3]));    

        GAMMA.GAMMA5=Ix(:,ipermute(M,[ 2,3,1 ]))+Ix(:,ipermute(M,[ 1,3,2 ]))+Ix(:,ipermute(M,[ 1,2,3 ]));    
    end
end

end


