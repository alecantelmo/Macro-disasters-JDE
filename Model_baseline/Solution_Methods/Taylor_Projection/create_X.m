function [X,Xx,Xxx,Xxxx,Xxxxx]=create_X(N,n_x,x,c0,W2,unique2,W3,unique3)
% compute the basis function and its derivatives
%
% © Copyright, Oren Levintal, June 13, 2016.

Xxx=[]; Xxxx=[]; Xxxxx=[];
x_c0=sparse(x-c0);
X=[1;x_c0];
Xx=[sparse(1,n_x);speye(n_x)];
if N==1
    Xxx=sparse(1+n_x,n_x^2);
elseif N>=2
    tempx=reshape(reshape(W2,unique2*n_x,n_x)*x_c0,unique2,n_x);
    temp=tempx*x_c0;
    X=[X;temp];
    Xx=[Xx;2*tempx];
    Xxx=[sparse(1+n_x,n_x^2);2*W2];
end
if N==2
    Xxxx=sparse(1+n_x+unique2,n_x^3);
elseif N>=3
    tempxx=reshape(reshape(W3,unique3*n_x^2,n_x)*x_c0,unique3,n_x^2);
    tempx=reshape(reshape(tempxx,unique3*n_x,n_x)*x_c0,unique3,n_x);
    temp=tempx*x_c0;
    X=[X;temp];
    Xx=[Xx;3*tempx];
    Xxx=[Xxx;6*tempxx];
    Xxxx=[sparse(1+n_x+unique2,n_x^3);6*W3];
end
if N==3
    Xxxxx=sparse(1+n_x+unique2+unique3,n_x^4);
end
