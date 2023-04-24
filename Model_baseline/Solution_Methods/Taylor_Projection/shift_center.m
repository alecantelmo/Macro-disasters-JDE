function [ coeffs ] = shift_center( coeffs,c0,newc0,order,model )
%The function shifts the center of the approximating power series from c0
%to newc0.
%   coeffs: vector of Polynomial coefficients
%   c0: old center of power series
%   newc0: new center of power series
%   n_f: number of power series
%   n_b: number of terms in power series (excluding symmetric monomials)
%   n_x: number of variables
%   order: Polynomial order of power series
%   U: cell variable with compression matrices
%   W: cell variable with uncompression matrices
%
% © Copyright, Oren Levintal, June 13, 2016.

n_f=model.n_f;
n_b=model.n_b;
n_x=model.n_x;
U=model.U;
W=model.W;

coeffs=reshape(coeffs,n_f,n_b);
GH0=coeffs(:,1);
if order==1
    GH1=coeffs(:,2:1+n_x);
    [ GH0,GH1 ] = shiftpoly( newc0,c0,[],GH0,GH1 );
    coeffs=[GH0,GH1];
    coeffs=coeffs(:);
elseif order==2
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,2+n_x:1+n_x+model.unique2)*W{2};
    [ GH0,GH1,GH2 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2 );
    coeffs=[GH0,GH1,GH2*U{2}];
    coeffs=coeffs(:);
elseif order==3
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,2+n_x:1+n_x+model.unique2)*W{2};
    GH3=coeffs(:,2+n_x+model.unique2:1+n_x+model.unique2+model.unique3)*W{3};
    [ GH0,GH1,GH2,GH3 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2,GH3 );
    coeffs=[GH0,GH1,GH2*U{2},GH3*U{3}];
    coeffs=coeffs(:);
else 
    error('order must be 1, 2, or 3')
end

coeffs=full(coeffs);

end

