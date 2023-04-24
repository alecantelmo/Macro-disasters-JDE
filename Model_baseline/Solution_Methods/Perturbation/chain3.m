function [ fxxx ] = chain3(fv,fvv,fvvv,vx,vxx,vxxx,varargin)
%chain3 calculates the matrix of third derivatives of the composite
%function f(v(x)) with respect to x, using the third order multivariate chain
%rule derived in Levintal, Oren, "Fifth Order Perturbation Solution
%to DSGE Models", July 6, 2014.
%
%   Input arguments: 
%   fv is the Jacobian matrix of f(v) wrt v.
%   fvv is an array of the second derivatives of f wrt v.
%   fvvv is an array of the third derivatives of f wrt v.
%   vx is the Jacobian matrix of v(x) wrt x.
%   vxx is an array of the second derivatives of v wrt x.
%   vxxx is an array of the third derivatives of v wrt x.
%
%   All arguments, except fv and vx, can be reshaped in any form by the
%   reshape.m function. fv and vx must be in matrix form, where the rows of
%   fv correspond to the rows of f(v), and the row of vx correspond to the
%   rows of v(x).
%   
%   [ fxxx ] = chain3(fv,fvv,fvvv,vx,vxx,vxxx,OMEGA1) uses the user-defined
%   OMEGA1 matrix.
%
% © Copyright, Oren Levintal, June 13, 2016.


[n_f,n_v]=size(fv);
n_x=size(vx,2);

if isempty(varargin)==0
    OMEGA1=varargin{1};
end

% fvvv(vx^3)

if isempty(fvvv)
    term1=sparse(1,1);
else
    term1=innerkron(n_f,n_v,fvvv,vx,vx,vx);    
end

% fvv(vx*vxx)

term2=innerkron(n_f,n_v,fvv,vx,vxx);    

if isempty(varargin)
    term2=reshape(full(term2),n_f,n_x,n_x,n_x);
    term2_OMEGA1=ipermute(term2,[1,3,4,2])+ipermute(term2,[1,2,4,3])+ipermute(term2,[1,2,3,4]);
    term2_OMEGA1=reshape(term2_OMEGA1,n_f,n_x^3);
else
    term2_OMEGA1=term2*OMEGA1;
end

% fv(vxxx)

if isempty(vxxx)
    term3=sparse(1,1);
else
    term3=innerkron(n_f,n_v,fv,vxxx);
end

fxxx=term1+term2_OMEGA1+term3;


end

